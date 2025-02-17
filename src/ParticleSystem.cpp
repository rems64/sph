#include "ParticleSystem.hpp"
#include "Globals.hpp"
#include "Kernels.hpp"
#include "ResourceManager.hpp"
#include "Transform.hpp"
#include "Window.hpp"
#include "glm/common.hpp"
#include "glm/ext/quaternion_transform.hpp"
#include "glm/geometric.hpp"
#include "glm/gtc/quaternion.hpp"
#include "mem.hpp"
#include "typedefs.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <thread>
#include <utility>
#include <vector>

#include "imgui.h"

#include <threads.h>

#define TIME(name, x) \
    { \
        float _start_time = get_time(); \
        x; \
        float _delta_time = get_time() - _start_time; \
    }
// std::cout << "[time] " name << " " << 1000 * _delta_time << std::endl; \

static real_t density_kernel(real_t dst, real_t radius) { return base_kernel(dst, radius); }

static real_t density_derivative(real_t dst, real_t radius) {
    return base_kernel_derivative(dst, radius);
}

static real_t near_density_kernel(real_t dst, real_t radius) {
    if (dst < radius) {
        // float scale = 15 / (M_PI * pow(radius, 6));
        float scale = 2.f / (M_PI * M_PI * powf(radius, 4));
        float v = radius - dst;
        return v * v * v * scale;
    }
    return 0;
}

static real_t near_density_derivative(real_t dst, real_t radius) {
    if (dst <= radius) {
        float scale = 45 / (pow(radius, 6) * M_PI);
        float v = radius - dst;
        return -v * v * scale;
    }
    return 0;
}

const cellhash_t ParticleSystem::hash_cell(cell_t cell) const {
    const auto hash = std::hash<cell_t>()(cell);
    return hash % m_particles_count;
    // return (std::hash<size_t>()(cell.x) ^ std::hash<size_t>()(cell.y) ^
    // std::hash<size_t>()(cell.z)) % m_particles_count;
}

const void ParticleSystem::clear_updated() {
    // const size_t updated_particles_size = std::ceil((float)m_particles_count / 32.f);
    // m_updated_particles.resize(updated_particles_size);
    // std::fill_n(m_updated_particles.begin(), updated_particles_size, 0);
    // m_updated_particles.resize(m_particles_count);
    for (index_t i = 0; i < m_particles_count; i++) {
        if (m_updated_particles.size() < i)
            m_updated_particles.emplace_back(false);
        m_updated_particles[i].store(false);
    }
}

static const vec3 random_dir() {
    return glm::normalize(vec3(drand48() - 0.5f, drand48() - 0.5f, drand48() - 0.5f));
}

const void ParticleSystem::set_updated(index_t index) {
    // const auto info = div(index, 32);
    // m_updated_particles[info.quot] |= (1u << info.rem);

    // To prevent unecessary substraction, we store each integer reversed
    m_updated_particles[index].store(true);
}

const bool ParticleSystem::is_updated(index_t index) {
    // const auto info = div(index, 32);
    // return m_updated_particles[info.quot] & (1u << info.rem);
    return m_updated_particles[index].load();
}

const bool ParticleSystem::is_wall(index_t index) { return index < m_first_simulated_particle; }

void ParticleSystem::build_cells() {
    m_active_cells.clear();
    m_spatial_lookup.resize(m_particles_count);
    m_cell_start_index.resize(m_particles_count);
    for (index_t i = 0; i < m_particles_count; i++) {
        const auto position = m_predicted_positions[i];
        const cell_t cell = get_cell(position);
        m_active_cells.insert(cell);
        const cellhash_t hash = hash_cell(cell);
        m_spatial_lookup[i] = std::make_pair(hash, i);
    }
    std::stable_sort(m_spatial_lookup.begin(), m_spatial_lookup.end());
    index_t last_hash = 0;
    for (index_t i = 0; i < m_particles_count; i++) {
        const auto hash_index = m_spatial_lookup[i];
        if (i == 0 || hash_index.first != last_hash)
            m_cell_start_index[hash_index.first] = i;
        last_hash = hash_index.first;
    }
}

cell_t ParticleSystem::get_cell(const vec3 position) const {
    const auto grid_size = m_smoothing_radius;
    const auto x = std::floor(position.x / grid_size);
    const auto y = std::floor(position.y / grid_size);
    const auto z = std::floor(position.z / grid_size);
    return cell_t(x, y, z);
}

const std::vector<index_t> ParticleSystem::get_particles_in_cell(cell_t cell) const {
    const auto hash = hash_cell(cell);
    return get_particles_in_cell_from_hash(hash);
}

const std::vector<index_t> ParticleSystem::get_particles_in_cell_from_hash(cellhash_t hash) const {
    std::vector<index_t> particle_indices = {};
    index_t index = m_cell_start_index[hash];
    while (index < m_particles_count) {
        const auto hash_index = m_spatial_lookup[index];
        if (hash_index.first != hash)
            break;
        particle_indices.push_back(hash_index.second);
        index++;
    }
    return particle_indices;
}

std::vector<index_t> ParticleSystem::get_neighbors_particles(const index_t particle_index) const {
    std::vector<index_t> particle_indices = {};
    const auto position = m_predicted_positions[particle_index];
    const cell_t cell = get_cell(position);
    const cell_t neighbor_cells[] = {
        cell_t(cell.x + 1, cell.y + 1, cell.z + 1), cell_t(cell.x + 0, cell.y + 1, cell.z + 1),
        cell_t(cell.x - 1, cell.y + 1, cell.z + 1), cell_t(cell.x + 1, cell.y + 0, cell.z + 1),
        cell_t(cell.x + 0, cell.y + 0, cell.z + 1), cell_t(cell.x - 1, cell.y + 0, cell.z + 1),
        cell_t(cell.x + 1, cell.y - 1, cell.z + 1), cell_t(cell.x + 0, cell.y - 1, cell.z + 1),
        cell_t(cell.x - 1, cell.y - 1, cell.z + 1), cell_t(cell.x + 1, cell.y + 1, cell.z + 0),
        cell_t(cell.x + 0, cell.y + 1, cell.z + 0), cell_t(cell.x - 1, cell.y + 1, cell.z + 0),
        cell_t(cell.x + 1, cell.y + 0, cell.z + 0), cell_t(cell.x + 0, cell.y + 0, cell.z + 0),
        cell_t(cell.x - 1, cell.y + 0, cell.z + 0), cell_t(cell.x + 1, cell.y - 1, cell.z + 0),
        cell_t(cell.x + 0, cell.y - 1, cell.z + 0), cell_t(cell.x - 1, cell.y - 1, cell.z + 0),
        cell_t(cell.x + 1, cell.y + 1, cell.z - 1), cell_t(cell.x + 0, cell.y + 1, cell.z - 1),
        cell_t(cell.x - 1, cell.y + 1, cell.z - 1), cell_t(cell.x + 1, cell.y + 0, cell.z - 1),
        cell_t(cell.x + 0, cell.y + 0, cell.z - 1), cell_t(cell.x - 1, cell.y + 0, cell.z - 1),
        cell_t(cell.x + 1, cell.y - 1, cell.z - 1), cell_t(cell.x + 0, cell.y - 1, cell.z - 1),
        cell_t(cell.x - 1, cell.y - 1, cell.z - 1),
    };
    for (const auto cell : neighbor_cells) {
        const auto cell_particles = get_particles_in_cell(cell);
        for (const auto neighbor_particle_index : cell_particles) {
            particle_indices.push_back(neighbor_particle_index);
        }
    }
    return particle_indices;
}

void ParticleSystem::spawn_shell(vec3 extent) {
    const real_t wall_spacing = 0.1f;
    const vec3 half_extent = extent / 2;

    for (real_t x = -half_extent.x; x <= half_extent.x; x += wall_spacing) {
        for (real_t y = -half_extent.y; y <= half_extent.y; y += wall_spacing) {
            m_positions.push_back(vec3(x, y, -half_extent.z));
            m_positions.push_back(vec3(x, y, half_extent.z));
        }
    }

    for (real_t x = -half_extent.x; x <= half_extent.x; x += wall_spacing) {
        for (real_t z = -half_extent.z; z <= half_extent.z; z += wall_spacing) {
            m_positions.push_back(vec3(x, -half_extent.z, z));
            m_positions.push_back(vec3(x, half_extent.z, z));
        }
    }

    for (real_t y = -half_extent.y; y <= half_extent.y; y += wall_spacing) {
        for (real_t z = -half_extent.z; z <= half_extent.z; z += wall_spacing) {
            m_positions.push_back(vec3(-half_extent.z, y, z));
            m_positions.push_back(vec3(half_extent.z, y, z));
        }
    }
}

void ParticleSystem::spawn_walls() {
    const size_t layers_count = 0;
    const real_t layer_width = 0.05f;
    for (size_t layer = 0; layer < layers_count; layer++) {
        spawn_shell(m_extents - vec3(layer_width) * layer);
    }

    m_first_simulated_particle = m_positions.size();
}

ParticleSystem::ParticleSystem(handle<Transform> transform)
    : m_transform(transform), m_simulated_particles_count(500), m_hist_active_cells_count(100),
      m_hist_updated_particles_count(100), m_first_simulated_particle(0) {

    // spawn_walls();

    m_particles_count = m_positions.size() + m_simulated_particles_count;

    m_positions.resize(m_particles_count);
    m_velocities.resize(m_particles_count);
    m_densities.resize(m_particles_count);
    m_near_densities.resize(m_particles_count);
    m_predicted_positions.resize(m_particles_count);
    m_is_neighbor.resize(m_particles_count);

    // IISPH
    m_rho.resize(m_particles_count);
    m_v.resize(m_particles_count);
    m_v_advection.resize(m_particles_count);
    m_rho_advection.resize(m_particles_count);
    m_dii.resize(m_particles_count);
    m_p_current.resize(m_particles_count);
    m_p_new.resize(m_particles_count);
    m_aii.resize(m_particles_count);
    m_sums.resize(m_particles_count);

    for (index_t i = m_first_simulated_particle; i < m_particles_count; i++) {
        m_positions[i] = (m_extents - 0.05f * vec3(1)) *
                         vec3(drand48() - .5f, drand48() - .5f, drand48() - .5f);
        m_velocities[i] = vec3(0);
    }
}

const vec3 ParticleSystem::calculate_pressure(index_t i) {
    const vec3 sample_point = m_predicted_positions[i];
    vec3 pressure_force = vec3(0);
    const real_t density = m_densities[i];
    // const real_t near_density = m_near_densities[i];
    const real_t pressure = equation_of_state(density);
    // const real_t near_pressure = equation_of_state(density);
    const auto neighbor_particles = get_neighbors_particles(i);
    if (i == m_tracked_particle)
        m_is_neighbor[i] = 1.f;
    for (const auto j : neighbor_particles) {
        // for (index_t j = 0; j < m_particles_count; j++) {
        if (i == j)
            continue;
        real_t dst = glm::distance(sample_point, m_predicted_positions[j]);

        if (dst > m_smoothing_radius)
            continue;
        vec3 dir;
        if (dst < 0.00001f) {
            m_predicted_positions[j] += 0.0001f * random_dir();
            dst += 0.00001f;
        } else {
            dir = (m_predicted_positions[j] - sample_point) / dst;
        }
        if (i == m_tracked_particle)
            m_is_neighbor[j] = 0.5f;
        // const real_t neighbor_density = m_densities[j];
        // const real_t near_neighbor_density = m_densities[j];
        // const real_t shared_near_pressure =
        //     0.5f * (near_pressure +
        //     convert_near_density_to_near_pressure(near_neighbor_density));
        // if (neighbor_density > 0.001f)
        //     pressure_force += dir * density_derivative(dst, m_smoothing_radius) *
        //     shared_pressure /
        //                       neighbor_density;
        // if (near_neighbor_density > 0.001f)
        //     pressure_force += dir * near_density_derivative(dst, m_smoothing_radius) *
        //                       shared_near_pressure / near_neighbor_density;
        const real_t neighbor_density = m_densities[j];
        const real_t shared_pressure = 0.5f * (pressure + equation_of_state(neighbor_density));
        const real_t pressure = shared_pressure / neighbor_density *
                                density_derivative(dst, m_smoothing_radius);
        if (std::abs(neighbor_density) > 0.001f)
            pressure_force += pressure * dir;
    }
    return pressure_force;
}

const vec3 ParticleSystem::calculate_viscosity(index_t i) {
    const vec3 sample_point = m_predicted_positions[i];
    const vec3 velocity = m_velocities[i];
    vec3 viscosity_force = vec3(0);
    const real_t density = m_densities[i];

    for (const auto j : get_neighbors_particles(i)) {
        if (i == j)
            continue;
        const vec3 neighbor_position = m_predicted_positions[j];
        real_t dst = glm::distance(sample_point, neighbor_position);
        vec3 dir;
        if (dst < 0.00001f) {
            m_predicted_positions[j] += 0.0001f * random_dir();
            dst += 0.00001f;
        } else {
            dir = (m_predicted_positions[j] - sample_point) / dst;
        }
        viscosity_force += (m_velocities[j] - velocity) * density_kernel(dst, m_smoothing_radius);
    }
    return viscosity_force;
}

const real_t ParticleSystem::equation_of_state(real_t density) {
    return std::max<real_t>(0, (density - m_target_density) * m_pressure_multiplier);
    // return std::max<real_t>(
    //     0, (std::pow(density / m_target_density, 7.f) - 1) * m_pressure_multiplier);
}

const real_t ParticleSystem::convert_near_density_to_near_pressure(real_t density) {
    return std::max<real_t>(0, (density)*m_near_pressure_multiplier);
}

const void ParticleSystem::compute_densities(cellhash_t hash) {
    // Compute densities
    if (m_densities.size() != m_particles_count)
        m_densities.resize(m_particles_count);
    const auto particles_indices = get_particles_in_cell_from_hash(hash);
    for (const auto i : particles_indices) {
        if (is_updated(i))
            continue;
        set_updated(i);
        // for (index_t i = 0; i < m_particles_count; i++) {
        m_densities[i] = 0;
        const auto neighbors = get_neighbors_particles(i);
        for (const auto j : neighbors) {
            real_t dst = glm::distance(m_predicted_positions[i], m_predicted_positions[j]);
            m_densities[i] += density_kernel(dst, m_smoothing_radius);
            m_near_densities[i] += near_density_kernel(dst, m_smoothing_radius);
        }
    }
}

void ParticleSystem::set_container_speed(const vec3 &speed) {
    const auto dt = G.t.dt;
    const auto speed_local = glm::inverse(m_transform.lock()->rotation()) * speed;
    for (index_t i = 0; i < m_particles_count; i++) {
        m_positions[i] -= speed_local * dt;
    }
}

void ParticleSystem::set_container_delta_angle(const quat &delta_angle) {
    const auto delta_angle_inv = glm::inverse(m_transform.lock()->rotation()) *
                                 glm::inverse(delta_angle) * m_transform.lock()->rotation();
    for (index_t i = 0; i < m_particles_count; i++) {
        m_positions[i] = delta_angle_inv * m_positions[i];
    }
};

const void ParticleSystem::apply_gravity(cellhash_t hash) {
    vec3 gravity = glm::inverse(m_transform.lock()->rotation()) * vec3(0, 0, m_gravity);
    const auto particles_indices = get_particles_in_cell_from_hash(hash);
    for (const auto i : particles_indices) {
        if (is_updated(i) || is_wall(i))
            continue;
        set_updated(i);
        m_velocities[i] += gravity * dt_scaled;
    }
}

const void ParticleSystem::compute_predicted_position() {
    vec3 gravity = glm::inverse(m_transform.lock()->rotation()) * vec3(0, 0, m_gravity);
    // const auto particles_indices = get_particles_in_cell_from_hash(hash);
    for (index_t i = 0; i < m_particles_count; i++) {
        // for (const auto i : particles_indices) {
        // if (is_updated(i))
        //     continue;
        // set_updated(i);
        m_predicted_positions[i] = m_positions[i] + m_velocities[i] * dt_scaled;
        m_is_neighbor[i] = 0.f;
    }
}

const void ParticleSystem::apply_pressure(cellhash_t hash) {
    // Apply pressure
    // for (index_t i = 0; i < m_particles_count; i++) {
    const auto particles_indices = get_particles_in_cell_from_hash(hash);
    for (const auto i : particles_indices) {
        if (is_updated(i) || is_wall(i))
            continue;
        set_updated(i);
        vec3 pressure_force = calculate_pressure(i);
        vec3 viscosity_force = calculate_viscosity(i);
        vec3 pressure_acceleration = pressure_force / m_densities[i];
        vec3 viscosity_acceleration = viscosity_force * m_viscosity_strength;
        m_velocities[i] += (pressure_acceleration + viscosity_acceleration) * dt_scaled;
    }
}

const void ParticleSystem::integrate_and_collide(cellhash_t hash) {
    // Integrate and resolve collisions
    const auto colliders = G.resource_manager.lock()->colliders();
    // for (index_t i = 0; i < m_particles_count; i++) {
    const auto particles_indices = get_particles_in_cell_from_hash(hash);
    for (const auto i : particles_indices) {
        if (is_updated(i) || is_wall(i))
            continue;
        set_updated(i);
        G.debug.particles_updated++;
        vec3 previous = m_positions[i];
        m_velocities[i] = glm::clamp(m_velocities[i], -m_max_velocity, m_max_velocity);
        m_positions[i] += m_velocities[i] * dt_scaled;

        resolve_collision(i);

        for (const auto &collider : colliders) {
            const vec3 world_position = m_transform.lock()->world_matrix() *
                                        vec4(m_positions[i], 1);
            const auto collision_result = collider->collide(world_position, previous);
            if (collision_result.collided) {
                const vec3 mtv = collision_result.mtv;
                const vec3 mtv_local = glm::inverse(m_transform.lock()->rotation()) * mtv;
                m_positions[i] += mtv_local;
                m_velocities[i] = m_collision_damping * glm::reflect(m_velocities[i], mtv_local);
            }
        }

        resolve_collision(i);
    }
}

void ParticleSystem::join_threads() {
    for (auto &thread : m_threads) {
        thread.join();
    }
}

// IISPH procedures
const void ParticleSystem::predict_advection(cellhash_t hash) {
    vec3 gravity = glm::inverse(m_transform.lock()->rotation()) * vec3(0, 0, m_gravity);
    const auto particles_indices = get_particles_in_cell_from_hash(hash);

    for (const auto i : particles_indices) {
        m_rho[i] = 0;
        const auto neighbors = get_neighbors_particles(i);
        for (const auto j : neighbors) {
            if (i == j)
                continue;
            m_rho[i] += density_kernel(
                glm::distance(m_predicted_positions[i], m_predicted_positions[j]),
                m_smoothing_radius);
        }
        m_v_advection[i] = m_v[i] + gravity * dt_scaled;
        m_dii[i] = vec3(0.);
        for (const auto j : neighbors) {
            if (i == j)
                continue;
            real_t dst = glm::distance(m_predicted_positions[i], m_predicted_positions[j]);

            if (dst > m_smoothing_radius)
                continue;
            vec3 dir;
            if (dst < 0.00001f) {
                m_predicted_positions[j] += 0.0001f * random_dir();
                dst += 0.00001f;
            }
            dir = (m_predicted_positions[j] - m_predicted_positions[i]) / dst;

            m_dii[i] += -dir * density_derivative(dst, m_smoothing_radius);
        }
        m_dii[i] *= dt_scaled * dt_scaled;
    }

    for (const auto i : particles_indices) {
        m_rho_advection[i] = m_rho[i];
        const auto neighbors = get_neighbors_particles(i);
        for (const auto j : neighbors) {
            const vec3 relative_speed = m_v_advection[i] - m_v_advection[j];

            if (i == j)
                continue;
            real_t dst = glm::distance(m_predicted_positions[i], m_predicted_positions[j]);

            if (dst > m_smoothing_radius)
                continue;
            vec3 dir;
            if (dst < 0.00001f) {
                m_predicted_positions[j] += 0.0001f * random_dir();
                dst += 0.00001f;
            }
            dir = (m_predicted_positions[j] - m_predicted_positions[i]) / dst;

            m_rho_advection[i] += dt_scaled * glm::dot(relative_speed, dir) *
                                  density_derivative(dst, m_smoothing_radius);
        }
        m_p_current[i] = 0.5f * m_p_current[i];
        m_aii[i] = 0.;
        for (const auto j : neighbors) {
            if (i == j)
                continue;
            real_t dst = glm::distance(m_predicted_positions[i], m_predicted_positions[j]);

            if (dst > m_smoothing_radius)
                continue;
            vec3 dir;
            if (dst < 0.00001f) {
                m_predicted_positions[j] += 0.0001f * random_dir();
                dst += 0.00001f;
            }
            dir = (m_predicted_positions[j] - m_predicted_positions[i]) / dst;

            if (m_rho[i] > 0.001f) {
                const vec3 d_ji = -(dt_scaled * dt_scaled) / (m_rho[i] * m_rho[i]) *
                                  density_derivative(dst, m_smoothing_radius) * dir;

                m_aii[i] += glm::dot(m_dii[i] - d_ji,
                                     density_derivative(dst, m_smoothing_radius) * dir);
            }
        }
    }
}
const void ParticleSystem::pressure_solve(cellhash_t hash) {
    const auto particles_indices = get_particles_in_cell_from_hash(hash);

    // Relaxation factor
    const real_t omega = 0.5f;

    for (index_t l = 0; l < 10; l++) {
        for (const auto i : particles_indices) {
            m_sums[i] = vec3(0.);
            const auto neighbors = get_neighbors_particles(i);
            for (const auto j : neighbors) {
                if (i == j)
                    continue;
                real_t dst = glm::distance(m_predicted_positions[i], m_predicted_positions[j]);

                if (dst > m_smoothing_radius)
                    continue;
                vec3 dir;
                if (dst < 0.00001f) {
                    m_predicted_positions[j] += 0.0001f * random_dir();
                    dst += 0.00001f;
                }
                dir = (m_predicted_positions[j] - m_predicted_positions[i]) / dst;

                if (m_rho[i] > 0.001f)
                    m_sums[i] += -1. / (m_rho[j] * m_rho[j]) * m_p_current[j] *
                                 density_derivative(dst, m_smoothing_radius) * dir;
            }
            m_sums[i] *= dt_scaled * dt_scaled;
        }

        for (const auto i : particles_indices) {
            m_p_new[i] = (1 - omega) * m_p_current[i];
            real_t sum = 0.;

            const auto neighbors = get_neighbors_particles(i);
            for (const auto j : neighbors) {
                if (i == j)
                    continue;
                real_t dst = glm::distance(m_predicted_positions[i], m_predicted_positions[j]);

                if (dst > m_smoothing_radius)
                    continue;
                vec3 dir;
                if (dst < 0.00001f) {
                    m_predicted_positions[j] += 0.0001f * random_dir();
                    dst += 0.00001f;
                }
                dir = (m_predicted_positions[j] - m_predicted_positions[i]) / dst;

                if (m_rho[i] > 0.001f) {
                    const vec3 d_ji = -(dt_scaled * dt_scaled) / (m_rho[i] * m_rho[i]) *
                                      density_derivative(dst, m_smoothing_radius) * dir;
                    vec3 other_sum = m_sums[j] - d_ji * m_p_current[i];

                    const vec3 gradient = density_derivative(dst, m_smoothing_radius) * dir;
                    sum += glm::dot(m_sums[i] - m_dii[j] * m_p_current[j] - other_sum, gradient);
                }

                // const auto neighbors_j = get_neighbors_particles(j);
                // for (const auto k : neighbors) {
                //     if (i == k)
                //         continue;

                //     real_t dst2 = glm::distance(m_predicted_positions[j],
                //                                 m_predicted_positions[k]);

                //     if (dst2 > m_smoothing_radius)
                //         continue;
                //     vec3 dir2;
                //     if (dst2 < 0.00001f) {
                //         m_predicted_positions[k] += 0.0001f * random_dir();
                //     } else {
                //         dir2 = (m_predicted_positions[k] - m_predicted_positions[j]) / dst2;
                //     }

                //     const vec3 d_jk = -(dt_scaled * dt_scaled) / (m_rho[k] * m_rho[k]) *
                //                       density_derivative(dst2, m_smoothing_radius) * dir2;

                //     other_sum += d_jk * m_p_current[k];
                // }
            }
            real_t aii = 1.f;
            if (m_aii[i] > 1.f)
                aii = m_aii[i];
            m_p_new[i] += omega / aii * (m_target_density - m_rho_advection[i] - sum);
        }

        for (const auto i : particles_indices) {
            m_p_current[i] = m_p_new[i];
        }
    }
}

const void ParticleSystem::integration(cellhash_t hash) {
    const auto particles_indices = get_particles_in_cell_from_hash(hash);

    for (const auto i : particles_indices) {
        const vec3 F_ip = 1. / (dt_scaled * dt_scaled) * (m_dii[i] * m_p_current[i] + m_sums[i]);
        m_v[i] = m_v_advection[i] + dt_scaled * F_ip;
        m_velocities[i] = m_v[i];
        m_densities[i] = m_rho[i];
        // m_positions[i] = m_positions[i] + dt_scaled * m_v[i];
    }
}

static bool isnan(vec3 &vec) {
    return std::isnan(vec.x) || std::isnan(vec.y) || std::isnan(vec.z);
}

void ParticleSystem::update(real_t dt) {
    m_tracked_particle = G.simulation.highlight;
    G.debug.missed_cells = 0;
    const real_t t = G.t.time;
    dt_scaled = std::min(dt * m_simulation_speed, 1.f * 60.f);
    // dt_scaled = 0.01f;

    for (size_t particle_index; particle_index < m_particles_count; particle_index++) {
        if (isnan(m_positions[particle_index]))
            m_positions[particle_index] = vec3(0.f);
        if (isnan(m_velocities[particle_index]))
            m_velocities[particle_index] = vec3(0.f);
    }

    compute_predicted_position();
    // std::cout << std::endl;
    build_cells();
    std::vector<cellhash_t> hashes;

    index_t last_hash = 0;
    for (index_t i = 0; i < m_particles_count; i++) {
        const auto hash_index = m_spatial_lookup[i];
        if (i == 0 || hash_index.first != last_hash)
            hashes.push_back(hash_index.first);
        last_hash = hash_index.first;
    }

    clear_updated();
    G.debug.particles_updated = 0;

    const auto active_cells_count = hashes.size();
    m_threads.resize(active_cells_count);

    // for (const auto cell : m_active_cells) {
    //     apply_gravity(cell);
    // }
    std::cout << "predict advection" << std::endl;

    size_t thread_index = 0;
    for (auto hash : hashes) {
        m_threads[thread_index++] = std::thread([this, hash]() { this->predict_advection(hash); });
    }
    join_threads();
    clear_updated();

    std::cout << "pressure solve" << std::endl;
    thread_index = 0;
    for (auto hash : hashes) {
        m_threads[thread_index++] = std::thread([this, hash]() { this->pressure_solve(hash); });
    }
    join_threads();
    clear_updated();

    std::cout << "integration" << std::endl;
    thread_index = 0;
    for (auto hash : hashes) {
        m_threads[thread_index++] = std::thread([this, hash]() { this->integration(hash); });
    }
    join_threads();
    clear_updated();

    // std::cout << "launched " << thread_index << " threads" << std::endl;
    // std::cout << "should have launched " << active_cells_count << std::endl;

    // for (const auto cell : m_active_cells) {
    //     compute_predicted_position(cell);
    // }
    // join_threads();
    // clear_updated();
    // thread_index = 0;
    // for (auto hash : hashes) {
    //     m_threads[thread_index++] = std::thread(
    //         [this, hash]() { this->compute_predicted_position(hash); });
    // }

    // for (const auto cell : m_active_cells) {
    //     compute_densities(cell);
    // }
    // join_threads();
    // clear_updated();
    // thread_index = 0;
    // for (auto hash : hashes) {
    //     m_threads[thread_index++] = std::thread([this, hash]() {
    //     this->compute_densities(hash);
    //     });
    // }

    // for (const auto cell : m_active_cells) {
    //     apply_pressure(cell);
    // }
    // join_threads();
    // clear_updated();
    // thread_index = 0;
    // for (auto hash : hashes) {
    //     m_threads[thread_index++] = std::thread([this, hash]() { this->apply_pressure(hash);
    //     });
    // }

    // for (const auto cell : m_active_cells) {
    //     integrate_and_collide(cell);
    // }
    // join_threads();
    // clear_updated();
    thread_index = 0;
    for (auto hash : hashes) {
        m_threads[thread_index++] = std::thread(
            [this, hash]() { this->integrate_and_collide(hash); });
    }

    join_threads();

    m_hist_active_cells_count.push_back(m_active_cells.size());
    m_hist_updated_particles_count.push_back(G.debug.particles_updated);
}

bool ParticleSystem::resolve_collision(index_t i) {
    bool collide = false;
    const auto extent = m_extents / 2 - vec3(1) * 0.05f;
    if (std::abs(m_positions[i].x) > extent.x) {
        m_positions[i].x = glm::sign(m_positions[i].x) * extent.x;
        m_velocities[i].x *= -1.f * m_collision_damping;
        m_velocities[i].y *= m_collision_tangent_damping;
        m_velocities[i].z *= m_collision_tangent_damping;
        collide = true;
    }
    if (std::abs(m_positions[i].y) > extent.y) {
        m_positions[i].y = glm::sign(m_positions[i].y) * extent.y;
        m_velocities[i].y *= -1.f * m_collision_damping;
        m_velocities[i].x *= m_collision_tangent_damping;
        m_velocities[i].z *= m_collision_tangent_damping;
        collide = true;
    }
    if (std::abs(m_positions[i].z) > extent.z) {
        m_positions[i].z = glm::sign(m_positions[i].z) * extent.z;
        m_velocities[i].z *= -1.f * m_collision_damping;
        m_velocities[i].x *= m_collision_tangent_damping;
        m_velocities[i].y *= m_collision_tangent_damping;
        collide = true;
    }
    return collide;
}

const float ParticleSystem::max_velocity() { return m_max_velocity; }

void ParticleSystem::imgui_controls() {
    ImGui::Begin("Particle System Controls");
    ImGui::DragFloat("Simulation Speed", &m_simulation_speed, 0.1f, 0.1f, 10.f);
    int count = m_particles_count;
    ImGui::DragInt("Particles count", &count);
    m_particles_count = count;
    ImGui::DragFloat("Extent X", &m_extents.x, 0.1f, 0.f, 10.f);
    ImGui::DragFloat("Extent Y", &m_extents.y, 0.1f, 0.f, 10.f);
    ImGui::DragFloat("Extent Z", &m_extents.z, 0.1f, 0.f, 10.f);

    ImGui::DragFloat("Max Velocity", &m_max_velocity, 0.1f, 0.1f, 100.f);

    ImGui::DragFloat("Collision Damping", &m_collision_damping, 0.1f, 0.1f, 1.f);
    ImGui::DragFloat("Collision Tangent Damping", &m_collision_tangent_damping, 0.1f, 0.1f, 1.f);

    ImGui::DragFloat("Gravity", &m_gravity, 0.1f, -20.f, 20.f);

    ImGui::DragFloat("Smoothing Radius", &m_smoothing_radius, 0.01f, 1.f);
    ImGui::DragFloat("Target Density", &m_target_density, 1.f, 0.f);
    ImGui::DragFloat("Pressure Multiplier", &m_pressure_multiplier, 0.1f, 0.f, 10.f);
    ImGui::DragFloat("Near Pressure Multiplier", &m_near_pressure_multiplier, 0.001f, 0.f, 1.f);
    ImGui::DragFloat("Viscosity strength", &m_viscosity_strength, 0.1f, 0.f, 10.f);

    ImGui::DragFloat("Density Error Offset", &G.debug.density_error_offset);
    ImGui::DragFloat("Density Color Range", &G.debug.density_color_range);

    ImGui::Separator();

    ImGui::Text("Particles count : %lu", m_particles_count);
    ImGui::Text("Missed cells : %u (%.1f%%)",
                G.debug.missed_cells,
                100.f * (float)G.debug.missed_cells / m_particles_count);

    ImGui::Text("Particles updated : %u", G.debug.particles_updated);
    ImGui::Text("Active cells count : %lu", m_active_cells.size());
    // ImGui::PlotHistogram("Active cells count", m_hist_active_cells_count.data(),
    // m_hist_active_cells_count.size());
    ImGui::PlotLines("Active cells count",
                     m_hist_active_cells_count.ptr(),
                     m_hist_active_cells_count.count(),
                     m_hist_active_cells_count.offset());

    ImGui::PlotLines("Updated particle count",
                     m_hist_updated_particles_count.ptr(),
                     m_hist_updated_particles_count.count(),
                     m_hist_updated_particles_count.offset());

    ImGui::End();
}

const std::vector<vec3> &ParticleSystem::positions() const { return m_positions; }
const std::vector<vec3> &ParticleSystem::velocities() const { return m_velocities; }
const std::vector<real_t> &ParticleSystem::densities() const { return m_densities; }
vec3 &ParticleSystem::extent_ref() { return m_extents; }
handle<Transform> ParticleSystem::transform() const { return m_transform; }
