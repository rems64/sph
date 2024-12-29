#include "ParticleSystem.hpp"
#include "Globals.hpp"
#include "Window.hpp"
#include "glm/common.hpp"
#include "glm/geometric.hpp"
#include "typedefs.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>

#include "imgui.h"

#define TIME(name, x) \
    { \
        float _start_time = get_time(); \
        x; \
        float _delta_time = get_time() - _start_time; \
        std::cout << "[time] " name << " " << 1000 * _delta_time << std::endl; \
    }

static real_t density_kernel(real_t dst, real_t radius) {
    if (dst < radius) {
        float scale = 15 / (2 * M_PI * pow(radius, 5));
        float v = radius - dst;
        return v * v * scale;
    }
    return 0;
}

static real_t density_derivative(real_t dst, real_t radius) {
    if (dst <= radius) {
        float scale = 15 / (pow(radius, 5) * M_PI);
        float v = radius - dst;
        return -v * scale;
    }
    return 0;
}

static real_t near_density_kernel(real_t dst, real_t radius) {
    if (dst < radius) {
        float scale = 15 / (M_PI * pow(radius, 6));
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
    // TODO: use actual prime numbers
    return (0x123456 * cell.x + 0x654232367 * cell.y + 0x734457 * cell.z) % m_particles_count;
}

void ParticleSystem::build_cells() {
    m_spatial_lookup.resize(m_particles_count);
    m_cell_start_index.resize(m_particles_count);
    for (index_t i = 0; i < m_particles_count; i++) {
        const auto position = m_predicted_positions[i];
        const cell_t cell = get_cell(position);
        const cellhash_t hash = hash_cell(cell);
        m_spatial_lookup[i] = std::make_pair(hash, i);
    }
    std::sort(m_spatial_lookup.begin(), m_spatial_lookup.end());
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
        const auto hash = hash_cell(cell);
        index_t index = m_cell_start_index[hash];
        while (index < m_particles_count) {
            const auto hash_index = m_spatial_lookup[index];
            if (hash_index.first != hash)
                break;
            particle_indices.push_back(hash_index.second);
            index++;
        }
    }
    return particle_indices;
}

ParticleSystem::ParticleSystem() : m_particles_count(500) {
    m_positions.resize(m_particles_count);
    m_velocities.resize(m_particles_count);
    m_densities.resize(m_particles_count);
    m_near_densities.resize(m_particles_count);
    m_predicted_positions.resize(m_particles_count);

    for (index_t i = 0; i < m_particles_count; i++) {
        m_positions[i] = m_extents * vec3(((float)rand() / RAND_MAX - .5f),
                                          ((float)rand() / RAND_MAX - .5f),
                                          ((float)rand() / RAND_MAX - .5f));
        m_velocities[i] = vec3(0);
    }
}

const vec3 ParticleSystem::calculate_pressure(index_t i) {
    const vec3 sample_point = m_predicted_positions[i];
    vec3 pressure_force = vec3(0);
    const real_t density = m_densities[i];
    const real_t near_density = m_near_densities[i];
    const real_t pressure = convert_density_to_pressure(density);
    const real_t near_pressure = convert_density_to_pressure(density);

    const auto neighbors = get_neighbors_particles(i);
    bool missed = false;
    for (const auto j : neighbors) {
        // for (index_t j = 0; j < m_particles_count; j++) {
        if (i == j)
            continue;
        real_t dst = glm::distance(sample_point, m_predicted_positions[j]);
        if (dst > (sqrt(3.) * m_smoothing_radius) && !missed && (missed |= true))
            G.debug.missed_cells++;
        vec3 dir;
        if (dst > 0.00001f)
            dir = (m_predicted_positions[j] - sample_point) / dst;
        else
            dir = vec3(1, 0, 0);
        const real_t neighbor_density = m_densities[j];
        const real_t near_neighbor_density = m_densities[j];
        const real_t shared_pressure = (pressure + convert_density_to_pressure(neighbor_density));
        const real_t shared_near_pressure = (near_pressure + convert_near_density_to_near_pressure(
                                                                 near_neighbor_density));
        pressure_force += dir * density_derivative(dst, m_smoothing_radius) * shared_pressure /
                          neighbor_density;
        pressure_force += dir * near_density_derivative(dst, m_smoothing_radius) *
                          shared_near_pressure / near_neighbor_density;
    }
    return pressure_force;
}

const real_t ParticleSystem::convert_density_to_pressure(real_t density) {
    return (density - m_target_density) * m_pressure_multiplier;
}

const real_t ParticleSystem::convert_near_density_to_near_pressure(real_t density) {
    return density * m_near_pressure_multiplier;
}

const void ParticleSystem::compute_densities() {
    // Compute densities
    for (index_t i = 0; i < m_particles_count; i++) {
        m_densities[i] = 0;
        const auto neighbors = get_neighbors_particles(i);
        // Calculate density
        // for (index_t j = 0; j < m_particles_count; j++) {
        for (const auto j : neighbors) {
            real_t dst = glm::distance(m_predicted_positions[i], m_predicted_positions[j]);
            m_densities[i] += density_kernel(dst, m_smoothing_radius);
            m_near_densities[i] += near_density_kernel(dst, m_smoothing_radius);
        }
    }
}

const void ParticleSystem::apply_gravity(real_t dt_scaled) {
    // Apply gravity
    m_densities.resize(m_particles_count);
    for (index_t i = 0; i < m_particles_count; i++) {
        m_velocities[i] += vec3(0, 0, 1) * m_gravity * dt_scaled;
        m_predicted_positions[i] = m_positions[i] + m_velocities[i] * dt_scaled;
    }
}

const void ParticleSystem::apply_pressure(real_t dt_scaled) {
    // Apply pressure
    for (index_t i = 0; i < m_particles_count; i++) {
        vec3 pressure_force = calculate_pressure(i);
        vec3 pressure_acceleration = pressure_force / m_densities[i];
        m_velocities[i] += pressure_acceleration * dt_scaled;
    }
}

void ParticleSystem::update(real_t dt) {
    G.debug.missed_cells = 0;
    const real_t t = G.t.time;
    const real_t dt_scaled = dt * m_simulation_speed;

    std::cout << std::endl;
    TIME("apply_gravity", apply_gravity(dt_scaled))
    TIME("build_cells", build_cells())
    TIME("compute_densities", compute_densities())
    TIME("apply_pressure", apply_pressure(dt_scaled))

    // Integrate and resolve collisions
    for (index_t i = 0; i < m_particles_count; i++) {
        m_positions[i] += m_velocities[i] * dt_scaled;

        resolve_collision(i);
    }
}

void ParticleSystem::resolve_collision(index_t i) {
    if (std::abs(m_positions[i].x) > m_extents.x / 2.f) {
        m_positions[i].x = glm::sign(m_positions[i].x) * m_extents.x / 2.f;
        m_velocities[i].x *= -1.f * m_collision_damping;
        m_velocities[i].y *= m_collision_tangent_damping;
        m_velocities[i].z *= m_collision_tangent_damping;
    }
    if (std::abs(m_positions[i].y) > m_extents.y / 2.f) {
        m_positions[i].y = glm::sign(m_positions[i].y) * m_extents.y / 2.f;
        m_velocities[i].y *= -1.f * m_collision_damping;
        m_velocities[i].x *= m_collision_tangent_damping;
        m_velocities[i].z *= m_collision_tangent_damping;
    }
    if (std::abs(m_positions[i].z) > m_extents.z / 2.f) {
        m_positions[i].z = glm::sign(m_positions[i].z) * m_extents.z / 2.f;
        m_velocities[i].z *= -1.f * m_collision_damping;
        m_velocities[i].x *= m_collision_tangent_damping;
        m_velocities[i].y *= m_collision_tangent_damping;
    }
}

const float ParticleSystem::max_velocity() { return m_max_velocity; }

void ParticleSystem::imgui_controls() {
    ImGui::Begin("Particle System Controls");
    ImGui::DragFloat("Simulation Speed", &m_simulation_speed, 0.1f, 0.1f, 10.f);
    ImGui::DragFloat("Extent X", &m_extents.x, 0.1f, 0.f, 10.f);
    ImGui::DragFloat("Extent Y", &m_extents.y, 0.1f, 0.f, 10.f);
    ImGui::DragFloat("Extent Z", &m_extents.z, 0.1f, 0.f, 10.f);

    ImGui::DragFloat("Max Velocity", &m_max_velocity, 0.1f, 0.1f, 100.f);

    ImGui::DragFloat("Collision Damping", &m_collision_damping, 0.1f, 0.1f, 1.f);
    ImGui::DragFloat("Collision Tangent Damping", &m_collision_tangent_damping, 0.1f, 0.1f, 1.f);

    ImGui::DragFloat("Gravity", &m_gravity, 0.1f, -20.f, 20.f);

    ImGui::DragFloat("Smoothing Radius", &m_smoothing_radius, 0.01f, 1.f);
    ImGui::DragFloat("Target Density", &m_target_density, 0.001f, 10.f);
    ImGui::DragFloat("Pressure Multiplier", &m_pressure_multiplier, 0.1f, 0.001f, 10.f);
    ImGui::DragFloat("Near Pressure Multiplier", &m_near_pressure_multiplier, 0.1f, 0.001f, 10.f);

    ImGui::Text("Missed cells : %u (%.1f%%)", G.debug.missed_cells, (float)G.debug.missed_cells/m_particles_count);

    ImGui::End();
}

const std::vector<vec3> &ParticleSystem::positions() const { return m_positions; }
const std::vector<vec3> &ParticleSystem::velocities() const { return m_velocities; }
const std::vector<real_t> &ParticleSystem::densities() const { return m_densities; }
vec3 &ParticleSystem::extent_ref() { return m_extents; }
