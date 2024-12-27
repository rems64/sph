#include "ParticleSystem.hpp"
#include "Globals.hpp"
#include "glm/common.hpp"
#include "glm/geometric.hpp"
#include "typedefs.hpp"
#include <cmath>
#include <cstdlib>
#include <vector>

#include "imgui.h"

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

ParticleSystem::ParticleSystem() : m_particles_count(200) {
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
    for (index_t j = 0; j < m_particles_count; j++) {
        if (i == j)
            continue;
        real_t dst = glm::distance(sample_point, m_predicted_positions[j]);
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
        pressure_force += dir * density_derivative(dst, smoothing_radius) * shared_pressure /
                          neighbor_density;
        pressure_force += dir * near_density_derivative(dst, smoothing_radius) *
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

void ParticleSystem::update(real_t dt) {
    const real_t t = G.t.time;
    const real_t dt_scaled = dt * m_simulation_speed;

    // Apply gravity
    m_densities.resize(m_particles_count);
    for (index_t i = 0; i < m_particles_count; i++) {
        m_velocities[i] += vec3(0, 0, 1) * m_gravity * dt_scaled;
        m_predicted_positions[i] = m_positions[i] + m_velocities[i] * dt_scaled;
    }

    // Compute densities
    for (index_t i = 0; i < m_particles_count; i++) {
        m_densities[i] = 0;
        // Calculate density
        for (index_t j = 0; j < m_particles_count; j++) {
            real_t dst = glm::distance(m_predicted_positions[i], m_predicted_positions[j]);
            m_densities[i] += density_kernel(dst, smoothing_radius);
            m_near_densities[i] += near_density_kernel(dst, smoothing_radius);
        }
    }

    // Apply pressure
    for (index_t i = 0; i < m_particles_count; i++) {
        vec3 pressure_force = calculate_pressure(i);
        vec3 pressure_acceleration = pressure_force / m_densities[i];
        m_velocities[i] += pressure_acceleration * dt_scaled;
    }

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

    ImGui::DragFloat("Smoothing Radius", &smoothing_radius, 0.01f, 1.f);
    ImGui::DragFloat("Target Density", &m_target_density, 0.001f, 10.f);
    ImGui::DragFloat("Pressure Multiplier", &m_pressure_multiplier, 0.1f, 0.001f, 10.f);
    ImGui::DragFloat("Near Pressure Multiplier", &m_near_pressure_multiplier, 0.1f, 0.001f, 10.f);

    ImGui::End();
}

const std::vector<vec3> &ParticleSystem::positions() const { return m_positions; }
const std::vector<vec3> &ParticleSystem::velocities() const { return m_velocities; }
const std::vector<real_t> &ParticleSystem::densities() const { return m_densities; }
vec3 &ParticleSystem::extent_ref() { return m_extents; }
