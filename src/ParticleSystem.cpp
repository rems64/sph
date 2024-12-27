#include "ParticleSystem.hpp"
#include "Globals.hpp"
#include "glm/common.hpp"
#include "glm/geometric.hpp"
#include "typedefs.hpp"
#include <cmath>
#include <cstdlib>

#include "imgui.h"

static real_t smoothing_kernel(real_t dst, real_t radius) {
    if (dst >= radius)
        return 0;
    // volume in 2d
    // real_t volume = (M_PI * pow(radius, 4)) / 6;
    // volume in 3d
    real_t volume = (32 * M_PI * M_PI * pow(radius, 7)) / 35.f;
    real_t v = radius * radius - dst * dst;
    return v * v * v / volume;
}

ParticleSystem::ParticleSystem() : m_positions{}, m_extents{2.f, 2.f, 2.f} {
    m_positions.resize(10000);
    m_velocities.resize(m_positions.size());

    for (index_t i = 0; i < m_positions.size(); i++) {
        m_positions[i] = vec3(0);
        m_velocities[i] = ((float)rand() / RAND_MAX) *
                          glm::normalize(vec3(((float)rand() / RAND_MAX - .5f),
                                              ((float)rand() / RAND_MAX - .5f),
                                              ((float)rand() / RAND_MAX - .5f)));
    }
}

void ParticleSystem::update(real_t dt) {
    // Move the particles in a circle
    const real_t t = G.t.time;
    for (index_t i = 0; i < m_positions.size(); i++) {
        m_velocities[i] += vec3(0, 0, 1) * m_gravity * dt;
        m_positions[i] += m_velocities[i] * dt;

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
    ImGui::SliderFloat("Extent X", &m_extents.x, 0.1f, 10.f);
    ImGui::SliderFloat("Extent Y", &m_extents.y, 0.1f, 10.f);
    ImGui::SliderFloat("Extent Z", &m_extents.z, 0.1f, 10.f);

    ImGui::SliderFloat("Max Velocity", &m_max_velocity, 0.1f, 100.f);

    ImGui::SliderFloat("Collision Damping", &m_collision_damping, 0.1f, 1.f);
    ImGui::SliderFloat("Collision Tangent Damping", &m_collision_tangent_damping, 0.1f, 1.f);

    ImGui::SliderFloat("Gravity", &m_gravity, -20.f, 20.f);
    ImGui::End();
}

const std::vector<vec3> &ParticleSystem::positions() const { return m_positions; }
const std::vector<vec3> &ParticleSystem::velocities() const { return m_velocities; }
vec3 &ParticleSystem::extent_ref() { return m_extents; }
