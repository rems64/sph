#include "ParticleSystem.hpp"
#include "Globals.hpp"
#include "glm/common.hpp"
#include "glm/geometric.hpp"
#include "typedefs.hpp"
#include <cmath>
#include <cstdlib>

ParticleSystem::ParticleSystem() : m_positions{}, m_extents{2.f, 2.f, 2.f} {
    m_positions.resize(10000);
    m_velocities.resize(m_positions.size());

    // Arrange the particles in a circle
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
        m_positions[i] += m_velocities[i] * dt;
        if (std::abs(m_positions[i].x) > m_extents.x / 2.f) {
            m_positions[i].x = glm::sign(m_positions[i].x) * m_extents.x / 2.f;
            m_velocities[i] = glm::reflect(m_velocities[i], vec3(1, 0, 0));
        }
        if (std::abs(m_positions[i].y) > m_extents.y / 2.f) {
            m_positions[i].y = glm::sign(m_positions[i].y) * m_extents.y / 2.f;
            m_velocities[i] = glm::reflect(m_velocities[i], vec3(0, 1, 0));
        }
        if (std::abs(m_positions[i].z) > m_extents.z / 2.f) {
            m_positions[i].z = glm::sign(m_positions[i].z) * m_extents.z / 2.f;
            m_velocities[i] = glm::reflect(m_velocities[i], vec3(0, 0, 1));
        }
    }
}

const std::vector<vec3> &ParticleSystem::positions() const { return m_positions; }

const std::vector<vec3> &ParticleSystem::velocities() const { return m_velocities; }
