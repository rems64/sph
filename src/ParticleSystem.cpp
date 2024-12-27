#include "ParticleSystem.hpp"
#include "Globals.hpp"
#include "typedefs.hpp"
#include <cmath>

ParticleSystem::ParticleSystem() : m_positions{} {
    m_positions.resize(200);
    m_speeds.resize(m_positions.size());

    // Arrange the particles in a circle
    for (index_t i = 0; i < m_positions.size(); i++) {
        m_positions[i] = vec3(std::cos(2 * M_PI * i / m_positions.size()),
                              std::sin(2 * M_PI * i / m_positions.size()),
                              0);
        m_speeds[i] = vec3(0.f);
    }
}

void ParticleSystem::update(real_t dt) {
    // Move the particles in a circle
    const real_t t = G.t.time;
    const real_t speed = .05f;
    for (index_t i = 0; i < m_positions.size(); i++) {
        const real_t f = 2 * M_PI * i / m_positions.size();
        const vec3 new_position = vec3(
            std::cos(f + t * speed), std::sin(f + t * speed), std::sin(5 * f + 5 * t * speed));
        m_speeds[i] = (new_position - m_positions[i]) / dt;
        m_positions[i] = new_position;
    }
}

const std::vector<vec3> &ParticleSystem::positions() const { return m_positions; }

const std::vector<vec3> &ParticleSystem::speeds() const { return m_speeds; }