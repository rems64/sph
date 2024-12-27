#pragma once

#include "typedefs.hpp"
#include <vector>

class ParticleSystem {
public:
    ParticleSystem();

    void update(real_t dt);

    void resolve_collision(index_t i);

    const std::vector<vec3> &positions() const;
    const std::vector<vec3> &velocities() const;

    const float max_velocity();

    vec3 &extent_ref();

    void imgui_controls();

private:
    std::vector<vec3> m_positions;
    std::vector<vec3> m_velocities;

private:
    vec3 m_extents;
    float m_collision_damping = 1.f;
    float m_collision_tangent_damping = 1.f;
    float m_gravity = -9.81f;

    float m_max_velocity = 10.f;
};
