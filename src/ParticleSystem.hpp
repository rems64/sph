#pragma once

#include "typedefs.hpp"
#include <vector>

class ParticleSystem {
public:
    ParticleSystem();

    void update(real_t dt);

    void resolve_collision(index_t i);

    const vec3 calculate_pressure(index_t i);
    const real_t convert_density_to_pressure(real_t density);
    const real_t convert_near_density_to_near_pressure(real_t density);

    const std::vector<vec3> &positions() const;
    const std::vector<vec3> &velocities() const;
    const std::vector<real_t> &densities() const;

    const float max_velocity();

    vec3 &extent_ref();

    void imgui_controls();

    friend class ParticleRenderer;

private:
    size_t m_particles_count;
    std::vector<vec3> m_positions;
    std::vector<vec3> m_predicted_positions;
    std::vector<vec3> m_velocities;
    std::vector<real_t> m_near_densities;
    std::vector<real_t> m_densities;

private:
    vec3 m_extents = {.6, .6, 1.};
    float m_collision_damping = .05f;
    float m_collision_tangent_damping = 1.f;
    float m_gravity = -.0f;
    real_t smoothing_radius = 0.2f;
    real_t m_target_density = 400.f;
    real_t m_pressure_multiplier = .1f;
    real_t m_near_pressure_multiplier = .1f;

    real_t m_max_velocity = 10.f;

    real_t m_simulation_speed = 1.f;
};
