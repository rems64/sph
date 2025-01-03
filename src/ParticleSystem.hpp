#pragma once

#include "mem.hpp"
#include "typedefs.hpp"
#include <vector>

typedef ivec3 cell_t;
typedef index_t cellhash_t;

class Transform;

class ParticleSystem {
public:
    ParticleSystem(handle<Transform> transform);

    void update(real_t dt);

    bool resolve_collision(index_t i);

    inline const void apply_gravity(real_t dt);
    inline const void apply_pressure(real_t dt);
    inline const void compute_densities(real_t dt);
    inline const void integrate_and_collide(real_t dt);
    inline const vec3 calculate_pressure(index_t i);
    inline const real_t convert_density_to_pressure(real_t density);
    inline const real_t convert_near_density_to_near_pressure(real_t density);

    // Grid related optimisations
    cell_t get_cell(const vec3 position) const;
    void build_cells();
    std::vector<index_t> get_neighbors_particles(const index_t particle_index) const;

    const std::vector<vec3> &positions() const;
    const std::vector<vec3> &velocities() const;
    const std::vector<real_t> &densities() const;

    void set_container_speed(const vec3 &speed);
    void set_container_delta_angle(const quat &delta_angle);

    handle<Transform> transform() const;

    const float max_velocity();

    vec3 &extent_ref();

    void imgui_controls();

    friend class ParticleRenderer;

private:
    const cellhash_t hash_cell(cell_t cell) const;

private:
    size_t m_particles_count;
    handle<Transform> m_transform;
    std::vector<vec3> m_positions;
    std::vector<vec3> m_predicted_positions;
    std::vector<vec3> m_velocities;
    std::vector<real_t> m_near_densities;
    std::vector<real_t> m_densities;

    // Found in a video from Sebastian Lague
    std::vector<std::pair<cellhash_t, index_t>> m_spatial_lookup;
    std::vector<index_t> m_cell_start_index;

private:
    vec3 m_extents = {1., 1., 1.};
    float m_collision_damping = .05f;
    float m_collision_tangent_damping = 1.f;
    float m_gravity = -4.0f;
    real_t m_smoothing_radius = 0.2f;
    real_t m_target_density = 1200.f;
    real_t m_pressure_multiplier = 4.f;
    real_t m_near_pressure_multiplier = .1f;

    real_t m_max_velocity = 10.f;

    real_t m_simulation_speed = 1.f;
};
