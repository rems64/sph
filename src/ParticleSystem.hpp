#pragma once

#include "typedefs.hpp"
#include <vector>

class ParticleSystem {
public:
    ParticleSystem();

    void update(real_t dt);
    void render();

    const std::vector<vec3> &positions() const;
    const std::vector<vec3> &velocities() const;

private:
    std::vector<vec3> m_positions;
    std::vector<vec3> m_velocities;

private:
    vec3 m_extents;
};
