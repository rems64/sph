#pragma once

#include "typedefs.hpp"
#include <vector>

class ParticleSystem {
public:
    ParticleSystem();

    void update(real_t dt);
    void render();

    const std::vector<vec3> &positions() const;
    const std::vector<vec3> &speeds() const;

private:
    std::vector<vec3> m_positions;
    std::vector<vec3> m_speeds;
};
