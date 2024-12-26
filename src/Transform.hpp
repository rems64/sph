#pragma once

#include "typedefs.hpp"
#include <glm/ext/quaternion_float.hpp>

struct Transform {
    vec3 position = vec3{0.f};  // position
    vec3 scale = vec3{1.f};     // scale
    quat rotation = quat();     // rotation quaternion

    mat4 world_matrix() const;
};
