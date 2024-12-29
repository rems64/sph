#pragma once

#include <cstddef>
#include <glm/ext/vector_double3.hpp>
#include <glm/ext/vector_float3.hpp>

typedef float real_t;
typedef size_t index_t;
typedef glm::vec<2, real_t> vec2;
typedef glm::vec<3, real_t> vec3;
typedef glm::vec<4, real_t> vec4;
typedef glm::vec<3, index_t> ivec3;
typedef glm::mat<3, 3, real_t> mat3;
typedef glm::mat<4, 4, real_t> mat4;
typedef glm::qua<real_t, glm::defaultp> quat;

#define M_PI 3.141592653
