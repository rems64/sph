#include "Collider.hpp"
#include "Transform.hpp"
#include "typedefs.hpp"
#include <iostream>

Collider::Collider(handle<Transform> transform) : m_transform(transform) {}

handle<Transform> Collider::transform() const { return m_transform; }

BoxCollider::BoxCollider(handle<Transform> transform, vec3 extents)
    : Collider(transform), m_extents(extents) {}

vec3 &BoxCollider::extents() { return m_extents; }

Collision BoxCollider::collide(const vec3 &point, const vec3 &direction_point) const {
    vec4 tmp = m_transform.lock()->inverse_world_matrix() * vec4(point, 1);
    vec3 point_localspace = vec3(tmp);
    if (std::abs(tmp.w) > 0.0001f)
        point_localspace /= tmp.w;
    tmp = m_transform.lock()->inverse_world_matrix() * vec4(direction_point, 1);
    vec3 direction_point_localspace = vec3(tmp);
    if (std::abs(tmp.w) > 0.0001f)
        direction_point_localspace /= tmp.w;
    const vec3 min = -glm::abs(m_extents) / 2.f;
    const vec3 max = glm::abs(m_extents) / 2.f;

    // Test collision using SAT
    bool collided = (point_localspace.x > min.x && point_localspace.x < max.x) &&
                    (point_localspace.y > min.y && point_localspace.y < max.y) &&
                    (point_localspace.z > min.z && point_localspace.z < max.z);

    if (!collided) {
        return Collision{false, vec3(0)};
    }

    // Calculate minimum translation vector
    real_t diffs[3] = {
        std::min(std::abs(point_localspace.x - min.x), std::abs(point_localspace.x - max.x)),
        std::min(std::abs(point_localspace.y - min.y), std::abs(point_localspace.y - max.y)),
        std::min(std::abs(point_localspace.z - min.z), std::abs(point_localspace.z - max.z))};
    vec3 vecs[3] = {vec3(1, 0, 0) * diffs[0] * glm::sign(direction_point_localspace.x) *
                        glm::sign(point_localspace.x * direction_point_localspace.x),
                    vec3(0, 1, 0) * diffs[1] * glm::sign(direction_point_localspace.y) *
                        glm::sign(point_localspace.y * direction_point_localspace.y),
                    vec3(0, 0, 1) * diffs[2] * glm::sign(direction_point_localspace.z) *
                        glm::sign(point_localspace.z * direction_point_localspace.z)};

    index_t best_axis[3] = {0, 1, 2};
    // Order by penetration depth
    if (diffs[best_axis[0]] > diffs[best_axis[1]]) {
        std::swap(best_axis[0], best_axis[1]);
    }
    if (diffs[best_axis[0]] > diffs[best_axis[2]]) {
        std::swap(best_axis[0], best_axis[2]);
    }
    if (diffs[best_axis[1]] > diffs[best_axis[2]]) {
        std::swap(best_axis[1], best_axis[2]);
    }

    return Collision{true, vecs[best_axis[0]], vecs[best_axis[1]]};
}
