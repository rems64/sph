#include "Collider.hpp"
#include "Transform.hpp"
#include "typedefs.hpp"

Collider::Collider(handle<Transform> transform) : m_transform(transform) {}

handle<Transform> Collider::transform() const { return m_transform; }

BoxCollider::BoxCollider(handle<Transform> transform, vec3 extents)
    : Collider(transform), m_extents(extents) {}

vec3 &BoxCollider::extents() { return m_extents; }

Collision BoxCollider::collide(const vec3 &point) const {

    const vec3 point_localspace = m_transform.lock()->inverse_world_matrix() * vec4(point, 1);
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
    const real_t x_diff = std::min(std::abs(point_localspace.x - min.x),
                                   std::abs(point_localspace.x - max.x));
    const real_t y_diff = std::min(std::abs(point_localspace.y - min.y),
                                   std::abs(point_localspace.y - max.y));
    const real_t z_diff = std::min(std::abs(point_localspace.z - min.z),
                                   std::abs(point_localspace.z - max.z));

    if (x_diff < y_diff && x_diff < z_diff) {
        if (point_localspace.x < 0) {
            return Collision{collided, vec3(-x_diff, 0, 0)};
        } else {
            return Collision{collided, vec3(x_diff, 0, 0)};
        }
    } else if (y_diff < x_diff && y_diff < z_diff) {
        if (point_localspace.y < 0) {
            return Collision{collided, vec3(0, -y_diff, 0)};
        } else {
            return Collision{collided, vec3(0, y_diff, 0)};
        }
    } else {
        if (point_localspace.z < 0) {
            return Collision{collided, vec3(0, 0, -z_diff)};
        } else {
            return Collision{collided, vec3(0, 0, z_diff)};
        }
    }
}
