#include "Transform.hpp"
#include "glm/ext/quaternion_common.hpp"
#include "glm/gtc/quaternion.hpp"
#include "typedefs.hpp"

#include <glm/gtc/matrix_transform.hpp>

Transform::Transform()
    : m_parent(nullptr), m_position{0.f}, m_scale{1.f}, m_rotation(1, 0, 0, 0){};

Transform::Transform(Transform *parent)
    : m_parent{parent}, m_position{0.f}, m_scale{1.f}, m_rotation(1, 0, 0, 0) {
    if (m_parent) {
        m_parent->add_child(this);
    }
}

const mat4 Transform::world_matrix() {
    const auto local = local_matrix();
    if (m_world_dirty) {
        m_world_dirty = false;
        if (m_parent) {
            m_world_matrix = m_parent->world_matrix() * local;
        } else {
            m_world_matrix = local;
        }
    }
    return m_world_matrix;
}

const mat4 Transform::inverse_world_matrix() {
    const auto world = world_matrix();
    return glm::inverse(world);
}

const mat4 Transform::local_matrix() {
    if (m_local_dirty) {
        m_local_dirty = false;
        invalidate_children();
        mat4 translation = glm::translate(mat4(1.0), m_position);
        mat4 rotation_matrix = glm::mat4_cast(m_rotation);
        mat4 scaling = glm::scale(mat4(1.0f), m_scale);
        m_local_matrix = translation * rotation_matrix * scaling;
    }
    return m_local_matrix;
}

void Transform::invalidate_children() const {
    for (auto &child : m_children) {
        if (child) {
            child->m_world_dirty = true;
        }
    }
}

void Transform::set_position(const vec3 &position) {
    m_position = position;
    m_local_dirty = true;
    m_world_dirty = true;
}

void Transform::set_scale(const vec3 &scale) {
    m_scale = scale;
    m_local_dirty = true;
    m_world_dirty = true;
}

void Transform::set_rotation(const quat &rotation) {
    m_rotation = rotation;
    m_local_dirty = true;
    m_world_dirty = true;
}

void Transform::translate(const vec3 &translation) {
    m_position += translation;
    m_local_dirty = true;
    m_world_dirty = true;
}

void Transform::scale(const vec3 &scaling) {
    m_scale *= scaling;
    m_local_dirty = true;
    m_world_dirty = true;
}

void Transform::rotate(const quat &rotation) {
    m_rotation = rotation * m_rotation;
    m_local_dirty = true;
    m_world_dirty = true;
}

const vec3 Transform::position() const noexcept { return m_position; }
const vec3 Transform::scale() const noexcept { return m_scale; }
const quat Transform::rotation() const noexcept { return m_rotation; }

void Transform::set_parent(Transform *parent) {
    if (m_parent) {
        m_parent->m_children.erase(
            std::remove(m_parent->m_children.begin(), m_parent->m_children.end(), this),
            m_parent->m_children.end());
    }
    m_parent = parent;
    if (m_parent) {
        m_parent->add_child(this);
    }
}

void Transform::add_child(Transform *child) {
    if (child) {
        m_children.push_back(child);
        child->m_parent = this;
        child->m_world_dirty = true;
    }
}
