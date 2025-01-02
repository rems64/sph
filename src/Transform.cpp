#include "Transform.hpp"
#include "ImGuizmo.h"
#include "glm/ext/quaternion_common.hpp"
#include "glm/gtc/quaternion.hpp"
#include "glm/gtc/type_ptr.hpp"
#include "typedefs.hpp"
#include <iostream>

#include <glm/gtc/matrix_transform.hpp>

Transform::Transform()
    : m_parent(nullptr), m_position(0.f), m_scale(1.f), m_rotation(1, 0, 0, 0),
      m_local_dirty(true), m_world_dirty(true), m_local_matrix(1.f), m_world_matrix(1.f){};

Transform::Transform(Transform *parent)
    : m_parent(parent), m_position(0.f), m_scale(1.f), m_rotation(1, 0, 0, 0), m_local_dirty(true),
      m_world_dirty(true), m_local_matrix(1.f), m_world_matrix(1.f) {
    if (m_parent) {
        m_parent->add_child(this);
    }
}

mat4 Transform::world_matrix() {
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

mat4 Transform::inverse_world_matrix() {
    const auto world = world_matrix();
    return glm::inverse(world);
}

const mat4 Transform::local_matrix() {
    if (m_local_dirty) {
        m_world_dirty = true;
        m_local_dirty = false;
        invalidate_children();
        mat4 translation = glm::translate(mat4(1.0), m_position);
        mat4 rotation_matrix = glm::mat4_cast(m_rotation);
        mat4 scaling = glm::scale(mat4(1.0f), m_scale);
        m_local_matrix = translation * rotation_matrix * scaling;
    }
    return m_local_matrix;
}

real_t *Transform::local_matrix_ptr() {
    local_matrix();
    return glm::value_ptr(m_local_matrix);
}

void Transform::rebuild_from_matrix() {
    vec3 translation;
    vec3 rotation;
    vec3 scale;
    ImGuizmo::DecomposeMatrixToComponents(glm::value_ptr(m_local_matrix),
                                          glm::value_ptr(translation),
                                          glm::value_ptr(rotation),
                                          glm::value_ptr(scale));
    m_position = translation;
    m_rotation = glm::quat(rotation);
    m_scale = scale;
}

void Transform::invalidate() {
    m_local_dirty = true;
    m_world_dirty = true;
}

void Transform::invalidate_children() const {
    for (auto &child : m_children) {
        if (child) {
            child->invalidate();
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

vec3 &Transform::position_ref() { return m_position; }
vec3 &Transform::scale_ref() { return m_scale; }
quat &Transform::rotation_ref() { return m_rotation; }

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
