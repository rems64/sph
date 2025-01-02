#pragma once

#include "typedefs.hpp"
#include <glm/gtc/quaternion.hpp>

class Entity;

class Transform {
public:
    Transform();
    Transform(Transform *parent);

    void set_position(const vec3 &position);
    void set_scale(const vec3 &scale);
    void set_rotation(const quat &rotation);

    void translate(const vec3 &translation);
    void scale(const vec3 &scaling);
    void rotate(const quat &rotation);

    const vec3 position() const noexcept;
    const vec3 scale() const noexcept;
    const quat rotation() const noexcept;

    vec3 &position_ref();
    vec3 &scale_ref();
    quat &rotation_ref();

    const mat4 local_matrix();
    mat4 world_matrix();
    mat4 inverse_world_matrix();

    real_t *local_matrix_ptr();

    void set_parent(Transform *parent);
    void add_child(Transform *child);

    void invalidate();
    void invalidate_world();
    void rebuild_from_matrix();

private:
    void invalidate_children() const;

private:
    vec3 m_position;
    vec3 m_scale;
    quat m_rotation;

    bool m_local_dirty;
    bool m_world_dirty;

    mat4 m_world_matrix;
    mat4 m_local_matrix;

    Transform *m_parent;
    std::vector<Transform *> m_children;
};
