#pragma once

#include "mem.hpp"
#include "typedefs.hpp"

class Transform;

struct Collision {
    bool collided;
    vec3 mtv;
};

class Collider {
public:
    Collider(handle<Transform> transform);
    virtual ~Collider() = default;

    handle<Transform> transform() const;

    virtual Collision collide(const vec3 &point) const = 0;

protected:
    handle<Transform> m_transform;
};

class BoxCollider : public Collider {
public:
    BoxCollider(handle<Transform> transform, vec3 extents);

    vec3 &extents();

    Collision collide(const vec3 &point) const override;

private:
    vec3 m_extents;
};
