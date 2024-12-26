#include "Entity.hpp"
#include "BoundingBox.hpp"
#include "mem.hpp"
#include <limits>

Entity::Entity()
    : m_mesh{}, m_transform{}, m_material{},
      m_bounding_box{make_ref<BoundingBox>(vec3(std::numeric_limits<float>::infinity()),
                                           vec3(-std::numeric_limits<float>::infinity()))} {};

void Entity::set_mesh(const std::shared_ptr<Mesh> &mesh) { m_mesh = mesh; }

void Entity::set_material(const std::shared_ptr<Material> &material) { m_material = material; }

void Entity::set_transform(const std::shared_ptr<Transform> &transform) {
    m_transform = transform;
}

void Entity::set_bounding_box(const std::shared_ptr<BoundingBox> &bounding_box) {
    m_bounding_box = bounding_box;
}
