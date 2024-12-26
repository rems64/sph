#pragma once

#include <memory>
#include <optional>

#include "BoundingBox.hpp"
#include "Material.hpp"
#include "Mesh.hpp"
#include "Transform.hpp"
#include "mem.hpp"

template <typename T> using Component = std::optional<handle<T>>;

class Entity {
public:
    Entity();

    void set_mesh(const std::shared_ptr<Mesh> &mesh);
    void set_material(const std::shared_ptr<Material> &material);
    void set_transform(const std::shared_ptr<Transform> &transform);
    void set_bounding_box(const std::shared_ptr<BoundingBox> &bounding_box);

    [[nodiscard]] inline Component<Mesh> mesh() { return m_mesh; }
    [[nodiscard]] inline Component<Transform> transform() { return m_transform; }
    [[nodiscard]] inline Component<Material> material() { return m_material; }
    [[nodiscard]] inline Component<BoundingBox> bounding_box() { return m_bounding_box; }

private:
    Component<Mesh> m_mesh;
    Component<Material> m_material;

    Component<Transform> m_transform;
    Component<BoundingBox> m_bounding_box;
};
