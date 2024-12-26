#pragma once

#include "mem.hpp"
#include <utility>
#include <vector>

class Material;
class Mesh;

class Renderer {
public:
    Renderer();
    virtual ~Renderer() = default;

    virtual void render() = 0;
};

using MeshRendererItem = std::pair<handle<Material>, handle<Mesh>>;

class MeshRenderer : public Renderer {
public:
    MeshRenderer();
    virtual ~MeshRenderer() = default;

    void render() override;
    void compute_clusters();

    void add(const handle<Material> &material, const handle<Mesh> &mesh);

private:
    std::vector<MeshRendererItem> m_items;
    std::vector<std::vector<MeshRendererItem>> m_clusters;
    bool m_clusters_dirty = true;
};
