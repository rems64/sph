#pragma once

#include "Mesh.hpp"
#include "Transform.hpp"
#include "mem.hpp"
#include "typedefs.hpp"
#include <utility>
#include <vector>

class Material;
class Mesh;
class ParticleSystem;
class ShaderProgram;

class Renderer {
public:
    Renderer();
    virtual ~Renderer() = default;

    virtual void render() = 0;
};

using MeshRendererItem = std::tuple<handle<Material>, handle<Mesh>, handle<Transform>>;

class MeshRenderer : public Renderer {
public:
    MeshRenderer();
    virtual ~MeshRenderer() = default;

    void render() override;
    void compute_clusters();

    void add(const handle<Material> &material,
             const handle<Mesh> &mesh,
             const handle<Transform> &transform);

private:
    std::vector<MeshRendererItem> m_items;
    std::vector<std::vector<MeshRendererItem>> m_clusters;
    bool m_clusters_dirty = true;
};

class ParticleRenderer : public Renderer {
public:
    ParticleRenderer();
    virtual ~ParticleRenderer() = default;

    void render() override;

    void update_positions(const index_t index);
    void update_colors(const index_t index);
    void add(const handle<ParticleSystem> &particle_system);

private:
    handle<ShaderProgram> m_shader;
    std::vector<handle<ParticleSystem>> m_items;

private:
    handle<ParticleMesh> m_mesh;
    std::vector<vec3> m_colors;
    std::vector<vec3> m_positions;
    GLuint m_color_vbo = 0;
    GLuint m_position_vbo = 0;
};

glm::vec3 speed_map(const real_t speed);
