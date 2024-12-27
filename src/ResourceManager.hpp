#pragma once

#include "BoundingBox.hpp"
#include "Material.hpp"
#include "Mesh.hpp"
#include "ParticleSystem.hpp"
#include "Renderer.hpp"
#include "ShaderProgram.hpp"
#include "Transform.hpp"
#include "Window.hpp"
#include "mem.hpp"

class ResourceManager {
public:
    ResourceManager();

    handle<ShaderProgram> add(const ref<ShaderProgram> &resource);
    handle<Renderer> add(const ref<Renderer> &resource);
    handle<MeshRenderer> add(const ref<MeshRenderer> &resource);
    handle<ParticleRenderer> add(const ref<ParticleRenderer> &resource);
    handle<Window> add(const ref<Window> &resource);
    handle<Mesh> add(const ref<Mesh> &resource);
    handle<ParticleSystem> add(const ref<ParticleSystem> &resource);
    handle<Material> add(const ref<Material> &resource);
    handle<BoundingBox> add(const ref<BoundingBox> &resource);
    handle<Transform> add(const ref<Transform> &resource);

    template <typename... Args> handle<Window> build_window(Args &&...args);
    template <typename... Args> handle<Renderer> build_renderer(Args &&...args);
    template <typename... Args> handle<ShaderProgram> build_shader(Args &&...args);
    template <typename... Args> handle<Mesh> build_mesh(Args &&...args);
    template <typename... Args> handle<ParticleMesh> build_particle_mesh(Args &&...args);
    template <typename... Args> handle<ParticleSystem> build_particle_system(Args &&...args);
    template <typename... Args> handle<Material> build_material(Args &&...args);
    template <typename... Args> handle<BoundingBox> build_boundingbox(Args &&...args);
    template <typename... Args> handle<Transform> build_transform(Args &&...args);

    inline const std::vector<ref<Window>> &windows() const { return m_windows; }
    inline const std::vector<ref<Renderer>> &renderers() const { return m_renderers; }
    inline const std::vector<ref<ShaderProgram>> &shaders() const { return m_shaders; }
    inline const std::vector<ref<Mesh>> &meshes() const { return m_meshes; }
    inline const std::vector<ref<ParticleSystem>> &particle_systems() const {
        return m_particle_systems;
    }
    inline const std::vector<ref<Material>> &materials() const { return m_materials; }
    inline const std::vector<ref<BoundingBox>> &bounding_boxes() const { return m_bounding_boxes; }
    inline const std::vector<ref<Transform>> &transforms() const { return m_transforms; }

private:
    std::vector<ref<Window>> m_windows;
    std::vector<ref<Renderer>> m_renderers;
    std::vector<ref<ShaderProgram>> m_shaders;
    std::vector<ref<ParticleSystem>> m_particle_systems;
    std::vector<ref<Mesh>> m_meshes;
    std::vector<ref<Material>> m_materials;
    std::vector<ref<BoundingBox>> m_bounding_boxes;
    std::vector<ref<Transform>> m_transforms;

private:
    ref<ParticleMesh> m_particle_mesh;
};

ref<ResourceManager> make_resource_manager();

template <typename... Args> handle<Window> ResourceManager::build_window(Args &&...args) {
    const auto resource = make_ref<Window>(std::forward<Args>(args)...);
    m_windows.push_back(resource);
    return resource;
}
template <typename... Args> handle<Renderer> ResourceManager::build_renderer(Args &&...args) {
    const auto resource = make_ref<Renderer>(std::forward<Args>(args)...);
    m_renderers.push_back(resource);
    return resource;
}
template <typename... Args> handle<ShaderProgram> ResourceManager::build_shader(Args &&...args) {
    const auto resource = make_ref<ShaderProgram>(std::forward<Args>(args)...);
    m_shaders.push_back(resource);
    return resource;
}
template <typename... Args> handle<Mesh> ResourceManager::build_mesh(Args &&...args) {
    const auto resource = make_ref<Mesh>(std::forward<Args>(args)...);
    m_meshes.push_back(resource);
    return resource;
}
template <typename... Args>
handle<ParticleMesh> ResourceManager::build_particle_mesh(Args &&...args) {
    m_particle_mesh = make_ref<ParticleMesh>(std::forward<Args>(args)...);
    return m_particle_mesh;
}
template <typename... Args>
handle<ParticleSystem> ResourceManager::build_particle_system(Args &&...args) {
    const auto resource = make_ref<ParticleSystem>(std::forward<Args>(args)...);
    m_particle_systems.push_back(resource);
    return resource;
}
template <typename... Args> handle<Material> ResourceManager::build_material(Args &&...args) {
    const auto resource = make_ref<Material>(std::forward<Args>(args)...);
    m_materials.push_back(resource);
    return resource;
}
template <typename... Args>
handle<BoundingBox> ResourceManager::build_boundingbox(Args &&...args) {
    const auto resource = make_ref<BoundingBox>(std::forward<Args>(args)...);
    m_bounding_boxes.push_back(resource);
    return resource;
}
template <typename... Args> handle<Transform> ResourceManager::build_transform(Args &&...args) {
    const auto resource = make_ref<Transform>(std::forward<Args>(args)...);
    m_transforms.push_back(resource);
    return resource;
}
