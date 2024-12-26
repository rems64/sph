#include "ResourceManager.hpp"
#include "Mesh.hpp"
#include "ShaderProgram.hpp"
#include "mem.hpp"

ResourceManager::ResourceManager() {}

ref<ResourceManager> make_resource_manager() { return std::make_shared<ResourceManager>(); }

handle<Window> ResourceManager::add(const ref<Window> &resource) {
    m_windows.push_back(resource);
    return resource;
}
handle<Renderer> ResourceManager::add(const ref<Renderer> &resource) {
    m_renderers.push_back(resource);
    return resource;
}
handle<MeshRenderer> ResourceManager::add(const ref<MeshRenderer> &resource) {
    m_renderers.push_back(resource);
    return resource;
}
handle<ShaderProgram> ResourceManager::add(const ref<ShaderProgram> &resource) {
    m_shaders.push_back(resource);
    return resource;
}
handle<Mesh> ResourceManager::add(const ref<Mesh> &resource) {
    m_meshes.push_back(resource);
    return resource;
}
handle<Material> ResourceManager::add(const ref<Material> &resource) {
    m_materials.push_back(resource);
    return resource;
}
handle<BoundingBox> ResourceManager::add(const ref<BoundingBox> &resource) {
    m_bounding_boxes.push_back(resource);
    return resource;
}
