#include "Renderer.hpp"

#include "Material.hpp"
#include "Mesh.hpp"
#include "ShaderProgram.hpp"
#include "glm/ext/matrix_clip_space.hpp"
#include "glm/ext/matrix_transform.hpp"
#include "glm/trigonometric.hpp"
#include <memory>
#include <iostream>

Renderer::Renderer() {}

MeshRenderer::MeshRenderer() {}

void MeshRenderer::render() {
    if (m_clusters_dirty) {
        std::cout << "rebuilt clusters" << std::endl;
        compute_clusters();
    }
    for (const auto &cluster : m_clusters) {
        const auto shader = cluster.begin()->first.lock()->shader().lock();
        shader->use();

        shader->set("viewMat",
                    glm::lookAt(glm::vec3(-3., -3., 3.), glm::vec3(0), glm::vec3(0, 0, 1)));
        shader->set("projMat", glm::perspective(glm::radians(60.f), 1.f, 0.1f, 100.f));
        for (const auto &item : cluster) {
            const auto material = item.first.lock();
            const auto mesh = item.second.lock();
            const auto world_matrix = glm::mat4(1.f);
            
            shader->set("material.albedo", material->albedo());
            shader->set("modelMat", world_matrix);
            shader->set("normMat", glm::mat3(glm::inverseTranspose(world_matrix)));

            mesh->render();
        }
    }
}

void MeshRenderer::compute_clusters() {
    std::vector<handle<ShaderProgram>> shaders;
    m_clusters.clear();
    for (const auto &item : m_items) {
        const auto shader = item.first.lock()->shader();
        index_t i;
        for (i = 0; i < shaders.size(); i++) {
            // if (!shaders[i].owner_before(shader) && !shader.owner_before(shaders[i])) {
            if (!shaders[i].owner_before(shader) && !shader.owner_before(shaders[i])) {
                m_clusters[i].push_back(item);
                break;
            }
        }
        if (i == shaders.size()) {
            shaders.push_back(shader);
            m_clusters.push_back({item});
        }
    }
    m_clusters_dirty = false;
}

void MeshRenderer::add(const handle<Material> &material, const handle<Mesh> &mesh) {
    m_items.push_back(std::make_pair(material, mesh));
    m_clusters_dirty = true;
}
