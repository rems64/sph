#include "Renderer.hpp"

#include "Globals.hpp"
#include "Material.hpp"
#include "Mesh.hpp"
#include "ParticleSystem.hpp"
#include "ResourceManager.hpp"
#include "ShaderProgram.hpp"
#include "glm/ext/matrix_clip_space.hpp"
#include "glm/ext/matrix_transform.hpp"
#include "glm/trigonometric.hpp"
#include "typedefs.hpp"
#include <iostream>
#include <memory>

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
        shader->set("projMat",
                    glm::perspective(
                        glm::radians(60.f), G.window.lock()->get_aspect_ratio(), 0.1f, 100.f));
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

ParticleRenderer::ParticleRenderer() {
    // TODO: lazy create
    m_mesh = G.resource_manager.lock()->build_particle_mesh();
    m_mesh.lock()->init();
    m_shader = G.resource_manager.lock()->add(
        ShaderProgram::genBasicShaderProgram("shaders/particle.vert", "shaders/particle.frag"));
}

void ParticleRenderer::render() {
    if (m_mesh.expired())
        return;

    update_positions(0);
    update_colors(0);

    const auto shader = m_shader.lock();
    shader->use();
    shader->set("viewMat", glm::lookAt(glm::vec3(-3., -3., 3.), glm::vec3(0), glm::vec3(0, 0, 1)));
    shader->set(
        "projMat",
        glm::perspective(glm::radians(60.f), G.window.lock()->get_aspect_ratio(), 0.1f, 100.f));

    for (const auto &item : m_items) {
        const auto world_matrix = glm::scale(mat4(1.0f), 1.f * vec3(1.f));
        shader->set("particle_scale", 0.025f);
        shader->set("modelMat", world_matrix);
        shader->set("normMat", glm::mat3(glm::inverseTranspose(world_matrix)));

        const auto positions = item.lock()->positions();
        glBindVertexArray(m_mesh.lock()->vao());
        // glDrawArraysInstanced(GL_TRIANGLES, 0, 24, positions.size());
        glDrawElementsInstanced(GL_TRIANGLES,
                                static_cast<GLsizei>(m_mesh.lock()->triangles_count() * 3),
                                GL_UNSIGNED_INT,
                                0,
                                positions.size());
    }
}

void ParticleRenderer::update_positions(const index_t index) {
    const auto system = m_items[index].lock();
    const auto positions = system->positions();
    m_positions = positions;
    m_colors.resize(positions.size());
    glNamedBufferSubData(m_position_vbo, 0, sizeof(vec3) * m_positions.size(), m_positions.data());
}

void ParticleRenderer::update_colors(const index_t index) {
    const auto system = m_items[index].lock();
    const auto &speeds = system->velocities();
    const auto &densities = system->densities();
    m_colors.resize(speeds.size());
    real_t max_speed = system->max_velocity();
    for (index_t i = 0; i < speeds.size(); i++) {
        const real_t density_error = densities[i] - system->m_target_density;
        // const real_t normalized_speed = glm::length(speeds[i]) / max_speed;
        // m_colors[i] = glm::vec3(normalized_speed, .5f, .5f);
        m_colors[i] = speed_map(density_error/100.f);
    }
    glNamedBufferSubData(m_color_vbo, 0, sizeof(vec3) * m_colors.size(), m_colors.data());
}

void ParticleRenderer::add(const handle<ParticleSystem> &particle_system) {
    m_items.push_back(particle_system);

    const auto system = particle_system.lock();
    const auto positions = system->positions();
    m_positions.resize(positions.size());
    m_colors.resize(positions.size());

    m_positions = positions;

    glCreateBuffers(1, &m_color_vbo);
    glNamedBufferStorage(
        m_color_vbo, sizeof(vec3) * m_colors.size(), m_colors.data(), GL_DYNAMIC_STORAGE_BIT);

    glCreateBuffers(1, &m_position_vbo);
    glNamedBufferStorage(m_position_vbo,
                         sizeof(vec3) * m_positions.size(),
                         m_positions.data(),
                         GL_DYNAMIC_STORAGE_BIT);

    glBindVertexArray(m_mesh.lock()->vao());
    glEnableVertexAttribArray(2);
    glBindBuffer(GL_ARRAY_BUFFER, m_color_vbo);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(vec3), 0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glVertexAttribDivisor(2, 1);

    glEnableVertexAttribArray(3);
    glBindBuffer(GL_ARRAY_BUFFER, m_position_vbo);
    glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(vec3), 0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glVertexAttribDivisor(3, 1);
}

glm::vec3 lerp(const glm::vec3 &a, const glm::vec3 &b, const real_t t) { return a + t * (b - a); }

real_t clamp(const real_t x, const real_t a, const real_t b) {
    return glm::min(glm::max(x, a), b);
}

glm::vec3 hsv_to_rgb(const glm::vec3 &hsv) {
    const real_t c = hsv.y * hsv.z;
    const real_t h = hsv.x / 60.f;
    const real_t x = c * (1 - glm::abs(glm::mod(h, 2.f) - 1));
    const real_t m = hsv.z - c;

    if (h >= 0 && h < 1) {
        return glm::vec3(c + m, x + m, m);
    } else if (h >= 1 && h < 2) {
        return glm::vec3(x + m, c + m, m);
    } else if (h >= 2 && h < 3) {
        return glm::vec3(m, c + m, x + m);
    } else if (h >= 3 && h < 4) {
        return glm::vec3(m, x + m, c + m);
    } else if (h >= 4 && h < 5) {
        return glm::vec3(x + m, m, c + m);
    } else if (h >= 5 && h < 6) {
        return glm::vec3(c + m, m, x + m);
    } else {
        return glm::vec3(m);
    }
}

glm::vec3 speed_map(const real_t speed) {
    const real_t factor = clamp(speed, 0.f, 1.f);
    const std::vector<glm::vec3> colors_hsv = {{240, 1, 1}, {60, 1, 1}, {0, 1, 1}};
    // Lerp between the colors in hsv space
    const glm::vec3 color_hsv = lerp(
        lerp(colors_hsv[0], colors_hsv[1], factor), colors_hsv[2], factor);
    return hsv_to_rgb(color_hsv);
}
