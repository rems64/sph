#include <cstdlib>
#include <glm/ext/vector_float3.hpp>

#include "BoundingBox.hpp"
#include "Globals.hpp"
#include "Material.hpp"
#include "Mesh.hpp"
#include "Renderer.hpp"
#include "ResourceManager.hpp"
#include "ShaderProgram.hpp"
#include "UI.hpp"
#include "Window.hpp"
#include "imgui.h"

void build_scene() {
    auto resource_manager = G.resource_manager.lock();

    const handle<MeshRenderer> mesh_renderer = resource_manager->add(make_ref<MeshRenderer>());
    const handle<ParticleRenderer> particle_renderer = resource_manager->add(
        make_ref<ParticleRenderer>());
    const handle<ShaderProgram> shader = resource_manager->add(
        ShaderProgram::genBasicShaderProgram("shaders/basic.vert", "shaders/basic.frag"));
    const handle<Mesh> mesh = resource_manager->build_mesh();
    const handle<Material> material = resource_manager->build_material(
        shader, glm::vec3(0.5f, 0.3f, 0.1f));
    const handle<BoundingBox> bounding_box = resource_manager->build_boundingbox();
    const handle<Transform> transform = resource_manager->build_transform();
    const handle<ParticleSystem> particle_system = resource_manager->build_particle_system();

    mesh.lock()->addBox(1.f, 1.f, 1.f);
    // mesh.lock()->addPlane(1.f);
    mesh.lock()->init();

    mesh_renderer.lock()->add(material, mesh);
    particle_renderer.lock()->add(particle_system);
}

void init() {
    init_windowing();
    G.application_flags = 0;
    G.application_flags |= APP_RUNNING;

    G.window = G.resource_manager.lock()->build_window(800, 600, "Smooth Particle Hydrodynamics");

    build_scene();

    ui_init(G.window.lock()->raw_window());
}

void close() { close_windowing(); }

void update(real_t dt) {
    const auto resource_manager = G.resource_manager.lock();
    const auto particle_systems = resource_manager->particle_systems();
    for (const auto &particle_system : particle_systems) {
        particle_system->update(dt);
    }
}

void render() {
    ui_new_frame();

    const auto window = G.window.lock();
    auto resource_manager = G.resource_manager.lock();
    const auto renderers = resource_manager->renderers();
    glEnable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glViewport(0, 0, window->width(), window->height());
    for (const auto &renderer : renderers) {
        renderer->render();
    }

    G.resource_manager.lock()->particle_systems()[0]->imgui_controls();

    ui_render();
}

void frame(real_t dt) {
    update(dt);
    render();
}

int main(int argc, char **argv) {
    {
        const auto resource_manager = make_resource_manager();
        G.resource_manager = resource_manager;

        init();

        while (G.application_flags & APP_RUNNING) {
            auto window = G.window.lock();
            if (window->should_close()) {
                G.application_flags &= ~APP_RUNNING;
            }

            const real_t current_time = static_cast<real_t>(get_time());
            const real_t dt = current_time - G.t.time;
            G.t.time = current_time;

            frame(dt);

            window->swap_buffers();
            window->poll_events();
        }
    }

    close();
    return EXIT_SUCCESS;
}
