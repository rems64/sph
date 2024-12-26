#include <cstdlib>
#include <glm/ext/vector_float3.hpp>

#include "BoundingBox.hpp"
#include "Globals.hpp"
#include "Material.hpp"
#include "Renderer.hpp"
#include "ResourceManager.hpp"
#include "ShaderProgram.hpp"
#include "Window.hpp"

void build_scene() {
    auto resource_manager = G.resource_manager.lock();

    const handle<MeshRenderer> mesh_renderer = resource_manager->add(make_ref<MeshRenderer>());
    const handle<ShaderProgram> shader = resource_manager->add(
        ShaderProgram::genBasicShaderProgram("shaders/basic.vert", "shaders/basic.frag"));
    const handle<Mesh> mesh = resource_manager->build_mesh();
    const handle<Material> material = resource_manager->build_material(
        shader, glm::vec3(0.5f, 0.3f, 0.1f));
    const handle<BoundingBox> bounding_box = resource_manager->build_boundingbox();
    const handle<Transform> transform = resource_manager->build_transform();

    // mesh.lock()->addBox(1.f, 1.f, 1.f);
    mesh.lock()->addPlane(1.f);
    mesh.lock()->init();

    mesh_renderer.lock()->add(material, mesh);
}

void init() {
    init_windowing();
    G.application_flags = 0;
    G.application_flags |= APP_RUNNING;

    G.window = G.resource_manager.lock()->build_window(800, 600, "Smooth Particle Hydrodynamics");

    build_scene();
}

void close() { close_windowing(); }

void update(real_t dt) {}

void render() {
    auto resource_manager = G.resource_manager.lock();
    for (const auto &renderer : resource_manager->renderers()) {
        renderer->render();
    }
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
