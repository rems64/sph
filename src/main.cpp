#include <cstdlib>
#include <glm/ext/vector_float3.hpp>

#include "BoundingBox.hpp"
#include "Camera.hpp"
#include "Globals.hpp"
#include "Material.hpp"
#include "Mesh.hpp"
#include "Renderer.hpp"
#include "ResourceManager.hpp"
#include "ShaderProgram.hpp"
#include "UI.hpp"
#include "Window.hpp"
#include "glm/ext/quaternion_common.hpp"
#include "glm/ext/quaternion_transform.hpp"
#include "glm/ext/quaternion_trigonometric.hpp"
#include "glm/geometric.hpp"
#include "glm/gtc/quaternion.hpp"
#include "glm/trigonometric.hpp"
#include "imgui.h"
#include "typedefs.hpp"
#include <GLFW/glfw3.h>
#include <iostream>

void build_scene() {
    auto resource_manager = G.resource_manager.lock();
    const handle<PerspectiveCamera> camera = resource_manager->build_perspective_camera(
        glm::radians(60.f), 0.1f, 100.f);
    const handle<Transform> camera_transform = resource_manager->build_transform();
    const handle<ParticleRenderer> particle_renderer = resource_manager->add(
        make_ref<ParticleRenderer>());
    const handle<ShaderProgram> shader = resource_manager->add(
        ShaderProgram::genBasicShaderProgram("shaders/basic.vert", "shaders/basic.frag"));

    const handle<BoundingBox> bounding_box = resource_manager->build_boundingbox();
    const handle<Transform> transform = resource_manager->build_transform();
    const handle<ParticleSystem> particle_system = resource_manager->build_particle_system();

    camera_transform.lock()->translate(vec3(3.f, -3.f, 3.f));
    camera_transform.lock()->rotate(glm::angleAxis(glm::radians(45.f), vec3(1, 0, 0)));
    camera_transform.lock()->rotate(glm::angleAxis(glm::radians(45.f), vec3(0, 0, 1)));

    G.camera = {.camera = camera, .camera_transform = camera_transform};

    particle_renderer.lock()->add(particle_system);
}

real_t lerp(real_t a, real_t b, real_t s) { return (1 - s) * a + s * b; }

void update_camera(real_t dt) {
    PerspectiveCamera *cam = (PerspectiveCamera *)G.camera.camera.lock().get();
    real_t &target_radius = cam->target_distance();
    real_t &target_pitch = cam->target_pitch();
    real_t &target_yaw = cam->target_yaw();
    real_t &radius = cam->distance();
    real_t &pitch = cam->pitch();
    real_t &yaw = cam->yaw();

    const real_t lerp_factor = 10.f * dt;

    radius = lerp(radius, target_radius, lerp_factor);
    pitch = lerp(pitch, target_pitch, lerp_factor);
    yaw = lerp(yaw, target_yaw, lerp_factor);

    const auto camera_transform = G.camera.camera_transform.lock();
    const vec3 view = vec3(cos(yaw) * cos(pitch), sin(yaw) * cos(pitch), sin(pitch));
    const vec3 camera_position = -radius * view;
    camera_transform->set_position(camera_position);
    camera_transform->set_rotation(glm::quatLookAt(view, vec3(0, 0, 1)));
}

void cursor_pos_callback(GLFWwindow *window, double xpos, double ypos) {
    glm::vec2 new_pos = glm::vec2(xpos, ypos);
    glm::vec2 delta = new_pos - G.input.cursor.position;

    if (G.input.cursor.pressed) {
        float dphi = -0.01f * delta.x;
        float dtheta = -0.01f * delta.y;

        PerspectiveCamera *cam = (PerspectiveCamera *)G.camera.camera.lock().get();

        real_t &pitch = cam->target_pitch();
        real_t &yaw = cam->target_yaw();

        pitch += dtheta;
        yaw += dphi;
    }
    G.input.cursor.position = new_pos;
}

void scroll_callback(GLFWwindow *window, double xoffset, double yoffset) {
    real_t factor = (pow(2, -yoffset) - 1) * 0.2f + 1;

    PerspectiveCamera *cam = (PerspectiveCamera *)G.camera.camera.lock().get();
    real_t &radius = cam->target_distance();

    radius *= factor;
}

void mouse_button_callback(GLFWwindow *window, int button, int action, int mods) {
    if (ImGui::GetIO().WantCaptureMouse)
        return;
    switch (button) {
    case GLFW_MOUSE_BUTTON_LEFT:
        G.input.cursor.pressed = action == GLFW_PRESS;
        break;

    default:
        break;
    }
}

void init() {
    init_windowing();
    G.application_flags = 0;
    G.application_flags |= APP_RUNNING;

    G.window = G.resource_manager.lock()->build_window(800, 600, "Smooth Particle Hydrodynamics");
    glfwSetCursorPosCallback(G.window.lock()->raw_window(), cursor_pos_callback);
    glfwSetMouseButtonCallback(G.window.lock()->raw_window(), mouse_button_callback);
    glfwSetScrollCallback(G.window.lock()->raw_window(), scroll_callback);

    build_scene();

    ui_init(G.window.lock()->raw_window());
}

void close() { close_windowing(); }

void update(real_t dt) {
    update_camera(dt);
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
    ImGui::Begin("Stats");
    ImGui::Text("FPS %.0f", 1.f / G.t.dt);
    ImGui::Text("dt %.1f", 1000 * G.t.dt);
    ImGui::End();

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
            real_t current_time = static_cast<real_t>(get_time());
            real_t dt = current_time - G.t.time;
            while (dt < 1.f / 30.f) {
                current_time = static_cast<real_t>(get_time());
                dt = current_time - G.t.time;
            }
            G.t.time = current_time;
            G.t.dt = dt;

            auto window = G.window.lock();
            if (window->should_close()) {
                G.application_flags &= ~APP_RUNNING;
            }

            frame(dt);

            window->swap_buffers();
            window->poll_events();
        }
    }

    close();
    return EXIT_SUCCESS;
}
