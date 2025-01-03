#include <cstdlib>
#include <glm/ext/vector_float3.hpp>

#include "BoundingBox.hpp"
#include "Camera.hpp"
#include "Collider.hpp"
#include "Globals.hpp"
#include "ImGuizmo.h"
#include "Material.hpp"
#include "Mesh.hpp"
#include "Renderer.hpp"
#include "ResourceManager.hpp"
#include "ShaderProgram.hpp"
#include "UI.hpp"
#include "Window.hpp"
#include "glm/common.hpp"
#include "glm/ext/matrix_transform.hpp"
#include "glm/ext/quaternion_common.hpp"
#include "glm/ext/quaternion_transform.hpp"
#include "glm/ext/quaternion_trigonometric.hpp"
#include "glm/fwd.hpp"
#include "glm/geometric.hpp"
#include "glm/gtc/quaternion.hpp"
#include "glm/gtc/type_ptr.hpp"
#include "glm/matrix.hpp"
#include "glm/trigonometric.hpp"
#include "imgui.h"
#include "typedefs.hpp"
#include <GLFW/glfw3.h>
#include <memory>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/matrix_decompose.hpp>
#include <iostream>

const glm::mat4 yUpToZUp = glm::rotate(glm::mat4(1.0f), glm::radians(-90.0f), glm::vec3(1, 0, 0));
const glm::mat4 zUpToYUp = glm::rotate(glm::mat4(1.0f), glm::radians(90.0f), glm::vec3(1, 0, 0));

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
    const handle<ParticleSystem> particle_system = resource_manager->build_particle_system(
        transform);
    const handle<MeshRenderer> renderer = resource_manager->add(make_ref<MeshRenderer>());

    for (index_t i = 0; i < 2; i++) {
        const handle<Transform> collider_transform = resource_manager->build_transform();
        const handle<BoxCollider> collider = resource_manager->build_box_collider(
            collider_transform, vec3(1.f));
        const handle<Mesh> mesh = resource_manager->build_mesh();
        const handle<Material> material = resource_manager->build_material(shader,
                                                                           vec3(1.f, 1.f, 0.f));
        mesh.lock()->addBox(1.f, 1.f, 1.f);
        mesh.lock()->init();
        renderer.lock()->add(material, mesh, collider_transform);
    }

    resource_manager->transforms()[2]->translate(vec3(1.5f, 0.f, 0.f));
    resource_manager->transforms()[3]->translate(vec3(-1.5f, 0.f, 0.f));
    resource_manager->transforms()[2]->scale(vec3(1.1f));
    resource_manager->transforms()[3]->scale(vec3(1.1f));

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

    if (G.input.cursor.pressed && !ImGui::GetIO().WantCaptureMouse) {
        float dphi = -0.01f * delta.x;
        float dtheta = -0.01f * delta.y;

        PerspectiveCamera *cam = (PerspectiveCamera *)G.camera.camera.lock().get();

        real_t &pitch = cam->target_pitch();
        real_t &yaw = cam->target_yaw();

        pitch = glm::clamp(pitch + dtheta, (float)-M_PI / 2 + 0.01f, (float)M_PI / 2 - 0.01f);
        yaw = yaw + dphi;
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
    switch (button) {
    case GLFW_MOUSE_BUTTON_LEFT:
        G.input.cursor.pressed = action == GLFW_PRESS;
        break;

    default:
        break;
    }
}

void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods) {
    if (action == GLFW_PRESS) {
        switch (key) {
        case GLFW_KEY_ESCAPE:
            G.application_flags &= ~APP_RUNNING;
            break;
        case GLFW_KEY_SPACE:
            G.simulation.running = !G.simulation.running;
            break;
        case GLFW_KEY_C:
            G.resource_manager.lock()->particle_systems()[0]->transform().lock()->set_position(
                vec3(0.f));
            break;
        case GLFW_KEY_T:
            G.simulation.show_transform = !G.simulation.show_transform;
            break;
        case GLFW_KEY_B:
            G.simulation.show_bounds = !G.simulation.show_bounds;
            break;
        default:
            if (key <= GLFW_KEY_9 && key >= GLFW_KEY_1) {
                G.debug.active_transforms ^= 1 << (key - GLFW_KEY_1);
            } else if (key == GLFW_KEY_0)
                G.debug.active_transforms = 0;
            break;
        }
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
    glfwSetKeyCallback(G.window.lock()->raw_window(), key_callback);

    build_scene();

    ui_init(G.window.lock()->raw_window());
}

void close() { close_windowing(); }

void update(real_t dt) {
    update_camera(dt);
    const auto resource_manager = G.resource_manager.lock();
    if (G.simulation.running) {
        const auto particle_systems = resource_manager->particle_systems();
        for (const auto &particle_system : particle_systems) {
            particle_system->update(dt);
        }
    }
}

// void ToImGuizmo(float *dest, const mat4 &src) {
//     for (auto row = 0; row < 4; row++) {
//         for (auto col = 0; col < 4; col++)
//             dest[row * 4 + col] = src[col * 4 + row];
//     }
// }

// void FromImGuizmo(mat4 &dest, float *src) {
//     for (auto row = 0; row < 4; row++) {
//         for (auto col = 0; col < 4; col++)
//             (&dest.m00_)[col * 4 + row] = src[row * 4 + col];
//     }
// }

void render() {
    ui_new_frame();

    ImGuiIO &io = ImGui::GetIO();

    const auto window = G.window.lock();
    auto resource_manager = G.resource_manager.lock();
    const auto renderers = resource_manager->renderers();
    glEnable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glViewport(0, 0, window->width(), window->height());
    for (const auto &renderer : renderers) {
        renderer->render();
    }

    ImGuizmo::SetImGuiContext(ImGui::GetCurrentContext());

    bool useWindow = false;
    if (useWindow) {
        ImGui::Begin("Gizmo", 0);
        ImGuizmo::SetDrawlist();
    }
    float windowWidth = (float)ImGui::GetWindowWidth();
    float windowHeight = (float)ImGui::GetWindowHeight();

    if (!useWindow) {
        ImGuizmo::SetRect(0, 0, io.DisplaySize.x, io.DisplaySize.y);
    } else {
        ImGuizmo::SetRect(
            ImGui::GetWindowPos().x, ImGui::GetWindowPos().y, windowWidth, windowHeight);
    }

    ImGuizmo::BeginFrame();
    auto system = G.resource_manager.lock()->particle_systems()[0];
    auto transform = G.resource_manager.lock()->transforms()[1];
    transform->set_scale(system->extent_ref());
    auto local_matrix = transform->local_matrix();
    // float *trans = transform->local_matrix_ptr();
    const auto world_matrix = G.camera.camera_transform.lock()->inverse_world_matrix();
    const auto projection_matrix = G.camera.camera.lock()->projection();
    mat4 identity = mat4(1.f);
    index_t id = 42;
    ImGuizmo::PushID(id++);
    auto &extent = system->extent_ref();
    float bounds[] = {-0.5f, -0.5f, -0.5f, 0.5f, 0.5f, 0.5f};
    if (ImGuizmo::Manipulate(
            glm::value_ptr(world_matrix),
            glm::value_ptr(projection_matrix),
            (G.simulation.show_transform ? (ImGuizmo::TRANSLATE | ImGuizmo::ROTATE)
                                         : (ImGuizmo::OPERATION)0u) |
                (G.simulation.show_bounds ? ImGuizmo::BOUNDS : (ImGuizmo::OPERATION)0u),
            ImGuizmo::MODE::WORLD,
            glm::value_ptr(local_matrix),
            nullptr,
            nullptr,
            G.simulation.show_bounds ? bounds : nullptr)) {
        vec3 position;
        vec3 rotation;
        vec3 scale;
        ImGuizmo::DecomposeMatrixToComponents(glm::value_ptr(local_matrix),
                                              glm::value_ptr(position),
                                              glm::value_ptr(rotation),
                                              glm::value_ptr(scale));
        const auto old_position = transform->position();
        const auto old_rotation = transform->rotation();
        const auto new_rotation = glm::quat(glm::radians(rotation));
        transform->set_position(position);
        transform->set_rotation(new_rotation);
        extent = scale;
        system->set_container_speed((position - old_position) / G.t.dt);
        system->set_container_delta_angle(new_rotation * glm::inverse(old_rotation));
    }
    transform->set_scale(vec3(1.f));
    ImGuizmo::PopID();
    const auto colliders = resource_manager->colliders();
    index_t i = 0;
    for (const auto &collider : colliders) {
        if (!(G.debug.active_transforms & (1 << i++))) {
            continue;
        }
        const auto box_collider = std::dynamic_pointer_cast<BoxCollider>(collider);
        if (!box_collider) {
            continue;
        }
        auto transform = collider->transform().lock();
        // auto &extent = box_collider->extents();
        // transform->set_scale(extent);
        auto local_matrix = transform->local_matrix();
        float bounds[] = {-0.5f, -0.5f, -0.5f, 0.5f, 0.5f, 0.5f};
        ImGuizmo::PushID(id++);
        if (ImGuizmo::Manipulate(glm::value_ptr(world_matrix),
                                 glm::value_ptr(projection_matrix),
                                 ImGuizmo::TRANSLATE | ImGuizmo::ROTATE | ImGuizmo::SCALE,
                                 ImGuizmo::MODE::WORLD,
                                 glm::value_ptr(local_matrix),
                                 nullptr,
                                 nullptr,
                                 nullptr)) {
            vec3 position;
            vec3 rotation;
            vec3 scale;
            ImGuizmo::DecomposeMatrixToComponents(glm::value_ptr(local_matrix),
                                                  glm::value_ptr(position),
                                                  glm::value_ptr(rotation),
                                                  glm::value_ptr(scale));
            transform->set_position(position);
            transform->set_rotation(glm::quat(glm::radians(rotation)));
            transform->set_scale(scale);
            // box_collider->extents() = scale;
        }
        // transform->set_scale(vec3(1.f));
        ImGuizmo::PopID();
    }

    if (useWindow) {
        ImGui::End();
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
