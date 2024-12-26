#include "Window.hpp"

// clang-format off
#include <glad/glad.h>
#include <GLFW/glfw3.h>
// clang-format on

Window::Window(int width, int height, const char *title) {
    window = glfwCreateWindow(width, height, title, nullptr, nullptr);
    glfwMakeContextCurrent(window);
    gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);
}

bool Window::should_close() { return glfwWindowShouldClose(window); }

void Window::swap_buffers() const { glfwSwapBuffers(window); }

void Window::poll_events() const { glfwPollEvents(); }

const real_t get_time() { return static_cast<real_t>(glfwGetTime()); }

void init_windowing() {
    glfwInit();

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);
}

void close_windowing() { glfwTerminate(); }

ref<Window> make_window(int width, int height, const char *title) {
    return make_ref<Window>(width, height, title);
}
