#pragma once

#include "imgui.h"
#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"

void ui_init(GLFWwindow *window);
void ui_new_frame();
void ui_render();
void ui_shutdown();