#pragma once

#include "mem.hpp"
#include "typedefs.hpp"

class GLFWwindow;

class Window {
public:
    Window(int width, int height, const char *title);

    bool should_close();

    void swap_buffers() const;
    void poll_events() const;

private:
    GLFWwindow *window;
};

const real_t get_time();

ref<Window> make_window(int width, int height, const char *title);
void init_windowing();
void close_windowing();
