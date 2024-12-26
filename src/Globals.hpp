#pragma once

#include "mem.hpp"
#include "typedefs.hpp"

class Window;
class ResourceManager;

enum AppFlags {
    APP_RUNNING = 1 << 0,
};

struct Timings {
    real_t dt = 0.0;
    real_t time = 0.0;
};

// Globals doesn't own any memory, except static strings
struct Globals {
    const char *application_name = "OpenGL Boilerplate";
    uint32_t application_flags = APP_RUNNING;

    struct Timings t;

    handle<Window> window;
    handle<ResourceManager> resource_manager;
};

extern Globals G;
