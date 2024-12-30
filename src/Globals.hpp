#pragma once

#include "glm/ext/vector_float2.hpp"
#include "mem.hpp"
#include "typedefs.hpp"

class Window;
class ResourceManager;
class Transform;
class Camera;

enum AppFlags {
    APP_RUNNING = 1 << 0,
};

struct Timings {
    real_t dt = 0.0;
    real_t time = 0.0;
};

struct Debug {
    uint32_t missed_cells;
};

struct CameraAggregate {
    handle<Camera> camera;
    handle<Transform> camera_transform;
};

struct Cursor {
    bool pressed;
    vec2 position;
};

struct Input {
    struct Cursor cursor;
};

// Globals doesn't own any memory, except static strings
struct Globals {
    const char *application_name = "OpenGL Boilerplate";
    uint32_t application_flags = APP_RUNNING;

    struct Timings t;
    struct Debug debug;
    struct CameraAggregate camera;
    struct Input input;

    handle<Window> window;
    handle<ResourceManager> resource_manager;
};

extern Globals G;
