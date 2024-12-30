#include "Camera.hpp"

#include "Globals.hpp"
#include "Window.hpp"
#include "glm/ext/matrix_clip_space.hpp"
#include "glm/ext/matrix_common.hpp"
#include "typedefs.hpp"

PerspectiveCamera::PerspectiveCamera(real_t fov, real_t near, real_t far)
    : m_fov(fov), m_near(near), m_far(far), m_pitch(0), m_yaw(0), m_distance(4.f), m_target_pitch(0), m_target_yaw(0), m_target_distance(4.f) {}

const mat4 PerspectiveCamera::projection() const {
    const real_t aspect = G.window.lock()->get_aspect_ratio();
    return glm::perspective(m_fov, aspect, m_near, m_far);
}

real_t &PerspectiveCamera::target_distance() { return m_target_distance; }
real_t &PerspectiveCamera::target_pitch() { return m_target_pitch; }
real_t &PerspectiveCamera::target_yaw() { return m_target_yaw; }

real_t &PerspectiveCamera::distance() { return m_distance; }
real_t &PerspectiveCamera::pitch() { return m_pitch; }
real_t &PerspectiveCamera::yaw() { return m_yaw; }
