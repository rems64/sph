#pragma once

#include "typedefs.hpp"

class Camera {
public:
    Camera() = default;
    virtual ~Camera() = default;

    virtual const mat4 projection() const = 0;

private:
};

class PerspectiveCamera : public Camera {
public:
    PerspectiveCamera(real_t fov, real_t near, real_t far);
    ~PerspectiveCamera() override = default;

    const mat4 projection() const override;

    real_t &distance();
    real_t &pitch();
    real_t &yaw();
    real_t &target_distance();
    real_t &target_pitch();
    real_t &target_yaw();

private:
    real_t m_fov;
    real_t m_near;
    real_t m_far;

    real_t m_distance;
    real_t m_pitch;
    real_t m_yaw;

    real_t m_target_distance;
    real_t m_target_pitch;
    real_t m_target_yaw;
};
