#pragma once

#include "typedefs.hpp"
#include <cmath>

static real_t base_kernel(real_t dst, real_t radius) {
    if (dst < radius) {
        float scale = 15 / (2 * M_PI * pow(radius, 5));
        // float scale = 3.f / (2 * M_PI * M_PI * powf(radius, 3));
        float v = radius - dst;
        return v * v * scale;
    }
    return 0;
}

static real_t base_kernel_derivative(real_t dst, real_t radius) {
    if (dst <= radius) {
        float scale = 15 / (powf(radius, 5) * M_PI);
        // float scale = 3.f / (M_PI * M_PI * powf(radius, 3));
        float v = radius - dst;
        return -v * scale;
    }
    return 0;
}

static real_t base_near_kernel(real_t dst, real_t radius) {
    if (dst < radius) {
        // float scale = 15 / (M_PI * pow(radius, 6));
        float scale = 2.f / (M_PI * M_PI * powf(radius, 4));
        float v = radius - dst;
        return v * v * v * scale;
    }
    return 0;
}

static real_t base_near_kernel_derivative(real_t dst, real_t radius) {
    if (dst <= radius) {
        float scale = 2.f / (M_PI * M_PI * powf(radius, 4));
        float v = radius - dst;
        return -3.f * v * v * scale;
    }
    return 0;
}

static real_t cubic_spline_kernel(real_t r, real_t h) {
    const real_t q = r / h;
    const real_t alpha3 = 3.f / (2.f * M_PI * h * h * h);
    if (q >= 0 && q < 1) {
        return alpha3 * (2.f / 3.f - q * q + 0.5f * q * q * q);
    } else if (q >= 1 && q < 2) {
        return alpha3 * (1.f / 6.f) * (2.f - q) * (2.f - q) * (2.f - q);
    } else {
        return 0;
    }
}

static real_t cubic_spline_kernel_derivative(real_t r, real_t h) {
    const real_t q = r / h;
    const real_t alpha3 = 3.f / (2.f * M_PI * h * h * h);
    if (q >= 0 && q < 1) {
        return alpha3 * (-2.f * q / h + 3.f * 0.5f * q * q / h / h);
    } else if (q >= 1 && q < 2) {
        return alpha3 * (-1.f / 2.f) * (2.f - q / h) * (2.f - q / h);
    } else {
        return 0;
    }
}
