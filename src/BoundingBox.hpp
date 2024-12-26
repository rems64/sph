#pragma once

#include "typedefs.hpp"


struct Segment {
public:
    Segment();
    Segment(real_t a, real_t b);

    real_t overlap(Segment &other);
    void grow(real_t x);

    real_t min;
    real_t max;
};

struct BoundingBox {
public:
    BoundingBox();
    BoundingBox(vec3 min, vec3 max);

    void reset();
    void grow(vec3 v);
    bool intersect(BoundingBox &other);

    vec3 min;
    vec3 max;
};
