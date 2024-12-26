#include "BoundingBox.hpp"

#include <algorithm>
#include <limits>

Segment::Segment()
    : min(std::numeric_limits<real_t>::infinity()), max(-std::numeric_limits<real_t>::infinity()) {
}
Segment::Segment(real_t a, real_t b) : min(std::min(a, b)), max(std::max(a, b)) {}

// bool Segment::overlap(Segment &other) { return (min <= other.max && max >= other.min); }
real_t Segment::overlap(Segment &other) {
    // Return the length of the overlapping region
    if (min <= other.max && max >= other.min) {
        return std::min(max, other.max) - std::max(min, other.min);
    }
    return -1.f;
}

void Segment::grow(real_t x) {
    min = std::min(x, min);
    max = std::max(x, max);
}

BoundingBox::BoundingBox() { reset(); }

BoundingBox::BoundingBox(glm::vec<3, real_t> min, glm::vec<3, real_t> max) : min(min), max(max) {}

void BoundingBox::reset() {
    min = glm::vec<3, real_t>(std::numeric_limits<real_t>::infinity());
    max = glm::vec<3, real_t>(-std::numeric_limits<real_t>::infinity());
}

void BoundingBox::grow(glm::vec<3, real_t> v) {
    min.x = std::min(v.x, min.x);
    min.y = std::min(v.y, min.y);
    min.z = std::min(v.z, min.z);

    max.x = std::max(v.x, max.x);
    max.y = std::max(v.y, max.y);
    max.z = std::max(v.z, max.z);
}

bool BoundingBox::intersect(BoundingBox &other) {
    return (min.x <= other.max.x && max.x >= other.min.x) &&
           (min.y <= other.max.y && max.y >= other.min.y) &&
           (min.z <= other.max.z && max.z >= other.min.z);
}
