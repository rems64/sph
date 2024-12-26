#include "Transform.hpp"
#include "glm/gtc/quaternion.hpp"

#include <glm/gtc/matrix_transform.hpp>

mat4 Transform::world_matrix() const {
    mat4 translation = glm::translate(mat4(1.0), position);
    mat4 rotation_matrix = glm::mat4_cast(rotation);
    mat4 scaling = glm::scale(mat4(1.0f), scale);
    return translation * rotation_matrix * scaling;
}
