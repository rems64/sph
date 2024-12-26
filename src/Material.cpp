#include "Material.hpp"

Material::Material(handle<ShaderProgram> shader, glm::vec3 albedo)
    : m_shader(shader), m_albedo(albedo) {}
