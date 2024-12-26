#pragma once

#include "ShaderProgram.hpp"
#include "mem.hpp"
#include <glm/ext/vector_float3.hpp>

class Material {
public:
    Material(handle<ShaderProgram> shader, glm::vec3 albedo);

    [[nodiscard]] inline handle<ShaderProgram> shader() const { return m_shader; }
    [[nodiscard]] inline glm::vec3 albedo() const { return m_albedo; }

private:
    handle<ShaderProgram> m_shader;
    glm::vec3 m_albedo;
};
