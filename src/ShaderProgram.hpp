#pragma once

#include "mem.hpp"
#include <glad/glad.h>
#include <glm/ext.hpp>
#include <glm/glm.hpp>
#include <string>

class ShaderProgram {
public:
    // Create the program. A valid OpenGL context must be active.
    ShaderProgram();
    virtual ~ShaderProgram();

    // Generate a minimal shader program, made of one vertex shader and one fragment shader
    static ref<ShaderProgram> genBasicShaderProgram(const std::string &vertexShaderFilename,
                                                    const std::string &fragmentShaderFilename);

    // OpenGL identifier of the program
    GLuint id() const;

    // Loads and compile a shader from a text file, before attaching it to a program
    void loadShader(GLenum type, const std::string &shaderFilename);

    // The main GPU program is ready to be handle streams of polygons
    void link();

    // Activate the program
    void use();

    // Desactivate the current program
    static void stop();

    GLuint getLocation(const std::string &name) const;
    void set(const std::string &name, int value);
    void set(const std::string &name, float value);
    void set(const std::string &name, const glm::vec2 &value);
    void set(const std::string &name, const glm::vec3 &value);
    void set(const std::string &name, const glm::vec4 &value);
    void set(const std::string &name, const glm::mat4 &value);
    void set(const std::string &name, const glm::mat3 &value);

private:
    // Loads the content of an ASCII file in a standard C++ string
    std::string file2String(const std::string &filename);

    GLuint _id = 0;
};
