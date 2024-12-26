#pragma once

#include "typedefs.hpp"
#include <glad/glad.h>
#include <memory>
#include <string>
#include <vector>

#include <glm/ext.hpp>
#include <glm/glm.hpp>

class Mesh {
public:
    virtual ~Mesh();

    const std::vector<vec3> &vertexPositions() const { return _vertexPositions; }
    std::vector<vec3> &vertexPositions() { return _vertexPositions; }

    const std::vector<vec3> &vertexNormals() const { return _vertexNormals; }
    std::vector<vec3> &vertexNormals() { return _vertexNormals; }

    const std::vector<glm::vec2> &vertexTexCoords() const { return _vertexTexCoords; }
    std::vector<glm::vec2> &vertexTexCoords() { return _vertexTexCoords; }

    const std::vector<glm::uvec3> &triangleIndices() const { return _triangleIndices; }
    std::vector<glm::uvec3> &triangleIndices() { return _triangleIndices; }

    // Compute the parameters of a sphere which bounds the mesh
    void computeBoundingSphere(vec3 &center, real_t &radius) const;

    void recomputePerVertexNormals(bool angleBased = false);
    void recomputePerVertexTextureCoordinates();

    void init();
    void render();
    void clear();

    void addOffsetBox(const real_t x,
                      const real_t y,
                      const real_t z,
                      const real_t w,
                      const real_t h,
                      const real_t d);
    void addPlane(const real_t square_half_side = 1.0f);
    void addBox(const real_t w, const real_t h, const real_t d);

private:
    std::vector<vec3> _vertexPositions;
    std::vector<vec3> _vertexNormals;
    std::vector<glm::vec2> _vertexTexCoords;
    std::vector<glm::uvec3> _triangleIndices;

    GLuint _vao = 0;
    GLuint _posVbo = 0;
    GLuint _normalVbo = 0;
    GLuint _texCoordVbo = 0;
    GLuint _ibo = 0;
};

// utility: loader
void loadOFF(const std::string &filename, std::shared_ptr<Mesh> meshPtr);
