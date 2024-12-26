#version 330 core

struct Material {
    vec3 albedo;
};
uniform Material material;

uniform vec3 camPos;

in vec3 fPositionModel;
in vec3 fPosition;
in vec3 fNormal;
in vec2 fTexCoord;

out vec4 colorOut;

uniform mat4 modelMat;
uniform mat3 normMat;

void main() {
    vec3 albedo = material.albedo;

    vec3 radiance = vec3(0.0);
    radiance += albedo;

    // colorOut = vec4(radiance, 1.0);
    colorOut = vec4(radiance, 1.0);
}
