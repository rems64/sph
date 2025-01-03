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

const vec3 sun_dir = normalize(vec3(0.4, 0.2, -0.7));

void main() {
    vec3 albedo = material.albedo;

    vec3 radiance = vec3(0.0);
    radiance += albedo;

    float ambiant = 0.1;
    float luminosity = max(0., dot(normalize(fNormal), -sun_dir)) / (1+ambiant) + ambiant;

    // colorOut = vec4(radiance, 1.0);
    colorOut = vec4(radiance * luminosity, 1.0);
}
