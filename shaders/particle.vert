#version 330 core
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aNormal;
layout(location = 2) in vec3 aColor;
layout(location = 3) in vec3 aParticlePosition;

out vec3 fColor;

uniform float particle_scale;
uniform mat4 modelMat, viewMat, projMat;

void main() {
    gl_Position = projMat * viewMat * modelMat * vec4(particle_scale*aPos + aParticlePosition, 1.0);
    fColor = aColor;
}
