#version 330 core

in vec3 fColor;

out vec4 colorOut;

void main() {
    colorOut = vec4(fColor, 1.0);
}
