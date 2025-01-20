#version 330 core

in vec3 fColor;
in float wall;

out vec4 colorOut;

float map(float x, float a, float b, float c, float d) {
    return (x - a) / (b - a) * (d - c) + c;
}

void main() {
    colorOut = vec4(map(wall, 0, 1, 1, 0.3) * fColor, map(wall, 0, 1, 1, 0.5));
}
