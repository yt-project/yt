#version 330 core

uniform int channel;
out vec4 color;

void main() {
    color = vec4(0.0, 0.0, 0.0, 1.0);
    color[channel] = 1.0;
    color.rgb = color.rgb + vec3(float(channel == 3));
}
