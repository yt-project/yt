#version 330 core

in vec4 rgba_values;
in float x_coord;
uniform int channel;

void main() {
    // our *actual* viewport goes from -1 to 1 in x and y, with 0 everywhere for z.
    // our provided viewport describes the bottom left and top right of what we
    // actually what, but in 0..1 coordinates.

    gl_Position = vec4(x_coord * 2.0 - 1.0, rgba_values[channel] * 2.0 - 1.0, 0.0, 1.0);
}
