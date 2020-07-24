#version 330 core

in vec4 rgba_values;
in float x_coord;
uniform int channel;
uniform vec4 bounds;

void main() {
    // our *actual* viewport goes from -1 to 1 in x and y, with 0 everywhere for z.
    // our provided viewport describes the bottom left and top right of what we
    // actually what, but in 0..1 coordinates.

    float dx = (bounds[1] - bounds[0]);
    float dy = (bounds[3] - bounds[2]);

    gl_Position = vec4( dx * x_coord + bounds[0],
                        dy * rgba_values[channel] + bounds[2],
                        0.0, 1.0);
}
