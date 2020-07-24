#version 330 core

// Input vertex data, different for all executions of this shader.
in vec3 vertexPosition_modelspace;

uniform vec4 bounds;

// Output data ; will be interpolated for each fragment.
out vec2 UV;

void main()
{
    // Scale from -1..1 into bounds[0] .. bounds[1]
    UV = (vertexPosition_modelspace.xy + vec2(1,1))/2.0;
    float new_x = UV.x * (bounds[1] - bounds[0]) + bounds[0];
    float new_y = UV.y * (bounds[3] - bounds[2]) + bounds[2];
    gl_Position = vec4(new_x, new_y, vertexPosition_modelspace.z, 1.0);
}
