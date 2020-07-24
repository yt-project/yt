#version 330 core

// Input vertex data, different for all executions of this shader.
in vec3 vertexPosition_modelspace;

// Output data ; will be interpolated for each fragment.
out vec2 UV;

void main()
{
    gl_Position = vec4(vertexPosition_modelspace, 1.0);
    UV = (vertexPosition_modelspace.xy+vec2(1.0,1.0))/2.0;
}
