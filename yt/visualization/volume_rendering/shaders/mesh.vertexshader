#version 330 core

in vec3 vertexPosition_modelspace;
in float vertexData;
out float fragmentData;
uniform mat4 model_to_clip;
void main()
{
    gl_Position = model_to_clip * vec4(vertexPosition_modelspace, 1);
    fragmentData = vertexData;
}
