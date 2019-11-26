#version 330
in vec4 model_vertex; // The location of the vertex in model space
in vec3 in_dx;
in vec3 in_left_edge;
in vec3 in_right_edge;
out vec4 v_model;
out vec3 v_camera_pos;
flat out mat4 inverse_proj;
flat out mat4 inverse_mvm;
flat out mat4 inverse_pmvm;
out vec3 dx;
out vec3 left_edge;
out vec3 right_edge;

//// Uniforms
uniform vec3 camera_pos;
uniform mat4 modelview;
uniform mat4 projection;

uniform float near_plane;
uniform float far_plane;

void main()
{
    v_model = model_vertex;
    v_camera_pos = camera_pos;
    inverse_proj = inverse(projection);
    // inverse model-view-matrix
    inverse_mvm = inverse(modelview);
    inverse_pmvm = inverse(projection * modelview);
    gl_Position = projection * modelview * model_vertex;
    dx = vec3(in_dx);
    left_edge = vec3(in_left_edge);
    right_edge = vec3(in_right_edge);

}
