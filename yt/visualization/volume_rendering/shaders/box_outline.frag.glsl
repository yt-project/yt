#version 330
in vec4 v_model;
in vec3 v_camera_pos;
in vec3 dx;
in vec3 left_edge;
in vec3 right_edge;
flat in mat4 inverse_proj;
flat in mat4 inverse_mvm;
flat in mat4 inverse_pmvm;
out vec4 output_color;

uniform vec4 viewport; // (offset_x, offset_y, 1 / screen_x, 1 / screen_y)
uniform float box_width;
uniform vec3 box_color;
uniform float box_alpha;

//
//  THIS IS CURRENTLY UNUSED
//

void main()
{
    // Obtain screen coordinates
    // https://www.opengl.org/wiki/Compute_eye_space_from_window_space#From_gl_FragCoord
    vec4 ndcPos;
    ndcPos.xy = ((2.0 * gl_FragCoord.xy) - (2.0 * viewport.xy)) / (viewport.zw) - 1;
    ndcPos.z = (2.0 * gl_FragCoord.z - 1.0);
    ndcPos.w = 1.0;

    vec4 clipPos = ndcPos / gl_FragCoord.w;
    vec4 eyePos = inverse_proj * clipPos;
    eyePos /= eyePos.w;

    vec3 ray_position = (inverse_pmvm * clipPos).xyz;

    vec3 dist = min(abs(ray_position - right_edge),
                      abs(ray_position - left_edge));
    
    // We need to be close to more than one edge.

    int count = 0;
    count += int(dist.x < box_width * dx.x);
    count += int(dist.y < box_width * dx.y);
    count += int(dist.z < box_width * dx.z);

    if (count < 2) {
        discard;
    }

    output_color = vec4(box_color, box_alpha);
}
