#version 330
in vec4 v_model;
in vec3 v_camera_pos;
flat in mat4 inverse_proj;
flat in mat4 inverse_view;
out vec4 output_color;

uniform sampler3D ds_tex;
//layout (binding = 1) uniform sampler2D depth_tex;
uniform vec3 dx;
uniform vec3 left_edge;
uniform vec3 right_edge;
uniform vec4 viewport; // (offset_x, offset_y, 1 / screen_x, 1 / screen_y)

bool within_bb(vec3 pos)
{
    bvec3 left =  greaterThanEqual(pos, left_edge);
    bvec3 right = lessThanEqual(pos, right_edge);
    return all(left) && all(right);
}

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

    vec4 world_location = inverse_view * eyePos;

    float step_size = length(right_edge - left_edge) / dx.x / 100000.0;
    vec3 dir = normalize(world_location.xyz - v_camera_pos.xyz);
    vec3 curr_color = vec3(0.0);

    vec3 ray_position = world_location.xyz;

    vec3 tex_curr_pos = vec3(0.0);
    vec3 range = right_edge - left_edge;
    bool ray_in_bb = true;
    while (ray_in_bb) {
        tex_curr_pos = (ray_position - left_edge)/range;

        vec3 tex_sample = texture(ds_tex, tex_curr_pos).rgb;
        if (length(curr_color) < length(tex_sample)) {
            curr_color = tex_sample;
        }

        ray_position += dir * step_size;
        ray_in_bb = within_bb(ray_position);
    }
    output_color = vec4(curr_color.rrr, 1.0);
}
