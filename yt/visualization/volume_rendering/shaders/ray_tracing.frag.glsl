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

uniform vec3 camera_pos;
uniform mat4 modelview;
uniform mat4 projection;

uniform sampler3D ds_tex;
uniform sampler3D bitmap_tex;
//layout (binding = 1) uniform sampler2D depth_tex;
uniform vec4 viewport; // (offset_x, offset_y, 1 / screen_x, 1 / screen_y)

bool within_bb(vec3 pos)
{
    bvec3 left =  greaterThanEqual(pos, left_edge);
    bvec3 right = lessThanEqual(pos, right_edge);
    return all(left) && all(right);
}

bool sample_texture(vec3 tex_curr_pos, inout vec4 curr_color, float tdelta,
                    float t, vec3 dir);
vec4 cleanup_phase(in vec4 curr_color, in vec3 dir, in float t0, in float t1);

// This main() function will call a function called sample_texture at every
// step along the ray.  It must be of the form
//   void (vec3 tex_curr_pos, inout vec4 curr_color, float tdelta, float t,
//         vec3 direction);

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

    // Five samples
    vec3 step_size = dx/5.0;
    vec3 dir = normalize(v_camera_pos.xyz - ray_position);
    dir = max(abs(dir), 0.0001) * sign(dir);
    vec4 curr_color = vec4(0.0);

    // We need to figure out where the ray intersects the box, if it intersects the box.
    // This will help solve the left/right edge issues.

    vec3 idir = 1.0/dir;
    // These 't' prefixes actually mean 'parameter', as we use in grid_traversal.pyx.

    vec3 tl = (left_edge - ray_position)*idir;
    vec3 tr = (right_edge - ray_position)*idir;

    vec3 tmin = min(tl, tr);
    vec3 tmax = max(tl, tr); 

    vec2 temp_t;

    // biggest tmin
    temp_t = max(tmin.xx, tmin.yz);
    float t0 = max(temp_t.x, temp_t.y);

    // smallest tmax
    temp_t = min(tmax.xx, tmax.yz);
    float t1 = min(temp_t.x, temp_t.y);
    if (t1 <= t0) discard;

    // Some more discussion of this here:
    //  http://prideout.net/blog/?p=64

    t0 = max(t0, 0.0);
    t1 = max(t1, 0.0);

    vec3 p0 = v_camera_pos.xyz + dir * t0;
    vec3 p1 = v_camera_pos.xyz + dir * t1;

    vec3 dxidir = dx * abs(idir) / 10.0;

    temp_t = min(dxidir.xx, dxidir.yz);

    float tdelta = min(temp_t.x, temp_t.y);
    float t = t0;

    vec3 range = (right_edge + dx/2.0) - (left_edge - dx/2.0);
    vec3 nzones = range / dx;
    vec3 ndx = 1.0/nzones;

    temp_t = max(nzones.xx, nzones.yz);

    tdelta = max((t1 - t0)/(5.0*max(temp_t.x, temp_t.y)), tdelta);

    vec3 tex_curr_pos = vec3(0.0);

    vec3 step = normalize(p1 - p0) * step_size;
    bool sampled;
    bool ever_sampled = false;
    vec3 last_sampled;

    vec4 v_clip_coord;
    float f_ndc_depth;
    float depth = 1.0;


    while(t <= t1) {
        tex_curr_pos = (ray_position - left_edge) / range;  // Scale from 0 .. 1
        // But, we actually need it to be 0 + normalized dx/2 to 1 - normalized dx/2
        tex_curr_pos = (tex_curr_pos * (1.0 - ndx)) + ndx/2.0;

        sampled = sample_texture(tex_curr_pos, curr_color, tdelta, t, dir);

        if (sampled) {
            ever_sampled = true;
            v_clip_coord = projection * modelview * vec4(ray_position, 1.0);
            f_ndc_depth = v_clip_coord.z / v_clip_coord.w;
            depth = min(depth, (1.0 - 0.0) * 0.5 * f_ndc_depth + (1.0 + 0.0) * 0.5);

        }

        t += tdelta;
        ray_position += tdelta * dir;

    }

    output_color = cleanup_phase(curr_color, dir, t0, t1);

    if (ever_sampled) {
        gl_FragDepth = depth;
    }
}
