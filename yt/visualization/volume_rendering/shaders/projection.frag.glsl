bool sample_texture(vec3 tex_curr_pos, inout vec4 curr_color, float tdelta,
                    float t, vec3 dir) {

    vec3 tex_sample = texture(ds_tex, tex_curr_pos).rgb;
    float map_sample = texture(bitmap_tex, tex_curr_pos).r;
    if (map_sample > 0.0) {
        float val = length(tdelta * dir) * tex_sample.r + curr_color.r;
        curr_color = vec4(val, val, val, 1.0);
    }
    return bool(map_sample > 0.0);
}

vec4 cleanup_phase(in vec4 curr_color, in vec3 dir, in float t0, in float t1) 
{
  return vec4(curr_color);
}
