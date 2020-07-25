uniform float tf_min;
uniform float tf_max;
uniform float tf_log;
uniform sampler2D tf_tex;

bool sample_texture(vec3 tex_curr_pos, inout vec4 curr_color, float tdelta,
                    float t, vec3 dir)
{
    float tm = tf_min;
    float tp = tf_max;
    vec4 tf_sample;

    float map_sample = texture(bitmap_tex, tex_curr_pos).r;
    if (!(map_sample > 0.0)) return false;
 
    float tex_sample = texture(ds_tex, tex_curr_pos).r;
 
    if (tf_log > 0.5) {
       if(tex_sample <= 0.0) return false;
       tex_sample = log(tex_sample);
       tm = log(tm);
       tp = log(tp);
    }

    if(tex_sample < tm) return false;
    if(tex_sample > tp) return false;
    vec2 tex_sample_norm = vec2((tex_sample - tm)/(tp - tm), 0.5);
    tf_sample = texture(tf_tex, tex_sample_norm);
    float dt = length(tdelta * dir);
    float ta = max((1.0f - dt * tf_sample.a), 0.0);

    curr_color = dt * tf_sample + ta * curr_color;
    return true;
}

vec4 cleanup_phase(in vec4 curr_color, in vec3 dir, in float t0, in float t1) 
{
  // It's possible we need to scale, in which case this may be useful:
  // vec3 p0 = v_camera_pos.xyz + dir * t0;
  // vec3 p1 = v_camera_pos.xyz + dir * t1;
  //return vec4(curr_color.r, curr_color.g, curr_color.b,
  //            curr_color.a * length(p1-p0));
  return curr_color;
}

