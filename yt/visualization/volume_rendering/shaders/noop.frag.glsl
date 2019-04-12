#version 330 core

in vec2 UV;

out vec4 color;

uniform sampler2D fb_texture;
uniform sampler1D cmap;
uniform float min_val;
uniform float scale;
uniform float cmap_min;
uniform float cmap_max;
uniform float cmap_log;

void main(){
   color = texture(fb_texture, UV);
}
