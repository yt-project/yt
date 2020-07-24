#version 330 core

in vec2 UV;

out vec4 color;
uniform sampler1D cm_texture;

void main(){
   color = texture(cm_texture, UV.x);
   color.a = 1.0;
}
