#version 330 core

in vec2 UV;

out vec4 color;
uniform sampler2D fb_texture;
uniform sampler2D db_texture;

void main(){
   color = texture(fb_texture, UV);
   color.a = 1.0;
   gl_FragDepth = texture(db_texture, UV).r;
}
