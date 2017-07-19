#version 330 core

in vec2 UV;

out vec4 color;
uniform sampler2D fb_texture;
uniform sampler2D db_texture;

void main(){
   float val = texture(fb_texture, UV).r;
   gl_FragDepth = 0.0;
   if(val == 0) discard;
   color = vec4(val);
}

