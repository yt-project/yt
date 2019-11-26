#version 330 core

in vec2 UV;

out vec4 color;

uniform sampler2D fb_texture;
uniform sampler2D db_texture;
uniform sampler1D cmap;
uniform float cmap_min;
uniform float cmap_max;
uniform float cmap_log;

void main(){
   float scaled = texture(fb_texture, UV).x;
   float alpha = texture(fb_texture, UV).a;
   if (alpha == 0.0) discard;
   float cm = cmap_min;
   float cp = cmap_max;

   if (cmap_log > 0.5) {
      scaled = log(scaled);
      cm = log(cm);
      cp = log(cp);
   }
   color = texture(cmap, (scaled - cm) / (cp - cm));
   gl_FragDepth = texture(db_texture, UV).r;
}
