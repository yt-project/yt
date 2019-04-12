#version 330
in vec4 quad_vertex; // The location of the vertex in model space
uniform float x_offset;
uniform float y_offset;
uniform vec4 viewport;
uniform float scale;
uniform float y_origin;
uniform float x_origin;

out vec2 UV;

void main()
{
    // What we're passed in is all in 64/th of pixels.  Our quad_vertex 
    // runs from 0 .. width, 0 .. height.  Our x_offset is also in 64/th 
    // of pixels, and same for y_offset.
    // Note also that our offset is from the top.  So what we need to compute
    // is the dx and dy, for coordinates per pixel.
    float dx = scale * 2.0 / (viewport.z * 64);
    float dy = scale * 2.0 / (viewport.w * 64);
    // We set as our origin (-0.9, 0.9) (for the bottom left of the
    // character) in clip space coordinates.  So we offset from that.
    // To shift into clip space, we figure out how big our object needs to be
    // clip space.  We take our clip per pixel, multiply by the pixels, and
    // then set the origin to be 0.0 (clip) and offset by - y_offset * dy.
    gl_Position = vec4( dx * (quad_vertex.x + x_offset) + x_origin,
                        dy * (quad_vertex.y + y_offset) + y_origin,
                        0.0, 1.0); 
    UV = vec2(quad_vertex.wz);
    
}
