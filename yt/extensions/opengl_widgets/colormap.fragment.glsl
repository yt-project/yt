uniform sampler1D colormap;
uniform sampler2D buffer;
void main()
{
    vec2 pos = vec3(gl_TexCoord[0]).xy;
    vec4 value = texture2D(buffer, pos);
    //gl_FragColor = texture1D(colormap, pos[1]);
    gl_FragColor = texture1D(colormap, value[0]);
    //gl_FragColor = vec4(1.0, 0.5, 0.3, 1.0);
    //gl_FragColor = vec4(value[0], value[1], value[2], 1.0);
    //gl_FragColor = vec4(pos[1], 0.0, 0.0, 1.0);
}
