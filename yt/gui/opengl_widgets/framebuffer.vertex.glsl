/* Vertex shader to calculate the ray to cast through the volume data.
 * The result is passed to the fragment shader using a varying.
 *
 * This file is part of Visvis.
 * Copyright 2009 Almar Klein
 */

void main()
{    
    gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
    gl_TexCoord[0].xyz = gl_MultiTexCoord0.xyz;
}
