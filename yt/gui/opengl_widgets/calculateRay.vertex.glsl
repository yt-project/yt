/* Vertex shader to calculate the ray to cast through the volume data.
 * The result is passed to the fragment shader using a varying.
 *
 * This file is part of Visvis.
 * Copyright 2009 Almar Klein
 */

// the dimensions of the data, to determine stepsize
uniform vec3 shape;

// ratio to tune the number of steps
uniform float stepRatio;

// varyings to pass to fragment shader
varying vec3 ray;

void main()
{    
    
    // First of all, set position.
    // (We need to do this because this shader replaces the original shader.)
    gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
    
    // Store texture coordinate (also a default thing).
    gl_TexCoord[0].xyz = gl_MultiTexCoord0.xyz;
    
    // Calculate the scaling of the modelview matrix so we can correct
    // for axes.daspect and scale transforms of the wobject (in case
    // of anisotropic data).
    // We go from world coordinates to eye coordinates.
    vec4 p0 = gl_ModelViewMatrix * vec4(0.0,0.0,0.0,1.0);
    vec4 px = gl_ModelViewMatrix * vec4(1.0,0.0,0.0,1.0);
    vec4 py = gl_ModelViewMatrix * vec4(0.0,1.0,0.0,1.0);
    vec4 pz = gl_ModelViewMatrix * vec4(0.0,0.0,1.0,1.0);
    float sx = length(p0.xyz - px.xyz);
    float sy = length(p0.xyz - py.xyz);
    float sz = length(p0.xyz - pz.xyz);
    
    // Create a (diagonal) matrix to correct for the scaling
    mat4 Ms = mat4(0.0);
    Ms[0][0] = 1.0/(sx*sx);
    Ms[1][1] = 1.0/(sy*sy);
    Ms[2][2] = 1.0/(sz*sz);
    Ms[3][3] = 1.0;
    
    // Calculate ray direction. By correcting for the scaling, the ray is
    // expressed in textute coordinates.
    // We go from eye coordinates to world/texture coordinates.
    vec4 p1 = vec4(0.0, 0.0, 0.0, 1.0) * gl_ModelViewProjectionMatrix * Ms;
    vec4 p2 = vec4(0.0, 0.0, 1.0, 1.0) * gl_ModelViewProjectionMatrix * Ms;
    ray = (p2.xyz/p2[3]) - (p1.xyz/p1[3]);
    
    // Normalize ray to unit length.    
    ray = normalize(ray);
    
    // Make the ray represent the length of a single voxel.
    ray = ray / shape;
    ray = ray * 0.58; // 1 over root of three = 0.577
    
    // Scale ray to take smaller steps.
    ray = ray / stepRatio;
    
}
