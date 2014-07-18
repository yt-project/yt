//-----------------------------------------------------------------------------
// Copyright (c) 2014, yt Development Team.
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file COPYING.txt, distributed with this software.
//-----------------------------------------------------------------------------


#include <helper_cuda.h>
#include <helper_math.h>

#define BLOCK_X (blockIdx.x*blockDim.x + threadIdx.x)
#define BLOCK_Y (blockIdx.y*blockDim.y + threadIdx.y)
#define SCREEN_XY(a, b, w) (b*w+a)
#define COORD_STD(x, w) ((x / (float) w)*2.0f-1.0f)

#define COLOR_BLACK make_float3(0.0f, 0.0f, 0.0f)

typedef struct {
    float4 m[3];
} float3x4;

typedef struct {
    float4 m[4];
} float4x4;


struct Ray {
    float3 o;   // origin
    float3 d;   // direction
};



__device__ float4   mul(const float3x4 &M, const float4 &v);
__device__ float3   mul(const float3x4 &M, const float3 &v);
__device__ float3x4 mat_array_to_3x4(float *m);
__device__ Ray      make_ray_from_eye(float3x4 m, float u, float v);
__device__ int      intersect_box(Ray r, float3 boxmin, float3 boxmax,
                            float *tnear, float *tfar);
__device__ int      intersect_std_box(Ray r, float *tnear, float *tfar);

__device__ float    mag(const float3 &v);
__device__ float    avg(const float3 &v);
__device__ uint rgbaFloatToInt(float4 rgba);
__device__ uint rgbFloatToInt(float3 rgb);
__device__ float4   intToRGBAFloat(uint irgba);
__device__ float3   intToRGBFloat(uint irgba);



__device__ uint rgbFloatToInt(float3 rgba)
{
    rgba.x = __saturatef(rgba.x);   // clamp to [0.0, 1.0]
    rgba.y = __saturatef(rgba.y);
    rgba.z = __saturatef(rgba.z);
    return (uint(rgba.z*255)<<16) | (uint(rgba.y*255)<<8) | uint(rgba.x*255);
}

__device__ uint rgbaFloatToInt(float4 rgba)
{
    rgba.x = __saturatef(rgba.x);   // clamp to [0.0, 1.0]
    rgba.y = __saturatef(rgba.y);
    rgba.z = __saturatef(rgba.z);
    rgba.w = __saturatef(rgba.w);
    return (uint(rgba.w*255)<<24) | (uint(rgba.z*255)<<16) | (uint(rgba.y*255)<<8) | uint(rgba.x*255);
}

__device__ float4 intToRGBAFloat(uint irgba)
{
    float4 rgba;
    rgba.w = (irgba >> 24 & 0xFF)/255;
    rgba.z = (irgba >> 16 & 0xFF)/255;
    rgba.y = (irgba >> 8 & 0xFF)/255;
    rgba.x = (irgba & 0xFF)/255;
    return rgba;
}

__device__ float3 intToRGBFloat(uint irgba)
{
    float3 rgb;
    rgb.z = (irgba >> 16 & 0xFF)/255;
    rgb.y = (irgba >> 8 & 0xFF)/255;
    rgb.x = (irgba & 0xFF)/255;
    return rgb;
}


__device__ float3x4 mat_array_to_3x4(float *m) {
	float3x4 matrix;
  matrix.m[0] = make_float4(m[0], m[1], m[2], m[3]);
  matrix.m[1] = make_float4(m[4], m[5], m[6], m[7]);
  matrix.m[2] = make_float4(m[8], m[9], m[10], m[11]);
  return matrix;
}

__device__ Ray make_ray_from_eye(float3x4 m, float u, float v) {
  Ray eyeRay;
  eyeRay.o = make_float3(mul(m, make_float4(0.0f, 0.0f, 0.0f, 1.0f)));
  eyeRay.d = normalize(make_float3(u, v, -2.0f));
  eyeRay.d = mul(m, eyeRay.d);
  return eyeRay;
}

__device__
int intersect_box(Ray r, float3 boxmin, float3 boxmax, float *tnear, float *tfar)
{
    // compute intersection of ray with all six bbox planes
    float3 invR = make_float3(1.0f) / r.d;
    float3 tbot = invR * (boxmin - r.o);
    float3 ttop = invR * (boxmax - r.o);

    // re-order intersections to find smallest and largest on each axis
    float3 tmin = fminf(ttop, tbot);
    float3 tmax = fmaxf(ttop, tbot);

    // find the largest tmin and the smallest tmax
    float largest_tmin = fmaxf(fmaxf(tmin.x, tmin.y), fmaxf(tmin.x, tmin.z));
    float smallest_tmax = fminf(fminf(tmax.x, tmax.y), fminf(tmax.x, tmax.z));

    *tnear = largest_tmin;
    *tfar = smallest_tmax;

    return smallest_tmax > largest_tmin;
}

__device__
int intersect_std_box(Ray r, float *tnear, float *tfar)
{
    const float3 boxmin = make_float3(-1.0f, -1.0f, -1.0f);
    const float3 boxmax = make_float3( 1.0f,  1.0f,  1.0f);

    // compute intersection of ray with all six bbox planes
    float3 invR = make_float3(1.0f) / r.d;
    float3 tbot = invR * (boxmin - r.o);
    float3 ttop = invR * (boxmax - r.o);

    // re-order intersections to find smallest and largest on each axis
    float3 tmin = fminf(ttop, tbot);
    float3 tmax = fmaxf(ttop, tbot);

    // find the largest tmin and the smallest tmax
    float largest_tmin = fmaxf(fmaxf(tmin.x, tmin.y), fmaxf(tmin.x, tmin.z));
    float smallest_tmax = fminf(fminf(tmax.x, tmax.y), fminf(tmax.x, tmax.z));

    *tnear = largest_tmin;
    *tfar = smallest_tmax;

    return smallest_tmax > largest_tmin;
}

__device__
float4 mul(const float3x4 &M, const float4 &v)
{
    float4 r;
    r.x = dot(v, M.m[0]);
    r.y = dot(v, M.m[1]);
    r.z = dot(v, M.m[2]);
    r.w = 1.0f;
    return r;
}

// transform vector by matrix with translation

__device__
float3 mul(const float3x4 &M, const float3 &v)
{
    float3 r;
    r.x = dot(v, make_float3(M.m[0]));
    r.y = dot(v, make_float3(M.m[1]));
    r.z = dot(v, make_float3(M.m[2]));
    return r;
}


__device__
float mag(const float3 &v) {
	return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}


__device__
float avg(const float3 &v) {
	return (v.x + v.y + v.z) / 3.0;
}



