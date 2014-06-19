//-----------------------------------------------------------------------------
// Copyright (c) 2014, yt Development Team.
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file COPYING.txt, distributed with this software.
//-----------------------------------------------------------------------------


#include <helper_cuda.h>
#include <helper_math.h>


typedef float3 vec3;
typedef float4 vec4;
typedef float2 vec2;


typedef struct {
    vec4 m[4];
} mat4x4;

typedef struct {
    vec3 m[3];
} mat3x3;

typedef struct {
    vec2 m[2];
} mat2x2;

typedef struct {
    vec3 m[4];
} mat3x4;

typedef struct {
    vec4 m[3];
} mat4x3;

typedef struct {
    vec4 m[2];
} mat4x2;

typedef struct {
    vec2 m[4];
} mat2x4;

typedef mat4 mat4x4;
typedef mat3 mat3x3;
typedef mat2 mat2x2;



// Column Major Always
//   All vectors are k x 1, and cannot be at the front of a multiplication


// -----------------------------------------------------------------------------
// Create Matrices, identity by default
// -----------------------------------------------------------------------------
__device__ __host__ mat4 make_mat4() {
	mat4 tmp;
	tmp.m[0] = make_float4(1.0f, 0.0f, 0.0f, 0.0f);
	tmp.m[1] = make_float4(0.0f, 1.0f, 0.0f, 0.0f);
	tmp.m[2] = make_float4(0.0f, 0.0f, 1.0f, 0.0f);
	tmp.m[3] = make_float4(0.0f, 0.0f, 0.0f, 1.0f);
	return tmp;
}

__device__ __host__ mat3 make_mat3() {
	mat3 tmp;
	tmp.m[0] = make_float3(1.0f, 0.0f, 0.0f);
	tmp.m[1] = make_float3(0.0f, 1.0f, 0.0f);
	tmp.m[2] = make_float3(0.0f, 0.0f, 1.0f);
	return tmp;
}

__device__ __host__ mat2 make_mat2() {
	mat3 tmp;
	tmp.m[0] = make_float2(1.0f, 0.0f);
	tmp.m[1] = make_float2(0.0f, 1.0f);
	return tmp;
}


// -----------------------------------------------------------------------------
// Get Transpose vectors, get column vector
// -----------------------------------------------------------------------------
// 4s
__device__ __host__ vec4 mvec(mat4x4 a, int col) {
	vec4 tmp;
	if     (col == 0) tmp = make_float4(a.m[0].x, a.m[1].x, a.m[2].x, a.m[3].x);
	else if(col == 1) tmp = make_float4(a.m[0].y, a.m[1].y, a.m[2].y, a.m[3].y);
	else if(col == 2) tmp = make_float4(a.m[0].z, a.m[1].z, a.m[2].z, a.m[3].z);
	else if(col == 3) tmp = make_float4(a.m[0].w, a.m[1].w, a.m[2].w, a.m[3].w);
	return tmp;
}

__device__ __host__ vec3 mvec(mat4x3 a, int col) {
	vec3 tmp;
	if     (col == 0) tmp = make_float3(a.m[0].x, a.m[1].x, a.m[2].x);
	else if(col == 1) tmp = make_float3(a.m[0].y, a.m[1].y, a.m[2].y);
	else if(col == 2) tmp = make_float3(a.m[0].z, a.m[1].z, a.m[2].z);
	else if(col == 3) tmp = make_float3(a.m[0].w, a.m[1].w, a.m[2].w);
	return tmp;
}

__device__ __host__ vec2 mvec(mat4x2 a, int col) {
	vec2 tmp;
	if     (col == 0) tmp = make_float2(a.m[0].x, a.m[1].x);
	else if(col == 1) tmp = make_float2(a.m[0].y, a.m[1].y);
	else if(col == 2) tmp = make_float2(a.m[0].z, a.m[1].z);
	else if(col == 3) tmp = make_float2(a.m[0].w, a.m[1].w);
	return tmp;
}
// 3s

__device__ __host__ vec3 mvec(mat3x4 a, int col) {
	vec3 tmp;
	if     (col == 0) tmp = make_float4(a.m[0].x, a.m[1].x, a.m[2].x, a.m[3].x);
	else if(col == 1) tmp = make_float4(a.m[0].y, a.m[1].y, a.m[2].y, a.m[3].y);
	else if(col == 2) tmp = make_float4(a.m[0].z, a.m[1].z, a.m[2].z, a.m[3].z);
	return tmp;
}

__device__ __host__ vec4 mvec(mat3x3 a, int col) {
	vec4 tmp;
	if     (col == 0) tmp = make_float3(a.m[0].x, a.m[1].x, a.m[2].x);
	else if(col == 1) tmp = make_float3(a.m[0].y, a.m[1].y, a.m[2].y);
	else if(col == 2) tmp = make_float3(a.m[0].z, a.m[1].z, a.m[2].z);
	return tmp;
}

__device__ __host__ vec2 mvec(mat3x2 a, int col) {
	vec2 tmp;
	if     (col == 0) tmp = make_float2(a.m[0].x, a.m[1].x);
	else if(col == 1) tmp = make_float2(a.m[0].y, a.m[1].y);
	else if(col == 2) tmp = make_float2(a.m[0].z, a.m[1].z);
	return tmp;
}
// 2s

__device__ __host__ vec2 mvec(mat2x4 a, int col) {
	vec4 tmp;
	if     (col == 0) tmp = make_float4(a.m[0].x, a.m[1].x, a.m[2].x, a.m[3].x);
	else if(col == 1) tmp = make_float4(a.m[0].y, a.m[1].y, a.m[2].y, a.m[3].y);
	return tmp;
}

__device__ __host__ vec3 mvec(mat2x3 a, int col) {
	vec3 tmp;
	if     (col == 0) tmp = make_float3(a.m[0].x, a.m[1].x, a.m[2].x);
	else if(col == 1) tmp = make_float3(a.m[0].y, a.m[1].y, a.m[2].y);
	return tmp;
}

__device__ __host__ vec2 mvec(mat2x2 a, int col) {
	vec2 tmp;
	if     (col == 0) tmp = make_float2(a.m[0].x, a.m[1].x);
	else if(col == 1) tmp = make_float2(a.m[0].y, a.m[1].y);
	return tmp;
}


// -----------------------------------------------------------------------------
// Transpose Matrix
// Cannot transpose anything other than square matrices
// -----------------------------------------------------------------------------
__device__ __host__ mat4 transpose(mat4 a) {
	mat4 tmp;
	for(int i = 0; i < 4; i++)
		tmp.m[i] = mvec(a, i);
	return tmp
}

__device__ __host__ mat3 transpose(mat3 a) {
	mat3 tmp;
	for(int i = 0; i < 3; i++)
		tmp.m[i] = mvec(a, i);
	return tmp
}

__device__ __host__ mat2 transpose(mat2 a) {
	mat2 tmp;
	for(int i = 0; i < 2; i++)
		tmp.m[i] = mvec(a, i);
	return tmp
}

// -----------------------------------------------------------------------------
// Dot Products
// -----------------------------------------------------------------------------
inline __device__ __host__ float dot(vec4 a, vec4 b) {
	return (a.x*b.x + a.y*b.y + a.z*b.z + a.w*b.w);
}

inline __device__ __host__ float dot(vec3 a, vec3 b) {
	return (a.x*b.x + a.y*b.y + a.z*b.z);
}

inline __device__ __host__ float dot(vec2 a, vec2 b) {
	return (a.x*b.x + a.y*b.y);
}

// -----------------------------------------------------------------------------
// Scale Matrix
// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
// Multiply Matrix
// -----------------------------------------------------------------------------
//4s
__device__ __host__ mat4x4 operator*(mat4x4 a, mat4x4 b) {
	mat4x4 tmp;
	for(int i = 0; i < 4; i++){
		tmp.m[i].x = dot(a.m[i], mvec(b, 0));
		tmp.m[i].y = dot(a.m[i], mvec(b, 1));
		tmp.m[i].z = dot(a.m[i], mvec(b, 2));
		tmp.m[i].w = dot(a.m[i], mvec(b, 3));
	}
	return tmp;
}

__device__ __host__ mat4x3 operator*(mat4x4 a, mat4x3 b) {
	mat4x3 tmp;
	for(int i = 0; i < 3; i++){
		tmp.m[i].x = dot(a.m[i], mvec(b, 0));
		tmp.m[i].y = dot(a.m[i], mvec(b, 1));
		tmp.m[i].z = dot(a.m[i], mvec(b, 2));
		tmp.m[i].w = dot(a.m[i], mvec(b, 3));
	}
	return tmp;
}

__device__ __host__ mat4x2 operator*(mat4x4 a, mat4x2 b) {
	mat4x2 tmp;
	for(int i = 0; i < 2; i++){
		tmp.m[i].x = dot(a.m[i], mvec(b, 0));
		tmp.m[i].y = dot(a.m[i], mvec(b, 1));
		tmp.m[i].z = dot(a.m[i], mvec(b, 2));
		tmp.m[i].w = dot(a.m[i], mvec(b, 3));
	}
	return tmp;
}

__device__ __host__ vec4 operator*(mat4x4 a, vec4 b) {
	vec4 tmp;
	tmp.x = dot(a.m[0], b);
	tmp.y = dot(a.m[1], b);
	tmp.z = dot(a.m[2], b);
	tmp.w = dot(a.m[3], b);
	return tmp;
}

__device__ __host__ mat3x4 operator*(mat3x3 a, mat3x4 b);
__device__ __host__ mat3x3 operator*(mat3x3 a, mat3x3 b);
__device__ __host__ mat3x2 operator*(mat3x3 a, mat3x2 b);
__device__ __host__ vec3   operator*(mat3x3 a, vec3   b);

__device__ __host__ mat2x4 operator*(mat2x2 a, mat2x4 b);
__device__ __host__ mat2x3 operator*(mat2x2 a, mat2x3 b);
__device__ __host__ mat2x2 operator*(mat2x2 a, mat2x2 b);
__device__ __host__ vec2   operator*(mat2x2 a, vec2   b);


