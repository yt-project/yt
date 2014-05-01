//-----------------------------------------------------------------------------
// Copyright (c) 2014, yt Development Team.
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file COPYING.txt, distributed with this software.
//-----------------------------------------------------------------------------

#include <helpers/cuda.cu>

//A 3d texture holding scalar values to be volume rendered
texture<float , cudaTextureType3D, cudaReadModeElementType> volume;
//A 1d texture holding color transfer function values to act on raycaster results
texture<float4, cudaTextureType1D, cudaReadModeElementType> transfer;


extern "C"
__global__ void front_to_back(uint *buffer, float *modelviewmatrix, int buffer_w,
                            int buffer_h, int max_steps, float density,
                            float offset, float scale,
                            float brightness, float step_size
                            ) {
  //Rays will terminate when opacity_threshold is reached
  const float opacity_threshold = 0.95f; 

  //We don't use the far and near clipping information from the modelviewmatrix
  float3x4 matrix = mat_array_to_3x4(modelviewmatrix); 

  //The X,Y coordinate of the surface (image) array
  int x = BLOCK_X;
  int y = BLOCK_Y;

  if ((x >= buffer_w) || (y >= buffer_h)) return;

  float u = COORD_STD(x, buffer_w);
  float v = COORD_STD(y, buffer_h);

  // calculate eye ray in world space
  Ray eye_ray = make_ray_from_eye(matrix, u, v);

  // find intersection with box
  float tnear, tfar;
  int hit = intersect_std_box(eye_ray, &tnear, &tfar);

  // If the ray misses the box set the coloro to black
  if (!hit) {
    buffer[y*buffer_w + x] = 0;
    return;
  }

  if (tnear < 0.0f) tnear = 0.0f;     // clamp to near plane

  // march along ray from front to back, accumulating color
  float4 sum = make_float4(0.0f);
  float t = tnear;
  float3 pos = eye_ray.o + eye_ray.d*tnear;
  float3 step = eye_ray.d*step_size;

  // Cast Rays
  float4 col =  make_float4(0.0,0.0,0.0,0.0);
  float top = (1.0/scale) + offset;
  for (int i=0; i<=max_steps; i++) {
    float sample = tex3D(volume, pos.x*0.5f+0.5f, pos.y*0.5f+0.5f, pos.z*0.5f+0.5f);
    

    if (sample > offset) {
        col = tex1D(transfer, (sample-offset)*scale);
    }
    else {
        col = make_float4(0.0,0.0,0.0,0.0);
    }

    col.w *= density;

    // pre-multiply alpha
    col.x *= col.w;
    col.y *= col.w;
    col.z *= col.w;

    // "over" operator for front-to-back blending
    sum = sum + col*(1.0f - sum.w);

    // exit early if opaque
    if (sum.w > opacity_threshold)
      break;

    // Increment step size and position
    t += step_size;
    pos += step;

    // If ray has cast too far, end
    if (t > tfar) break;
  }

  sum *= brightness;

  buffer[SCREEN_XY(x,y,buffer_w)] = rgbaFloatToInt(sum);
}

