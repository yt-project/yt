/***********************************************************************
An attempt at putting the ray-casting operation into CUDA
An attempt at putting the ray-casting operation into CUDA

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2009 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/

//extern __shared__ float array[];

#define NUM_SAMPLES 5
#define VINDEX(A,B,C) data[((((A)+ci[0])*(ds[1]+1)+((B)+ci[1]))*(ds[2]+1)+ci[2]+(C))]

#define fmin(A, B) ( (A < B) ? A : B )
#define fmax(A, B) ( (A > B) ? A : B )
#define fclip(A, B, C) ( fmax( fmin(A, C), B) )

struct transfer_function
{
    float *vs[4];
    float dbin;
    float bounds[2];
};

struct grid
{
    float left_edge[3];
    float right_edge[3];
    float dds[3];
    int dims[3];
    float *data;
};

__device__ float interpolate(float *data, int ds[3], int *ci, float *dp)
{
    int i;
    float dv, dm[3];
    for(i=0;i<3;i++)dm[i] = (1.0 - dp[i]);
    dv  = 0.0;
    dv += VINDEX(0,0,0) * (dm[0]*dm[1]*dm[2]);
    dv += VINDEX(0,0,1) * (dm[0]*dm[1]*dp[2]);
    dv += VINDEX(0,1,0) * (dm[0]*dp[1]*dm[2]);
    dv += VINDEX(0,1,1) * (dm[0]*dp[1]*dp[2]);
    dv += VINDEX(1,0,0) * (dp[0]*dm[1]*dm[2]);
    dv += VINDEX(1,0,1) * (dp[0]*dm[1]*dp[2]);
    dv += VINDEX(1,1,0) * (dp[0]*dp[1]*dm[2]);
    dv += VINDEX(1,1,1) * (dp[0]*dp[1]*dp[2]);
    return dv;
}

__device__ void eval_transfer(float dt, float dv, float *rgba,
                               transfer_function &tf)
{
    int i, bin_id;
    float temp, bv, dy, dd, ta;

    bin_id = (int) ((dv - tf.bounds[0]) / tf.dbin);
    bv = tf.vs[3][bin_id  ];
    dy = tf.vs[3][bin_id+1] - bv;
    dd = dv - (tf.bounds[0] + bin_id*tf.dbin);
    temp = bv+dd*(dy/tf.dbin);
    ta = temp;
    for (i = 0; i < 3; i++)
    {
        bv = tf.vs[i][bin_id  ];
        dy = tf.vs[i][bin_id+1];
        dd = dv - (tf.bounds[0] + bin_id*tf.dbin);
        temp = bv+dd*(dy/tf.dbin);
        rgba[i] += (1.0 - rgba[3])*ta*temp*dt;
    }
    //rgba[3] += (1.0 - rgba[3])*ta*dt;
}

__device__ void sample_values(float v_pos[3], float v_dir[3],
                   float enter_t, float exit_t, int ci[3], float rgba[4],
                   transfer_function &tf, grid &tg)
{
    float cp[3], dp[3], dt, t, dv;
    int dti, i;
    dt = (exit_t - enter_t) / (NUM_SAMPLES-1);
    for (dti = 0; dti < NUM_SAMPLES - 1; dti++)
    {
        t = enter_t + dt*dti;
        for (i = 0; i < 3; i++)
        {
            cp[i] = v_pos[i] + t * v_dir[i];
            dp[i] = fclip(fmod(cp[i], tg.dds[i])/tg.dds[i], 0.0, 1.0);
        }
        dv = interpolate(tg.data, tg.dims, ci, dp);
        eval_transfer(dt, dv, rgba, tf);
    }
    rgba[3] += dt;
}
                   

/* We need to know several things if we want to ray cast through a grid.
   We need the grid spatial information, as well as its values.  We also need
   the transfer function, which defines what our image will look like. */

__global__ void ray_cast(int ngrids,
                         float *grid_data,
                         int *dims,
                         float *left_edge,
                         float *right_edge,
                         float *tf_r,
                         float *tf_g,
                         float *tf_b,
                         float *tf_a,
                         float *tf_bounds,
                         float *v_dir,
                         float *av_pos,
                         float *image_r,
                         float *image_g,
                         float *image_b,
                         float *image_a)
{

    int cur_ind[3], step[3], x, y, i, direction;
    float intersect_t = 1.0;
    float tmax[3];
    float ti, temp_x, temp_y;

    int offset;

    transfer_function tf;
    tf.vs[0] = tf_r;
    tf.vs[1] = tf_g;
    tf.vs[2] = tf_b;
    tf.vs[3] = tf_a;
    tf.bounds[0] = tf_bounds[0]; tf.bounds[1] = tf_bounds[1];
    tf.dbin = (tf.bounds[1] - tf.bounds[0])/255.0;

    /* Set up the grid, just for convenience */
    grid tg;
    int grid_i;

    int tidx = (blockDim.x * gridDim.x) * (
                    blockDim.y * blockIdx.y + threadIdx.y)
             + (blockDim.x * blockIdx.x + threadIdx.x);

    float rgba[4];
    rgba[0] = image_r[tidx];
    rgba[1] = image_g[tidx];
    rgba[2] = image_b[tidx];
    rgba[3] = image_a[tidx];

    float v_pos[3];
    v_pos[0] = av_pos[tidx + 0];
    v_pos[1] = av_pos[tidx + 1];
    v_pos[2] = av_pos[tidx + 2];

    tg.data = grid_data;

    for(grid_i = 0; grid_i < 500; grid_i++) {

        tg.dims[0] = dims[grid_i + 0];
        tg.dims[1] = dims[grid_i + 1];
        tg.dims[2] = dims[grid_i + 2];

        tg.left_edge[0] = left_edge[grid_i + 0];
        tg.left_edge[1] = left_edge[grid_i + 1];
        tg.left_edge[2] = left_edge[grid_i + 2];

        tg.right_edge[0] = right_edge[grid_i + 0];
        tg.right_edge[1] = right_edge[grid_i + 1];
        tg.right_edge[2] = right_edge[grid_i + 2];

        tg.dds[0] = (tg.right_edge[0] - tg.left_edge[0])/tg.dims[0];
        tg.dds[1] = (tg.right_edge[1] - tg.left_edge[1])/tg.dims[1];
        tg.dds[2] = (tg.right_edge[2] - tg.left_edge[2])/tg.dims[2];

        /* We integrate our ray */

        for (i = 0; i < 3; i++)
        {
            step[i] = ((v_dir[i] < 0) ? -1 : 1);
            x = (i + 1) % 3;
            y = (i + 2) % 3;

            ti = (tg.left_edge[i] - v_pos[i])/v_dir[i];
            temp_x = (v_pos[i] + ti*v_dir[x]);
            temp_y = (v_pos[i] + ti*v_dir[y]);

            if( (tg.left_edge[x] <= temp_x) &&
                (temp_x <= tg.right_edge[x]) &&
                (tg.left_edge[y] <= temp_y) &&
                (temp_y <= tg.right_edge[y]) &&
                (0.0 <= ti) && (ti < intersect_t) ) intersect_t = ti;

            ti = (tg.right_edge[i] - v_pos[i])/v_dir[i];
            temp_x = (v_pos[x] + ti*v_dir[x]);
            temp_y = (v_pos[y] + ti*v_dir[y]);

            if( (tg.left_edge[x] <= temp_x) &&
                (temp_x <= tg.right_edge[x]) &&
                (tg.left_edge[y] <= temp_y) &&
                (temp_y <= tg.right_edge[y]) &&
                (0.0 <= ti) && (ti < intersect_t) ) intersect_t = ti;

        }

        if ( (tg.left_edge[0] <= v_pos[0]) &&
             (v_pos[0] <= tg.right_edge[0]) &&
             (tg.left_edge[0] <= v_pos[0]) &&
             (v_pos[0] <= tg.right_edge[0]) &&
             (tg.left_edge[0] <= v_pos[0]) &&
             (v_pos[0] <= tg.right_edge[0]))
                intersect_t = 0.0;

        int skip = 0;
        if((intersect_t < 0) || (intersect_t > 1.0)) skip=1;

        for (i = 0; i < 3;  i++)
        {
            cur_ind[i] = (int) floor(((v_pos[i] + intersect_t * v_dir[i]) +
                        step[i]*1e-7*tg.dds[i] -
                        tg.left_edge[i])/tg.dds[i]);
            tmax[i] = (((cur_ind[i]+step[i])*tg.dds[i])+
                    tg.left_edge[i]-v_pos[i])/v_dir[i];
            if((cur_ind[i] == tg.dims[i]) && (step[i] < 0)) cur_ind[i] -= 1;
            if((cur_ind[i] < 0) || (cur_ind[i] >= tg.dims[i])) skip=1;
            if(step[i] > 0) offset = 1;
            if(step[i] < 0) offset = 0;
            tmax[i] = (((cur_ind[i]+offset)*tg.dds[i])+tg.left_edge[i]-v_pos[i])/v_dir[i];
        }

        /* This is the primary grid walking loop */
        while(!(skip) &&
              !( (cur_ind[0] < 0) || (cur_ind[0] >= tg.dims[0])
              || (cur_ind[1] < 0) || (cur_ind[1] >= tg.dims[1])
              || (cur_ind[2] < 0) || (cur_ind[2] >= tg.dims[2])))
        {
            if (tmax[0] < tmax[1]) {
                if (tmax[0] < tmax[2]) {
                    direction = 0;
                } else {
                    direction = 2;
                }
            } else {
                if (tmax[1] < tmax[2]) {
                    direction = 1;
                } else {
                    direction = 2;
                }
            }
            sample_values(v_pos, v_dir, intersect_t, tmax[direction],
                    cur_ind, rgba, tf, tg);
            cur_ind[direction] += step[direction];
            intersect_t = tmax[direction];
            tmax[direction] += abs(tg.dds[direction]/v_dir[direction]);
        }

        tg.data += (tg.dims[0] * tg.dims[1] * tg.dims[2]);

    }

    __syncthreads();
    image_r[tidx] = rgba[0];
    image_g[tidx] = rgba[1];
    image_b[tidx] = rgba[2];
    image_a[tidx] = rgba[3];
}
