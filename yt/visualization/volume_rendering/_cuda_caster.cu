/*******************************************************************************
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
*******************************************************************************/

// An attempt at putting the ray-casting operation into CUDA
//extern __shared__ float array[];


#define NUM_SAMPLES 25
#define VINDEX(A,B,C) tg.data[((((A)+ci[0])*(tg.dims[1]+1)+((B)+ci[1]))*(tg.dims[2]+1)+ci[2]+(C))]

#define fmin(A, B) ( A * (A < B) + B * (B < A) )
#define fmax(A, B) ( A * (A > B) + B * (B > A) )
#define fclip(A, B, C) ( fmax( fmin(A, C), B) )

struct transfer_function
{
    float vs[4][256];
    float dbin;
    float bounds[2];
};

__shared__ struct transfer_function tf;

struct grid
{
    float left_edge[3];
    float right_edge[3];
    float dds[3];
    int dims[3];
    float *data;
};

__shared__ struct grid tg;

__device__ float interpolate(float *cache, int *ds, int *ci, float *dp)
{
    int i;
    float dv, dm[3];
    for(i=0;i<3;i++)dm[i] = (1.0 - dp[i]);
    dv  = 0.0;
    dv += cache[0] * (dm[0]*dm[1]*dm[2]);
    dv += cache[1] * (dm[0]*dm[1]*dp[2]);
    dv += cache[2] * (dm[0]*dp[1]*dm[2]);
    dv += cache[3] * (dm[0]*dp[1]*dp[2]);
    dv += cache[4] * (dp[0]*dm[1]*dm[2]);
    dv += cache[5] * (dp[0]*dm[1]*dp[2]);
    dv += cache[6] * (dp[0]*dp[1]*dm[2]);
    dv += cache[7] * (dp[0]*dp[1]*dp[2]);
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
    float cache[8];

    cache[0] = VINDEX(0,0,0);
    cache[1] = VINDEX(0,0,1);
    cache[2] = VINDEX(0,1,0);
    cache[3] = VINDEX(0,1,1);
    cache[4] = VINDEX(1,0,0);
    cache[5] = VINDEX(1,0,1);
    cache[6] = VINDEX(1,1,0);
    cache[7] = VINDEX(1,1,1);

    dt = (exit_t - enter_t) / (NUM_SAMPLES-1);
    for (dti = 0; dti < NUM_SAMPLES - 1; dti++)
    {
        t = enter_t + dt*dti;
        for (i = 0; i < 3; i++)
        {
            cp[i] = v_pos[i] + t * v_dir[i];
            dp[i] = fclip(fmod(cp[i], tg.dds[i])/tg.dds[i], 0.0, 1.0);
        }
        dv = interpolate(cache, tg.dims, ci, dp);
        eval_transfer(dt, dv, rgba, tf);
    }
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
    float intersect_t = 1.0, intersect_ts[3];
    float tmax[3];
    float tl, tr, temp_xl, temp_yl, temp_xr, temp_yr;

    int offset;

    //transfer_function tf;
    for (i = 0; i < 4; i++)
    {
        x = 4 * (8 * threadIdx.x + threadIdx.y) + i;
        tf.vs[0][x] = tf_r[x];
        tf.vs[1][x] = tf_g[x];
        tf.vs[2][x] = tf_b[x];
        tf.vs[3][x] = tf_a[x];
    }

    tf.bounds[0] = tf_bounds[0]; tf.bounds[1] = tf_bounds[1];
    tf.dbin = (tf.bounds[1] - tf.bounds[0])/255.0;

    /* Set up the grid, just for convenience */
    //grid tg;
    int grid_i;

    int tidx = (blockDim.x * gridDim.x) * (
                    blockDim.y * blockIdx.y + threadIdx.y)
             + (blockDim.x * blockIdx.x + threadIdx.x);

    float rgba[4];
    //rgba[0] = image_r[tidx];
    //rgba[1] = image_g[tidx];
    //rgba[2] = image_b[tidx];
    //rgba[3] = image_a[tidx];

    float v_pos[3];
    v_pos[0] = av_pos[tidx + 0];
    v_pos[1] = av_pos[tidx + 1];
    v_pos[2] = av_pos[tidx + 2];

    tg.data = grid_data;
    int skip;
    for (i = 0; i < 3; i++)
    {
        step[i] = 0;
        step[i] +=      (v_dir[i] > 0);
        step[i] += -1 * (v_dir[i] < 0);
    }

    for(grid_i = 0; grid_i < ngrids; grid_i++) {
        skip = 0;

        if (threadIdx.x == 0)
        {
            if (threadIdx.y == 0) tg.dims[0] = dims[3*grid_i + 0];
            if (threadIdx.y == 1) tg.dims[1] = dims[3*grid_i + 1];
            if (threadIdx.y == 2) tg.dims[2] = dims[3*grid_i + 2];
        }

        if (threadIdx.x == 1)
        {
            if (threadIdx.y == 0) tg.left_edge[0] = left_edge[3*grid_i + 0];
            if (threadIdx.y == 1) tg.left_edge[1] = left_edge[3*grid_i + 1];
            if (threadIdx.y == 2) tg.left_edge[2] = left_edge[3*grid_i + 2];
        }

        if (threadIdx.x == 2) {
            if (threadIdx.y == 0) tg.right_edge[0] = right_edge[3*grid_i + 0];
            if (threadIdx.y == 1) tg.right_edge[1] = right_edge[3*grid_i + 1];
            if (threadIdx.y == 2) tg.right_edge[2] = right_edge[3*grid_i + 2];
        }

        if (threadIdx.x == 3) {
            if (threadIdx.y == 0) tg.dds[0] = (tg.right_edge[0] - tg.left_edge[0])/tg.dims[0];
            if (threadIdx.y == 1) tg.dds[1] = (tg.right_edge[1] - tg.left_edge[1])/tg.dims[1];
            if (threadIdx.y == 2) tg.dds[2] = (tg.right_edge[2] - tg.left_edge[2])/tg.dims[2];
        }

        /* We integrate our ray */

        for (i = 0; i < 3; i++)
        {
            x = (i + 1) % 3;
            y = (i + 2) % 3;

            tl = (tg.left_edge[i] - v_pos[i])/v_dir[i];
            temp_xl = (v_pos[i] + tl*v_dir[x]);
            temp_yr = (v_pos[i] + tl*v_dir[y]);

            tr = (tg.right_edge[i] - v_pos[i])/v_dir[i];
            temp_xr = (v_pos[x] + tr*v_dir[x]);
            temp_yr = (v_pos[y] + tr*v_dir[y]);

            intersect_ts[i] = 1.0;

            intersect_ts[i] += 
              ( (tg.left_edge[x] <= temp_xl) &&
                (temp_xl <= tg.right_edge[x]) &&
                (tg.left_edge[y] <= temp_yl) &&
                (temp_yl <= tg.right_edge[y]) &&
                (0.0 <= tl) && (tl < intersect_ts[i]) && (tl < tr) ) * tl;

            intersect_ts[i] += 
              ( (tg.left_edge[x] <= temp_xr) &&
                (temp_xr <= tg.right_edge[x]) &&
                (tg.left_edge[y] <= temp_yr) &&
                (temp_yr <= tg.right_edge[y]) &&
                (0.0 <= tr) && (tr < intersect_ts[i]) && (tr < tl) ) * tr;

            intersect_t = ( intersect_ts[i] < intersect_t) * intersect_ts[i];

        }

        intersect_t *= (!( (tg.left_edge[0] <= v_pos[0]) &&
                           (v_pos[0] <= tg.right_edge[0]) &&
                           (tg.left_edge[0] <= v_pos[0]) &&
                           (v_pos[0] <= tg.right_edge[0]) &&
                           (tg.left_edge[0] <= v_pos[0]) &&
                           (v_pos[0] <= tg.right_edge[0])));

        skip = ((intersect_t < 0) || (intersect_t > 1.0));

        for (i = 0; i < 3;  i++)
        {
            cur_ind[i] = (int) floor(((v_pos[i] + intersect_t * v_dir[i]) +
                        step[i]*1e-7*tg.dds[i] -
                        tg.left_edge[i])/tg.dds[i]);
            tmax[i] = (((cur_ind[i]+step[i])*tg.dds[i])+
                    tg.left_edge[i]-v_pos[i])/v_dir[i];
            cur_ind[i] -= ((cur_ind[i] == tg.dims[i]) && (step[i] < 0));
            skip = ((cur_ind[i] < 0) || (cur_ind[i] >= tg.dims[i]));
            offset = (step[i] > 0);
            tmax[i] = (((cur_ind[i]+offset)*tg.dds[i])+tg.left_edge[i]-v_pos[i])/v_dir[i];
        }

        /* This is the primary grid walking loop */
        while(!( (skip) 
              ||((cur_ind[0] < 0) || (cur_ind[0] >= tg.dims[0])
              || (cur_ind[1] < 0) || (cur_ind[1] >= tg.dims[1])
              || (cur_ind[2] < 0) || (cur_ind[2] >= tg.dims[2]))))
        {
            direction = 0;
            direction += 2 * (tmax[0] <  tmax[1]) * (tmax[0] >= tmax[2]);
            direction += 1 * (tmax[0] >= tmax[1]) * (tmax[1] <  tmax[2]);
            direction += 2 * (tmax[0] >= tmax[1]) * (tmax[1] >= tmax[2]);
            sample_values(v_pos, v_dir, intersect_t, tmax[direction],
                    cur_ind, rgba, tf, tg);
            cur_ind[direction] += step[direction];
            intersect_t = tmax[direction];
            tmax[direction] += abs(tg.dds[direction]/v_dir[direction]);
        }

        tg.data += (tg.dims[0]+1) * (tg.dims[1]+1) * (tg.dims[2]+1);

    }

    int iy = threadIdx.y + blockDim.y * blockIdx.y;
    int ix = threadIdx.x + blockDim.x * blockIdx.x;
    __syncthreads();
    
    image_r[tidx] = rgba[0];
    image_g[tidx] = rgba[1];
    image_b[tidx] = rgba[2];
    image_a[tidx] = rgba[3];
}
