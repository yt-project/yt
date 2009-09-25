/************************************************************************
* Copyright (C) 2009 Matthew Turk.  All Rights Reserved.
*
* This file is part of yt.
*
* yt is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*
************************************************************************/


//
// A small, tiny, itty bitty module for computation-intensive interpolation
// that I can't seem to make fast in Cython
//

#include "FixedInterpolator.h"

#define VINDEX(C,B,A) data[(((ci[0]+A)*(ds[1]+1)+(ci[1]+B))*(ds[2]+1)+ci[2]+C)]

npy_float64 fast_interpolate(npy_float64 left_edge[3], npy_float64 dds[3],
               int *ds, int ci[3], npy_float64 cp[3], npy_float64 *data)
{
    npy_float64 xm, x, ym, y, zm, z, dv;

    x = (((ci[0]+1)*dds[0] + left_edge[0]) - cp[0]);
    y = (((ci[1]+1)*dds[1] + left_edge[1]) - cp[1]);
    z = (((ci[2]+1)*dds[2] + left_edge[2]) - cp[2]);
    assert(0<x);assert(0<y);assert(0<z);
    assert(x<1);assert(y<1);assert(z<1);
    xm = (1.0 - x); ym = (1.0 - y); zm = (1.0 - z);
    dv  = 0.0;
    dv += VINDEX(0,0,0) * (xm*ym*zm);
    dv += VINDEX(0,0,1) * (xm*ym*z );
    dv += VINDEX(0,1,0) * (xm*y *zm);
    dv += VINDEX(0,1,1) * (xm*y *z );
    dv += VINDEX(1,0,0) * (x *ym*zm);
    dv += VINDEX(1,0,1) * (x *ym*z );
    dv += VINDEX(1,1,0) * (x *y *zm);
    dv += VINDEX(1,1,1) * (x *y *z );

    return dv;
}

inline void eval_shells(int nshells, npy_float64 dv,
                    npy_float64 *shells, npy_float64 rgba[4],
                    npy_float64 dt)
{
    npy_float64 dist, alpha, blend;
    npy_float64 ir, ig, ib, ia;
    int n;
    for (n = 0; n < nshells; n++) {
        dist = shells[n*6+0] - dv;
        if (dist < 0.0) dist *= -1.0;
        if (dist < shells[n*6+1]) {
            dist = exp(-dist/8.0);
            ir = shells[n*6+2];
            ig = shells[n*6+3];
            ib = shells[n*6+4];
            ia = shells[n*6+5];
            alpha = dt;
            rgba[0] += ia*dist*ir;
            rgba[1] += ia*dist*ig;
            rgba[2] += ia*dist*ib;
            rgba[0] *= (1.0 - alpha);
            break;
        }
    }
}
