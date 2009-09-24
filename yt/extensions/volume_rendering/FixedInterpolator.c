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

npy_float64 fast_interpolate(npy_float64 left_edge[3], npy_float64 dds[3],
               int *ds, int ci[3], npy_float64 cp[3], npy_float64 *data)
{
    npy_float64 xm, xp, ym, yp, zm, zp, dv;

    xm = (((ci[0]+1)*dds[0] + left_edge[0]) - cp[0])/dds[0];
    ym = (((ci[1]+1)*dds[1] + left_edge[1]) - cp[1])/dds[1];
    zm = (((ci[2]+1)*dds[2] + left_edge[2]) - cp[2])/dds[2];
    xp = (1.0 - xm); yp = (1.0 - ym); zp = (1.0 - zm);
    dv = data[(((ci[0]+0)*(ds[1]+1)+(ci[1]+0))*(ds[2]+1)+ci[2]+0)]*(xm*ym*zm)
       + data[(((ci[0]+1)*(ds[1]+1)+(ci[1]+0))*(ds[2]+1)+ci[2]+0)]*(xp*ym*zm)
       + data[(((ci[0]+0)*(ds[1]+1)+(ci[1]+1))*(ds[2]+1)+ci[2]+0)]*(xm*yp*zm)
       + data[(((ci[0]+0)*(ds[1]+1)+(ci[1]+0))*(ds[2]+1)+ci[2]+1)]*(xm*ym*zp)
       + data[(((ci[0]+1)*(ds[1]+1)+(ci[1]+0))*(ds[2]+1)+ci[2]+1)]*(xp*ym*zp)
       + data[(((ci[0]+0)*(ds[1]+1)+(ci[1]+1))*(ds[2]+1)+ci[2]+1)]*(xm*yp*zp)
       + data[(((ci[0]+1)*(ds[1]+1)+(ci[1]+1))*(ds[2]+1)+ci[2]+0)]*(xp*yp*zm)
       + data[(((ci[0]+1)*(ds[1]+1)+(ci[1]+1))*(ds[2]+1)+ci[2]+1)]*(xp*yp*zp);
    return dv;
}

inline void eval_shells(int nshells, npy_float64 dv,
                    npy_float64 *shells, npy_float64 rgba[4],
                    npy_float64 dt)
{
    npy_float64 dist, alpha, blend;
    int n;
    for (n = 0; n < nshells; n++) {
        dist = shells[n*6+0] - dv;
        if (dist < 0.0) dist *= -1.0;
        if (dist < shells[n*6+1]) {
            blend = shells[n*6+5]*exp(-dist/8.0)*rgba[3];
            alpha = shells[n*6+5];
            rgba[0] += blend*shells[n*6+2];
            rgba[1] += blend*shells[n*6+3];
            rgba[2] += blend*shells[n*6+4];
            rgba[3] *= (1.0 - alpha);
            break;
        }
    }
}
