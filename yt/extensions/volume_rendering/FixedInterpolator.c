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

#define VINDEX(A,B,C) data[((((A)+ci[0])*(ds[1]+1)+((B)+ci[1]))*(ds[2]+1)+ci[2]+(C))]
//  (((C*ds[1])+B)*ds[0]+A)

npy_float64 fast_interpolate(int *ds, int *ci, npy_float64 *dp,
                             npy_float64 *data)
{
    int i;
    npy_float64 dv, dm[3];
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
    /*assert(dv < -20);*/
    return dv;
}

npy_float64 trilinear_interpolate(int *ds, int *ci, npy_float64 *dp,
				  npy_float64 *data)
{
    /* dims is one less than the dimensions of the array */
    int i;
    npy_float64 dm[3], vz[4];
  //dp is the distance to the plane.  dm is val, dp = 1-val
    for(i=0;i<3;i++)dm[i] = (1.0 - dp[i]);
    
  //First interpolate in z
    vz[0] = dm[2]*VINDEX(0,0,0) + dp[2]*VINDEX(0,0,1);
    vz[1] = dm[2]*VINDEX(0,1,0) + dp[2]*VINDEX(0,1,1);
    vz[2] = dm[2]*VINDEX(1,0,0) + dp[2]*VINDEX(1,0,1);
    vz[3] = dm[2]*VINDEX(1,1,0) + dp[2]*VINDEX(1,1,1);

  //Then in y
    vz[0] = dm[1]*vz[0] + dp[1]*vz[1];
    vz[1] = dm[1]*vz[2] + dp[1]*vz[3];

  //Then in x
    vz[0] = dm[0]*vz[0] + dp[0]*vz[1];
    /*assert(dv < -20);*/
    return vz[0];
}
