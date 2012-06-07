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

#include "Python.h"

#include <stdio.h>
#include <math.h>
#include <signal.h>
#include <ctype.h>

#include "numpy/ndarrayobject.h"

npy_float64 fast_interpolate(int ds[3], int ci[3], npy_float64 dp[3],
                             npy_float64 *data);

npy_float64 offset_interpolate(int ds[3], npy_float64 dp[3], npy_float64 *data);

npy_float64 trilinear_interpolate(int ds[3], int ci[3], npy_float64 dp[3],
				  npy_float64 *data);

void eval_gradient(int ds[3], npy_float64 dp[3], npy_float64 *data, npy_float64 *grad);

void vertex_interp(npy_float64 v1, npy_float64 v2, npy_float64 isovalue,
                   npy_float64 vl[3], npy_float64 dds[3],
                   npy_float64 x, npy_float64 y, npy_float64 z,
                   int vind1, int vind2);
