/*******************************************************************************
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
*******************************************************************************/
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

npy_float64 offset_interpolate(int ds[3], npy_float64 dp[3], npy_float64 *data);

npy_float64 trilinear_interpolate(int ds[3], int ci[3], npy_float64 dp[3],
				  npy_float64 *data);

void eval_gradient(int ds[3], npy_float64 dp[3], npy_float64 *data, npy_float64 *grad);
