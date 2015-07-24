/*******************************************************************************
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
*******************************************************************************/
//
// Some Cython versions don't like module-level constants, so we'll put them
// here.
//

#include "Python.h"

#include <stdio.h>
#include <math.h>
#include <signal.h>
#include <ctype.h>

#include "numpy/ndarrayobject.h"

#define MAX_NUM_FACES 16

#define HEX_IND     0
#define HEX_NF      6
#define TETRA_IND   1
#define TETRA_NF    4
#define WEDGE_IND   2
#define WEDGE_NF    5

extern const npy_uint8 hex_face_defs[MAX_NUM_FACES][2][2];
extern const npy_uint8 tetra_face_defs[MAX_NUM_FACES][2][2];
extern const npy_uint8 wedge_face_defs[MAX_NUM_FACES][2][2];
