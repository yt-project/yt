/************************************************************************
* Copyright (C) 2007-2009 Matthew Turk.  All Rights Reserved.
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
// data_point_utilities
//   A module for merging points from different grids, in various ways.
//   Used for projections, interpolations, and binning profiles.
//

#include "Python.h"

#include <stdio.h>
#include <math.h>
#include <signal.h>
#include <ctype.h>

#include "numpy/ndarrayobject.h"

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))

static PyObject *_combineGridsError;

static PyObject *
Py_CombineGrids(PyObject *obj, PyObject *args)
{
    PyObject    *ogrid_src_x, *ogrid_src_y, *ogrid_src_vals,
        *ogrid_src_mask, *ogrid_src_wgt, *ogrid_used_mask;
    PyObject    *ogrid_dst_x, *ogrid_dst_y, *ogrid_dst_vals,
        *ogrid_dst_mask, *ogrid_dst_wgt;

    PyArrayObject    *grid_src_x, *grid_src_y, **grid_src_vals,
            *grid_src_mask, *grid_src_wgt, *grid_used_mask;
    PyArrayObject    *grid_dst_x, *grid_dst_y, **grid_dst_vals,
            *grid_dst_mask, *grid_dst_wgt;

    grid_src_x = grid_src_y = //grid_src_vals =
            grid_src_mask = grid_src_wgt = grid_used_mask =
    grid_dst_x = grid_dst_y = //grid_dst_vals = 
            grid_dst_mask = grid_dst_wgt = NULL;

    int NumArrays, src_len, dst_len, refinement_factor;
    NumArrays = 0;

    if (!PyArg_ParseTuple(args, "OOOOOOOOOOiO",
            &ogrid_src_x, &ogrid_src_y, 
        &ogrid_src_mask, &ogrid_src_wgt, &ogrid_src_vals,
            &ogrid_dst_x, &ogrid_dst_y,
        &ogrid_dst_mask, &ogrid_dst_wgt, &ogrid_dst_vals,
        &refinement_factor, &ogrid_used_mask))
    return PyErr_Format(_combineGridsError,
            "CombineGrids: Invalid parameters.");

    /* First the regular source arrays */

    grid_src_x    = (PyArrayObject *) PyArray_FromAny(ogrid_src_x,
                    PyArray_DescrFromType(NPY_INT64), 1, 1,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    src_len = PyArray_SIZE(grid_src_x);

    grid_src_y    = (PyArrayObject *) PyArray_FromAny(ogrid_src_y,
                    PyArray_DescrFromType(NPY_INT64), 1, 1,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if(PyArray_SIZE(grid_src_y) != src_len) {
    PyErr_Format(_combineGridsError,
             "CombineGrids: src_x and src_y must be the same shape.");
    goto _fail;
    }

    grid_src_mask = (PyArrayObject *) PyArray_FromAny(ogrid_src_mask,
                    PyArray_DescrFromType(NPY_INT64), 1, 1,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if(PyArray_SIZE(grid_src_mask) != src_len) {
    PyErr_Format(_combineGridsError,
             "CombineGrids: src_x and src_mask must be the same shape.");
    goto _fail;
    }

    grid_src_wgt  = (PyArrayObject *) PyArray_FromAny(ogrid_src_wgt,
                    PyArray_DescrFromType(NPY_FLOAT64), 1, 1,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if((grid_src_wgt == NULL) || (PyArray_SIZE(grid_src_wgt) != src_len)) {
    PyErr_Format(_combineGridsError,
             "CombineGrids: src_x and src_wgt must be the same shape.");
    goto _fail;
    }

    grid_used_mask  = (PyArrayObject *) PyArray_FromAny(ogrid_used_mask,
                    PyArray_DescrFromType(NPY_INT64), 1, 1,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if((grid_used_mask == NULL) || (PyArray_SIZE(grid_used_mask) != src_len)) {
    PyErr_Format(_combineGridsError,
             "CombineGrids: src_x and used_mask must be the same shape.");
    goto _fail;
    }

    grid_dst_x    = (PyArrayObject *) PyArray_FromAny(ogrid_dst_x,
                    PyArray_DescrFromType(NPY_INT64), 1, 1,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    dst_len = PyArray_SIZE(grid_dst_x);

    grid_dst_y    = (PyArrayObject *) PyArray_FromAny(ogrid_dst_y,
                    PyArray_DescrFromType(NPY_INT64), 1, 1,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if(PyArray_SIZE(grid_dst_y) != dst_len) {
    PyErr_Format(_combineGridsError,
             "CombineGrids: dst_x and dst_y must be the same shape.");
    goto _fail;
    }

    grid_dst_mask = (PyArrayObject *) PyArray_FromAny(ogrid_dst_mask,
                    PyArray_DescrFromType(NPY_INT64), 1, 1,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if(PyArray_SIZE(grid_dst_mask) != dst_len) {
    PyErr_Format(_combineGridsError,
             "CombineGrids: dst_x and dst_mask must be the same shape.");
    goto _fail;
    }

    grid_dst_wgt  = (PyArrayObject *) PyArray_FromAny(ogrid_dst_wgt,
                    PyArray_DescrFromType(NPY_FLOAT64), 1, 1,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if((grid_dst_wgt == NULL) || (PyArray_SIZE(grid_dst_wgt) != dst_len)) {
    PyErr_Format(_combineGridsError,
             "CombineGrids: dst_x and dst_wgt must be the same shape.");
    goto _fail;
    }

    /* Now we do our lists of values */
    NumArrays = PySequence_Length(ogrid_src_vals);
    if (NumArrays < 1) {
    PyErr_Format(_combineGridsError,
             "CombineGrids: You have to pass me lists of things.");
    goto _fail;
    }
    if (!(PySequence_Length(ogrid_dst_vals) == NumArrays)) {
    PyErr_Format(_combineGridsError,
             "CombineGrids: Sorry, but your lists of values are different lengths.");
    goto _fail;
    }

    grid_src_vals = malloc(NumArrays * sizeof(PyArrayObject*));
    grid_dst_vals = malloc(NumArrays * sizeof(PyArrayObject*));
    npy_float64 **src_vals = malloc(NumArrays * sizeof(npy_float64*));
    npy_float64 **dst_vals = malloc(NumArrays * sizeof(npy_float64*));
    PyObject *temp_object;
    int i;
    for (i = 0; i < NumArrays; i++) {
      temp_object = PySequence_GetItem(ogrid_src_vals, i);
      grid_src_vals[i] = (PyArrayObject *) PyArray_FromAny(
          temp_object,
          PyArray_DescrFromType(NPY_FLOAT64), 1, 0,
          NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
      src_vals[i] = (npy_float64 *) PyArray_GETPTR1(grid_src_vals[i],0);
      Py_DECREF(temp_object);

      temp_object = PySequence_GetItem(ogrid_dst_vals, i);
      grid_dst_vals[i] = (PyArrayObject *) PyArray_FromAny(
          temp_object,
          PyArray_DescrFromType(NPY_FLOAT64), 1, 0,
          NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
      dst_vals[i] = (npy_float64 *) PyArray_GETPTR1(grid_dst_vals[i],0);
      Py_DECREF(temp_object);
    }

    /* Now we're all set to call our sub-function. */

    npy_int64     *src_x    = (npy_int64 *) PyArray_GETPTR1(grid_src_x,0);
    npy_int64     *src_y    = (npy_int64 *) PyArray_GETPTR1(grid_src_y,0);
    npy_float64 *src_wgt  = (npy_float64 *) PyArray_GETPTR1(grid_src_wgt,0);
    npy_int64     *src_mask = (npy_int64 *) PyArray_GETPTR1(grid_src_mask,0);
    npy_int64    *src_used_mask = (npy_int64 *) PyArray_GETPTR1(grid_used_mask,0);

    npy_int64     *dst_x    = (npy_int64 *) PyArray_GETPTR1(grid_dst_x,0);
    npy_int64     *dst_y    = (npy_int64 *) PyArray_GETPTR1(grid_dst_y,0);
    npy_float64 *dst_wgt  = (npy_float64 *) PyArray_GETPTR1(grid_dst_wgt,0);
    npy_int64     *dst_mask = (npy_int64 *) PyArray_GETPTR1(grid_dst_mask,0);

    int si, di, x_off, y_off;
    npy_int64  fine_x, fine_y, init_x, init_y;
    int num_found = 0;

    for (si = 0; si < src_len; si++) {
      if (src_used_mask[si] == 0) continue;
      init_x = refinement_factor * src_x[si];
      init_y = refinement_factor * src_y[si];
      for (x_off = 0; x_off < refinement_factor; x_off++) {
        for(y_off = 0; y_off < refinement_factor; y_off++) {
          fine_x = init_x + x_off;
          fine_y = init_y + y_off;
          for (di = 0; di < dst_len; di++) {
            if ((fine_x == dst_x[di]) &&
                (fine_y == dst_y[di])) {
              num_found++;
              dst_wgt[di] += src_wgt[si];
              dst_mask[di] = ((src_mask[si] && dst_mask[di]) ||
                  ((refinement_factor != 1) && (dst_mask[di])));
              // So if they are on the same level, then take the logical and
              // otherwise, set it to the destination mask
              src_used_mask[si] = 0;
              for (i = 0; i < NumArrays; i++) {
                dst_vals[i][di] += src_vals[i][si];
              }
              if (refinement_factor == 1) break;
            }
          }
        }
      }
    }

    Py_DECREF(grid_src_x);
    Py_DECREF(grid_src_y);
    Py_DECREF(grid_src_mask);
    Py_DECREF(grid_src_wgt);
    Py_DECREF(grid_used_mask);

    Py_DECREF(grid_dst_x);
    Py_DECREF(grid_dst_y);
    Py_DECREF(grid_dst_mask);
    Py_DECREF(grid_dst_wgt);

    if (NumArrays > 0){
      for (i = 0; i < NumArrays; i++) {
        Py_DECREF(grid_src_vals[i]);
        Py_DECREF(grid_dst_vals[i]);
      }
    }

    free(grid_src_vals);
    free(grid_dst_vals);
    free(src_vals);
    free(dst_vals);

    PyObject *onum_found = PyInt_FromLong((long)num_found);
    return onum_found;

_fail:
    Py_XDECREF(grid_src_x);
    Py_XDECREF(grid_src_y);
    Py_XDECREF(grid_src_wgt);
    Py_XDECREF(grid_src_mask);
    Py_XDECREF(grid_used_mask);

    Py_XDECREF(grid_dst_x);
    Py_XDECREF(grid_dst_y);
    Py_XDECREF(grid_dst_wgt);
    Py_XDECREF(grid_dst_mask);
    if (NumArrays > 0){
      for (i = 0; i < NumArrays; i++) {
        Py_XDECREF(grid_src_vals[i]);
        Py_XDECREF(grid_dst_vals[i]);
      }
    }
    return NULL;

}

static PyObject *_dataCubeError;

static PyObject *DataCubeGeneric(PyObject *obj, PyObject *args,
                   void (*to_call)(PyArrayObject* c_data, npy_int64 xc,
                                        npy_int64 yc, npy_int64 zc,
                                   PyArrayObject* g_data, npy_int64 xg,
                                        npy_int64 yg, npy_int64 zg))
{
    /* Standard boilerplate unpacking...  */

    /* 
       rf              (py_int)                 i
       grid_leftedge   (npy_float64 COERCE)     O
       dx_grid         (npy_float64 COERCE)     O
       griddata        (npy_float64 array)      O
       childmask       (npy_bool array)         O
       cube_leftedge   (npy_float64 COERCE)     O
       cube_rightedge  (npy_float64 COERCE)     O
       dx_cube         (npy_float64 COERCE)     O
       cubedata        (npy_float64 array)      O
       lastlevel       (py_int)                 i
    */


    int ll, i, n;

    PyObject *og_le, *og_dx, *og_data, *og_cm,
             *oc_le, *oc_re, *oc_dx, *oc_data, *odr_edge, *odl_edge;
    PyArrayObject *g_le, *g_dx, *g_cm,
                  *c_le, *c_re, *c_dx, *dr_edge, *dl_edge;
    g_dx=g_cm=c_le=c_re=c_dx=NULL;
    PyArrayObject **g_data, **c_data;
    g_data = c_data = NULL;
    npy_int *ag_cm;
    npy_float64 ag_le[3], ag_dx[3], 
                ac_le[3], ac_re[3], ac_dx[3],
                adl_edge[3], adr_edge[3];
    Py_ssize_t n_fields = 0;

    if (!PyArg_ParseTuple(args, "OOOOOOOOiOO",
            &og_le, &og_dx, &og_data, &og_cm,
            &oc_le, &oc_re, &oc_dx, &oc_data,
        &ll, &odl_edge, &odr_edge))
    return PyErr_Format(_dataCubeError,
            "DataCubeGeneric: Invalid parameters.");

    /* First the regular source arrays */

    g_le    = (PyArrayObject *) PyArray_FromAny(og_le,
                    PyArray_DescrFromType(NPY_FLOAT64), 1, 1,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if((g_le==NULL) || (PyArray_SIZE(g_le) != 3)) {
    PyErr_Format(_dataCubeError,
             "CombineGrids: Three values, one dimension required for g_le.");
    goto _fail;
    }
    for(i=0;i<3;i++)ag_le[i]=*(npy_float64*)PyArray_GETPTR1(g_le,i);

    g_dx    = (PyArrayObject *) PyArray_FromAny(og_dx,
                    PyArray_DescrFromType(NPY_FLOAT64), 1, 1,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if((g_dx==NULL) || (PyArray_SIZE(g_dx) != 3)) {
    PyErr_Format(_dataCubeError,
             "CombineGrids: Three values, one dimension required for g_dx.");
    goto _fail;
    }
    for(i=0;i<3;i++)ag_dx[i]=*(npy_float64*)PyArray_GETPTR1(g_dx,i);

    g_cm    = (PyArrayObject *) PyArray_FromAny(og_cm,
                    PyArray_DescrFromType(NPY_INT), 3, 3,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if((g_cm==NULL) || (g_cm->nd != 3)) {
    PyErr_Format(_dataCubeError,
             "CombineGrids: Three dimensions required for g_cm.");
    goto _fail;
    }
    ag_cm = (npy_int*) g_cm->data;

    /* Now the cube */
 
    c_le    = (PyArrayObject *) PyArray_FromAny(oc_le,
                    PyArray_DescrFromType(NPY_FLOAT64), 1, 1,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if((c_le==NULL) || (PyArray_SIZE(c_le) != 3)) {
    PyErr_Format(_dataCubeError,
             "CombineGrids: Three values, one dimension required for c_le.");
    goto _fail;
    }
    for(i=0;i<3;i++)ac_le[i]=*(npy_float64*)PyArray_GETPTR1(c_le,i);

    c_re    = (PyArrayObject *) PyArray_FromAny(oc_re,
                    PyArray_DescrFromType(NPY_FLOAT64), 1, 1,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if((c_re==NULL) || (PyArray_SIZE(c_re) != 3)) {
    PyErr_Format(_dataCubeError,
             "CombineGrids: Three values, one dimension required for c_re.");
    goto _fail;
    }
    for(i=0;i<3;i++)ac_re[i]=*(npy_float64*)PyArray_GETPTR1(c_re,i);

    c_dx    = (PyArrayObject *) PyArray_FromAny(oc_dx,
                    PyArray_DescrFromType(NPY_FLOAT64), 1, 1,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if((c_dx==NULL) || (PyArray_SIZE(c_dx) != 3)) {
    PyErr_Format(_dataCubeError,
             "CombineGrids: Three values, one dimension required for c_dx.");
    goto _fail;
    }
    for(i=0;i<3;i++)ac_dx[i]=*(npy_float64*)PyArray_GETPTR1(c_dx,i);

    if (!PyList_Check(oc_data)){
      PyErr_Format(_dataCubeError,
          "CombineGrids: c_data must be a list of arrays!");
      goto _fail;
    }

    dl_edge    = (PyArrayObject *) PyArray_FromAny(odl_edge,
                    PyArray_DescrFromType(NPY_FLOAT64), 1, 1,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if((dl_edge==NULL) || (PyArray_SIZE(dl_edge) != 3)) {
    PyErr_Format(_dataCubeError,
             "CombineGrids: Three values, one dimension required for dl_edge.");
    goto _fail;
    }
    for(i=0;i<3;i++)adl_edge[i]=*(npy_float64*)PyArray_GETPTR1(dl_edge,i);

    dr_edge    = (PyArrayObject *) PyArray_FromAny(odr_edge,
                    PyArray_DescrFromType(NPY_FLOAT64), 1, 1,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if((dr_edge==NULL) || (PyArray_SIZE(dr_edge) != 3)) {
    PyErr_Format(_dataCubeError,
             "CombineGrids: Three values, one dimension required for dr_edge.");
    goto _fail;
    }
    for(i=0;i<3;i++)adr_edge[i]=*(npy_float64*)PyArray_GETPTR1(dr_edge,i);

    n_fields = PyList_Size(oc_data);
    if(n_fields == 0) {
      PyErr_Format(_dataCubeError,
          "CombineGrids: Length zero for c_data is invalid.");
      goto _fail;
    }

    PyObject *tc_data;
    c_data = (PyArrayObject**)
             malloc(sizeof(PyArrayObject*)*n_fields);
    for (n=0;n<n_fields;n++)c_data[n]=NULL;
    for (n=0;n<n_fields;n++){
      tc_data = (PyObject *) PyList_GetItem(oc_data, n);
      c_data[n]    = (PyArrayObject *) PyArray_FromAny(tc_data,
          PyArray_DescrFromType(NPY_FLOAT64), 3, 3,
          NPY_UPDATEIFCOPY, NULL);
      if((c_data[n]==NULL) || (c_data[n]->nd != 3)) {
        PyErr_Format(_dataCubeError,
            "CombineGrids: Three dimensions required for c_data[%i].",n);
        goto _fail;
      }
    }

    if ((!PyList_Check(og_data)) ||
         (PyList_Size(og_data) != n_fields)){
      PyErr_Format(_dataCubeError,
          "CombineGrids: g_data must be a list of arrays same length as c_data!");
      goto _fail;
    }

    PyObject *tg_data;
    g_data = (PyArrayObject**)
             malloc(sizeof(PyArrayObject*)*n_fields);
    for (n=0;n<n_fields;n++)g_data[n]=NULL;
    for (n=0;n<n_fields;n++){
      /* Borrows a reference */
      tg_data = (PyObject *) PyList_GetItem(og_data, n);
      /* We set up an array so we only have to do this once 
         Note that this is an array in the C sense, not the NumPy sense */
      g_data[n]    = (PyArrayObject *) PyArray_FromAny(tg_data,
          PyArray_DescrFromType(NPY_FLOAT64), 3, 3,
          NPY_UPDATEIFCOPY, NULL);
      if((g_data[n]==NULL) || (g_data[n]->nd != 3)) {
        PyErr_Format(_dataCubeError,
            "CombineGrids: Three dimensions required for g_data[%i].",n);
        goto _fail;
      }
    }

    /* And let's begin */

    npy_int64 xg, yg, zg, xc, yc, zc, cmax_x, cmax_y, cmax_z,
              cmin_x, cmin_y, cmin_z, cm, pxl, pyl, pzl;
    long int total=0;

    int p_niter[3] = {1,1,1};
    int itc;
    npy_float64 ac_le_p[3][3];
    npy_float64 ac_re_p[3][3];
    npy_float64 ag_re[3];
    /* This is for checking for periodic boundary conditions.
       Manually set the right edge to be offset from the left. */
    for(i=0;i<3;i++){ag_re[i] = ag_le[i]+ag_dx[i]*(g_data[0]->dimensions[i]+1);}

    for(i=0;i<3;i++){ac_le_p[i][0] = ac_le[i]; ac_re_p[i][0] = ac_re[i];}
    for(i=0;i<3;i++) {
            itc = 1;
            if (ac_le[i] < adl_edge[i]) {
                ac_le_p[i][itc  ] = adr_edge[i] - (adl_edge[i] - ac_le[i]);
                ac_re_p[i][itc++] = adr_edge[i] + (ac_re[i] - adl_edge[i]);
            }
            if (ac_re[i] > adr_edge[i]) {
                ac_le_p[i][itc  ] = ac_le[i] - (adr_edge[i] - adl_edge[i]);
                ac_re_p[i][itc++] = adl_edge[i] + (ac_re[i] - adr_edge[i]);
            }
            p_niter[i] = itc;
    }
    npy_intp nx, ny, nz;
    /* This is easier than doing a lookup every loop */
    nx = PyArray_DIM(c_data[0], 0);
    ny = PyArray_DIM(c_data[0], 1);
    nz = PyArray_DIM(c_data[0], 2);
    npy_int64 xg_min, yg_min, zg_min;
    npy_int64 xg_max, yg_max, zg_max;

    /* Periodic iterations, *if necessary* */
    for (pxl = 0; pxl < p_niter[0]; pxl++) {
    xg_min = max(floor((ac_le_p[0][pxl]-ag_le[0])/ag_dx[0]) - 1, 0);
    xg_max = min(ceil((ac_re_p[0][pxl]-ag_le[0])/ag_dx[0]) + 1, 
                 g_data[0]->dimensions[0]);
    for (xg = xg_min; xg < xg_max; xg++) {
      /* If we're off the destination cell boundary, skip */
      if (ag_le[0]+ag_dx[0]*xg     > ac_re_p[0][pxl]) continue;
      if (ag_le[0]+ag_dx[0]*(xg+1) < ac_le_p[0][pxl]) continue;
      /* Floor to the source edge */
      cmin_x = max(floorl((ag_le[0]+ag_dx[0]*xg     - ac_le_p[0][pxl])/ac_dx[0]),0);
      cmax_x = min( ceill((ag_le[0]+ag_dx[0]*(xg+1) - ac_le_p[0][pxl])/ac_dx[0]),nx);
      if(cmin_x==cmax_x)continue;
      for (pyl = 0; pyl < p_niter[1]; pyl++) {
      yg_min = max(floor((ac_le_p[1][pyl]-ag_le[1])/ag_dx[1]) - 1, 0);
      yg_max = min(ceil((ac_re_p[1][pyl]-ag_le[1])/ag_dx[1]),
                   g_data[0]->dimensions[1]);
      for (yg = yg_min; yg < yg_max; yg++) {
        if (ag_le[1]+ag_dx[1]*yg     > ac_re_p[1][pyl]) continue;
        if (ag_le[1]+ag_dx[1]*(yg+1) < ac_le_p[1][pyl]) continue;
        cmin_y = max(floorl((ag_le[1]+ag_dx[1]*yg     - ac_le_p[1][pyl])/ac_dx[1]),0);
        cmax_y = min( ceill((ag_le[1]+ag_dx[1]*(yg+1) - ac_le_p[1][pyl])/ac_dx[1]),ny);
        if(cmin_y==cmax_y)continue;
        for (pzl = 0; pzl < p_niter[2]; pzl++) {
        zg_min = max(floor((ac_le_p[2][pzl]-ag_le[2])/ag_dx[2]) - 1, 0);
        zg_max = min(ceil((ac_re_p[2][pzl]-ag_le[2])/ag_dx[2]), 
                     g_data[0]->dimensions[2]);
        for (zg = zg_min; zg < zg_max; zg++) {
        cm = *(npy_int *)PyArray_GETPTR3(g_cm,xg,yg,zg);
        if ((!ll) && (cm == 0)) continue;
          if (ag_le[2]+ag_dx[2]*zg     > ac_re_p[2][pzl]) continue;
          if (ag_le[2]+ag_dx[2]*(zg+1) < ac_le_p[2][pzl]) continue;
          cmin_z = max(floorl((ag_le[2]+ag_dx[2]*zg     - ac_le_p[2][pzl])/ac_dx[2]),0);
          cmax_z = min( ceill((ag_le[2]+ag_dx[2]*(zg+1) - ac_le_p[2][pzl])/ac_dx[2]),nz);
          if(cmin_z==cmax_z)continue;
          for (xc = cmin_x; xc < cmax_x ; xc++) {
            for (yc = cmin_y; yc < cmax_y ; yc++) {
              for (zc = cmin_z; zc < cmax_z ; zc++) {
                for(n=0;n<n_fields;n++){
                  to_call(c_data[n], xc, yc, zc,
                          g_data[n], xg, yg, zg);
                }
                total += 1;
              }
            }
          }
        }}
      }}
    }}

    /* Cleanup time */

    Py_DECREF(g_le);
    Py_DECREF(g_dx);
    Py_DECREF(g_cm);
    Py_DECREF(c_le);
    Py_DECREF(c_re);
    Py_DECREF(c_dx);
    Py_DECREF(dl_edge);
    Py_DECREF(dr_edge);
    for(n=0;n<n_fields;n++) {
        Py_DECREF(g_data[n]);
        Py_DECREF(c_data[n]);
    }
    free(g_data);
    free(c_data);

    PyObject *status = PyInt_FromLong(total);
    return status;
    
_fail:
    Py_XDECREF(g_le);
    Py_XDECREF(g_dx);
    Py_XDECREF(g_cm);
    Py_XDECREF(c_le);
    Py_XDECREF(c_re);
    Py_XDECREF(c_dx);
    for(n=0;n<n_fields;n++) {
        if(g_data[n]!=NULL){Py_XDECREF(g_data[n]);}
        if(c_data[n]!=NULL){Py_XDECREF(c_data[n]);}
    }
    if(g_data!=NULL)free(g_data);
    if(c_data!=NULL)free(c_data);
    return NULL;

}

/* These functions are both called with
    func(cubedata, griddata) */

static void dcNothing(PyArrayObject* c_data, npy_int64 xc, npy_int64 yc, npy_int64 zc,
                     PyArrayObject* g_data, npy_int64 xg, npy_int64 yg, npy_int64 zg)
{
    return;
}

/* These functions are both called with
    func(cubedata, griddata) */

static void dcRefine(PyArrayObject* c_data, npy_int64 xc, npy_int64 yc, npy_int64 zc,
                     PyArrayObject* g_data, npy_int64 xg, npy_int64 yg, npy_int64 zg)
{
    // c_data used to be val1, g_data used to be val2
    // so we go val2 -> val1, or grid -> covering grid
    *(npy_float64*) PyArray_GETPTR3(c_data,xc,yc,zc) =
        *(npy_float64*) PyArray_GETPTR3(g_data,xg,yg,zg);
}

static void dcReplace(PyArrayObject* c_data, npy_int64 xc, npy_int64 yc, npy_int64 zc,
                      PyArrayObject* g_data, npy_int64 xg, npy_int64 yg, npy_int64 zg)
{
    // c_data used to be val1, g_data used to be val2
    // so we go val1 -> val2, or covering grid -> grid
    *(npy_float64*) PyArray_GETPTR3(g_data,xg,yg,zg) =
        *(npy_float64*) PyArray_GETPTR3(c_data,xc,yc,zc);
}

static PyObject *
Py_DataCubeRefine(PyObject *obj, PyObject *args)
{
    PyObject* to_return = DataCubeGeneric(obj, args, dcRefine);
    return to_return;
}

static PyObject *
Py_DataCubeReplace(PyObject *obj, PyObject *args)
{
    PyObject* to_return = DataCubeGeneric(obj, args, dcReplace);
    return to_return;
}

static PyObject *Py_FillRegion(PyObject *obj, PyObject *args)
{
    PyObject *oc_data, *og_data,
             *oc_start, *og_start,
             *oc_dims, *og_dims, *omask;
    PyObject *tg_data, *tc_data, *dw_data;
    oc_data = og_data = oc_start = og_start = oc_dims = og_dims = omask = NULL;
    tg_data = tc_data = dw_data = NULL;
    PyArrayObject **g_data, **c_data, *mask,
                  *g_start, *c_start, *c_dims, *g_dims, *dwa;
    mask = g_start = c_start = c_dims = g_dims = NULL;
    g_data = c_data = NULL;
    int refratio, ll, direction, n;
    npy_int64 gxs, gys, gzs, gxe, gye, gze;
    npy_int64 cxs, cys, czs, cxe, cye, cze;
    npy_int64 ixs, iys, izs, ixe, iye, ize;
    npy_int64 gxi, gyi, gzi, cxi, cyi, czi;
    npy_int64 cdx, cdy, cdz;
    npy_int64 dw[3];
    int i;
    npy_int64 ci, cj, ck, ri, rj, rk;
    int total = 0;
    void (*to_call)(PyArrayObject* c_data, npy_int64 xc,
                         npy_int64 yc, npy_int64 zc,
                    PyArrayObject* g_data, npy_int64 xg,
                         npy_int64 yg, npy_int64 zg);
    if (!PyArg_ParseTuple(args, "iOOOOOOOOii",
            &refratio, &og_start, &oc_start,
            &oc_data, &og_data,
            &oc_dims, &og_dims, &omask, &dw_data, &ll, &direction))
    return PyErr_Format(_dataCubeError,
            "DataCubeGeneric: Invalid parameters.");

    if (direction == 0) to_call = dcRefine;
    else if (direction == 1) to_call = dcReplace;

    g_start = (PyArrayObject *) PyArray_FromAny(og_start,
                PyArray_DescrFromType(NPY_INT64), 1, 1, 0, NULL);
    if(g_start == NULL){
      PyErr_Format(_dataCubeError, "FillRegion: g_start invalid.");
      goto _fail;
    }

    c_start = (PyArrayObject *) PyArray_FromAny(oc_start,
                PyArray_DescrFromType(NPY_INT64), 1, 1, 0, NULL);
    if(c_start == NULL){
      PyErr_Format(_dataCubeError, "FillRegion: c_start invalid.");
      goto _fail;
    }

    g_dims  = (PyArrayObject *) PyArray_FromAny(og_dims,
                PyArray_DescrFromType(NPY_INT32), 1, 1, 0, NULL);
    if(g_dims == NULL){
      PyErr_Format(_dataCubeError, "FillRegion: g_dims invalid.");
      goto _fail;
    }

    c_dims  = (PyArrayObject *) PyArray_FromAny(oc_dims,
                PyArray_DescrFromType(NPY_INT32), 1, 1, 0, NULL);
    if(c_dims == NULL){
      PyErr_Format(_dataCubeError, "FillRegion: c_dims invalid.");
      goto _fail;
    }

    mask    = (PyArrayObject *) PyArray_FromAny(omask,
                PyArray_DescrFromType(NPY_INT32), 3, 3, 0, NULL);
    if(mask == NULL){
      PyErr_Format(_dataCubeError, "FillRegion: mask invalid.");
      goto _fail;
    }

    dwa     = (PyArrayObject *) PyArray_FromAny(dw_data,
                PyArray_DescrFromType(NPY_INT64), 1, 1, 0, NULL);
    if(dwa == NULL){
      PyErr_Format(_dataCubeError, "FillRegion: domain width invalid.");
      goto _fail;
    }
    for (i=0;i<3;i++)dw[i] = *(npy_int64*) PyArray_GETPTR1(dwa, i);

    int n_fields = PyList_Size(oc_data);
    if(n_fields == 0) {
      /*PyErr_Format(_dataCubeError,
          "CombineGrids: Length zero for c_data is invalid.");
      goto _fail;*/
    }
    if (!PyList_Check(og_data) || (PyList_Size(og_data) != n_fields)){
      PyErr_Format(_dataCubeError,
          "CombineGrids: g_data must be a list of arrays same length as c_data!");
      goto _fail;
    }

    c_data = (PyArrayObject**)
             malloc(sizeof(PyArrayObject*)*n_fields);
    g_data = (PyArrayObject**)
             malloc(sizeof(PyArrayObject*)*n_fields);
    for (n=0;n<n_fields;n++)c_data[n]=g_data[n]=NULL;
    for (n=0;n<n_fields;n++){
      tg_data = (PyObject *) PyList_GetItem(og_data, n);
      /* We set up an array so we only have to do this once 
         Note that this is an array in the C sense, not the NumPy sense */
      g_data[n]    = (PyArrayObject *) PyArray_FromAny(tg_data,
          PyArray_DescrFromType(NPY_FLOAT64), 3, 3,
          NPY_UPDATEIFCOPY, NULL);
      if((g_data[n]==NULL) || (g_data[n]->nd != 3)) {
        PyErr_Format(_dataCubeError,
            "CombineGrids: Three dimensions required for g_data[%i].",n);
        goto _fail;
      }
      tc_data = (PyObject *) PyList_GetItem(oc_data, n);
      c_data[n]    = (PyArrayObject *) PyArray_FromAny(tc_data,
          PyArray_DescrFromType(NPY_FLOAT64), 3, 3,
          NPY_UPDATEIFCOPY, NULL);
      if((c_data[n]==NULL) || (c_data[n]->nd != 3)) {
        PyErr_Format(_dataCubeError,
            "CombineGrids: Three dimensions required for c_data[%i].",n);
        goto _fail;
      }
    }

    /* g[xyz][se] are the start and end index in integers 
       of the grid, at its refinement level               */
    gxs = *(npy_int64 *) PyArray_GETPTR1(og_start, 0);
    gys = *(npy_int64 *) PyArray_GETPTR1(og_start, 1);
    gzs = *(npy_int64 *) PyArray_GETPTR1(og_start, 2);
    gxe = gxs + *(npy_int32 *) PyArray_GETPTR1(og_dims, 0);
    gye = gys + *(npy_int32 *) PyArray_GETPTR1(og_dims, 1);
    gze = gzs + *(npy_int32 *) PyArray_GETPTR1(og_dims, 2);

    /* c[xyz][se] are the start and end index in integers
       of the covering grid, at its refinement level      */
    cxs = *(npy_int64 *) PyArray_GETPTR1(oc_start, 0);
    cys = *(npy_int64 *) PyArray_GETPTR1(oc_start, 1);
    czs = *(npy_int64 *) PyArray_GETPTR1(oc_start, 2);

    /* cd[xyz] are the dimensions of the covering grid */
    cdx = (*(npy_int32 *) PyArray_GETPTR1(oc_dims, 0));
    cdy = (*(npy_int32 *) PyArray_GETPTR1(oc_dims, 1));
    cdz = (*(npy_int32 *) PyArray_GETPTR1(oc_dims, 2));
    cxe = (cxs + cdx);
    cye = (cys + cdy);
    cze = (czs + cdz);

    /* It turns out that C89 doesn't define a mechanism for choosing the sign
       of the remainder.
    */
        //fprintf(stderr, "ci == %d, cxi == %d, dw[0] == %d\n", (int) ci, (int) cxi, (int) dw[0]);
    for(cxi=cxs;cxi<cxe;cxi++) {
        ci = (cxi % dw[0]);
        ci = (ci < 0) ? ci + dw[0] : ci;
        if ( ci < gxs*refratio || ci >= gxe*refratio) continue;
        gxi = floor(ci / refratio) - gxs;
        for(cyi=cys;cyi<cye;cyi++) {
            cj = cyi % dw[1];
            cj = (cj < 0) ? cj + dw[1] : cj;
            if ( cj < gys*refratio || cj >= gye*refratio) continue;
            gyi = floor(cj / refratio) - gys;
            for(czi=czs;czi<cze;czi++) {
                ck = czi % dw[2];
                ck = (ck < 0) ? ck + dw[2] : ck;
                if ( ck < gzs*refratio || ck >= gze*refratio) continue;
                gzi = floor(ck / refratio) - gzs;
                    if ((ll) || (*(npy_int32*)PyArray_GETPTR3(mask, gxi,gyi,gzi) > 0)) 
                {
                if (direction!=2)
                  for(n=0;n<n_fields;n++){
                      to_call(c_data[n],
                          cxi - cxs, cyi - cys, czi - czs,
                          g_data[n], gxi, gyi, gzi);
                  }
                total += 1;
                }
            }
        }
    }

    Py_DECREF(g_start);
    Py_DECREF(c_start);
    Py_DECREF(g_dims);
    Py_DECREF(c_dims);
    Py_DECREF(mask);
    for(n=0;n<n_fields;n++) {
        Py_DECREF(g_data[n]);
        Py_DECREF(c_data[n]);
    }
    free(g_data);
    free(c_data);
    PyObject *status = PyInt_FromLong(total);
    return status;

_fail:
    Py_XDECREF(g_start);
    Py_XDECREF(c_start);
    Py_XDECREF(g_dims);
    Py_XDECREF(c_dims);
    Py_XDECREF(mask);
    if(g_data != NULL)
      for(n=0;n<n_fields;n++)if(g_data[n]!=NULL){Py_XDECREF(g_data[n]);}
    if(c_data != NULL)
      for(n=0;n<n_fields;n++)if(c_data[n]!=NULL){Py_XDECREF(c_data[n]);}
    if(g_data!=NULL)free(g_data);
    if(c_data!=NULL)free(c_data);
    return NULL;
}

static PyObject *Py_FillBuffer(PyObject *obj, PyObject *args)
{
    PyObject *oc_data, *og_data,
             *oc_start, *og_start,
             *oc_dims, *og_dims, *omask, *odls;
    PyObject *tg_data, *tc_data, *dw_data;
    oc_data = og_data = oc_start = og_start = oc_dims = og_dims = omask = NULL;
    tg_data = tc_data = dw_data = odls = NULL;
    PyArrayObject **g_data, **c_data, *mask,
                  *g_start, *c_start, *c_dims, *g_dims, *dwa;
    mask = g_start = c_start = c_dims = g_dims = NULL;
    double *dls = NULL;
    int refratio, ll, direction, n;
    npy_int64 gxs, gys, gzs, gxe, gye, gze;
    npy_int64 cxs, cys, czs, cxe, cye, cze;
    npy_int64 ixs, iys, izs, ixe, iye, ize;
    npy_int64 gxi, gyi, gzi, cxi, cyi, czi;
    npy_int64 cdx, cdy, cdz;
    npy_int64 dw[3];
    int i, axis;
    int ci, cj, ck, ri, rj, rk;
    int total = 0;

    if (!PyArg_ParseTuple(args, "iOOOOOOOOOi",
            &refratio, &og_start, &oc_start,
            &oc_data, &og_data,
            &oc_dims, &og_dims, &omask, &dw_data, &odls, &axis))
    return PyErr_Format(_dataCubeError,
            "DataCubeGeneric: Invalid parameters.");

    g_start = (PyArrayObject *) PyArray_FromAny(og_start,
                PyArray_DescrFromType(NPY_INT64), 1, 1, 0, NULL);
    if(g_start == NULL){
      PyErr_Format(_dataCubeError, "FillRegion: g_start invalid.");
      goto _fail;
    }

    c_start = (PyArrayObject *) PyArray_FromAny(oc_start,
                PyArray_DescrFromType(NPY_INT64), 1, 1, 0, NULL);
    if(c_start == NULL){
      PyErr_Format(_dataCubeError, "FillRegion: c_start invalid.");
      goto _fail;
    }

    g_dims  = (PyArrayObject *) PyArray_FromAny(og_dims,
                PyArray_DescrFromType(NPY_INT32), 1, 1, 0, NULL);
    if(g_dims == NULL){
      PyErr_Format(_dataCubeError, "FillRegion: g_dims invalid.");
      goto _fail;
    }

    c_dims  = (PyArrayObject *) PyArray_FromAny(oc_dims,
                PyArray_DescrFromType(NPY_INT32), 1, 1, 0, NULL);
    if(c_dims == NULL){
      PyErr_Format(_dataCubeError, "FillRegion: c_dims invalid.");
      goto _fail;
    }

    mask    = (PyArrayObject *) PyArray_FromAny(omask,
                PyArray_DescrFromType(NPY_INT32), 3, 3, 0, NULL);
    if(mask == NULL){
      PyErr_Format(_dataCubeError, "FillRegion: mask invalid.");
      goto _fail;
    }

    dwa     = (PyArrayObject *) PyArray_FromAny(dw_data,
                PyArray_DescrFromType(NPY_INT64), 1, 1, 0, NULL);
    if(dwa == NULL){
      PyErr_Format(_dataCubeError, "FillRegion: domain width invalid.");
      goto _fail;
    }
    for (i=0;i<3;i++)dw[i] = *(npy_int64*) PyArray_GETPTR1(dwa, i);

    int n_fields = PyList_Size(oc_data);
    if(n_fields == 0) {
      PyErr_Format(_dataCubeError,
          "CombineGrids: Length zero for c_data is invalid.");
      goto _fail;
    }
    if (!PyList_Check(og_data) || (PyList_Size(og_data) != n_fields)){
      PyErr_Format(_dataCubeError,
          "CombineGrids: g_data must be a list of arrays same length as c_data!");
      goto _fail;
    }
    if (!PyList_Check(odls) || (PyList_Size(odls) != n_fields)){
      PyErr_Format(_dataCubeError,
          "CombineGrids: dls must be a list of arrays same length as c_data!");
      goto _fail;
    }

    c_data = (PyArrayObject**)
             malloc(sizeof(PyArrayObject*)*n_fields);
    g_data = (PyArrayObject**)
             malloc(sizeof(PyArrayObject*)*n_fields);
    dls = (double *) malloc(sizeof(double) * n_fields);
    PyObject *temp = NULL;
    for (n=0;n<n_fields;n++)c_data[n]=g_data[n]=NULL;
    for (n=0;n<n_fields;n++){
      /* Borrowed reference ... */
      temp = PyList_GetItem(odls, n);
      dls[n] = PyFloat_AsDouble(temp);

      tg_data = (PyObject *) PyList_GetItem(og_data, n);
      /* We set up an array so we only have to do this once 
         Note that this is an array in the C sense, not the NumPy sense */
      g_data[n]    = (PyArrayObject *) PyArray_FromAny(tg_data,
          PyArray_DescrFromType(NPY_FLOAT64), 3, 3,
          NPY_UPDATEIFCOPY, NULL);
      if((g_data[n]==NULL) || (g_data[n]->nd != 3)) {
        PyErr_Format(_dataCubeError,
            "CombineGrids: Three dimensions required for g_data[%i].",n);
        goto _fail;
      }
      tc_data = (PyObject *) PyList_GetItem(oc_data, n);
      c_data[n]    = (PyArrayObject *) PyArray_FromAny(tc_data,
          PyArray_DescrFromType(NPY_FLOAT64), 2, 2,
          NPY_UPDATEIFCOPY, NULL);
      if((c_data[n]==NULL) || (c_data[n]->nd != 2)) {
        PyErr_Format(_dataCubeError,
            "CombineGrids: Two dimensions required for c_data[%i].",n);
        goto _fail;
      }
    }

    /* g[xyz][se] are the start and end index in integers 
       of the grid, at its refinement level               */
    gxs = *(npy_int64 *) PyArray_GETPTR1(og_start, 0);
    gys = *(npy_int64 *) PyArray_GETPTR1(og_start, 1);
    gzs = *(npy_int64 *) PyArray_GETPTR1(og_start, 2);
    gxe = gxs + *(npy_int32 *) PyArray_GETPTR1(og_dims, 0);
    gye = gys + *(npy_int32 *) PyArray_GETPTR1(og_dims, 1);
    gze = gzs + *(npy_int32 *) PyArray_GETPTR1(og_dims, 2);

    /* c[xyz][se] are the start and end index in integers
       of the covering grid, at its refinement level      */
    cxs = *(npy_int64 *) PyArray_GETPTR1(oc_start, 0);
    cys = *(npy_int64 *) PyArray_GETPTR1(oc_start, 1);
    czs = *(npy_int64 *) PyArray_GETPTR1(oc_start, 2);

    /* cd[xyz] are the dimensions of the covering grid */
    cdx = (*(npy_int32 *) PyArray_GETPTR1(oc_dims, 0));
    cdy = (*(npy_int32 *) PyArray_GETPTR1(oc_dims, 1));
    cdz = (*(npy_int32 *) PyArray_GETPTR1(oc_dims, 2));
    cxe = (cxs + cdx - 1);
    cye = (cys + cdy - 1);
    cze = (czs + cdz - 1);

    /* It turns out that C89 doesn't define a mechanism for choosing the sign
       of the remainder.
    */
    int x_loc, y_loc; // For access into the buffer
    for(cxi=cxs;cxi<=cxe;cxi++) {
        ci = (cxi % dw[0]);
        ci = (ci < 0) ? ci + dw[0] : ci;
        if ( ci < gxs*refratio || ci >= gxe*refratio) continue;
        gxi = floor(ci / refratio) - gxs;
        for(cyi=cys;cyi<=cye;cyi++) {
            cj = cyi % dw[1];
            cj = (cj < 0) ? cj + dw[1] : cj;
            if ( cj < gys*refratio || cj >= gye*refratio) continue;
            gyi = floor(cj / refratio) - gys;
            for(czi=czs;czi<=cze;czi++) {
                ck = czi % dw[2];
                ck = (ck < 0) ? ck + dw[2] : ck;
                if ( ck < gzs*refratio || ck >= gze*refratio) continue;
                gzi = floor(ck / refratio) - gzs;
                    if (refratio == 1 || *(npy_int32*)PyArray_GETPTR3(mask, gxi,gyi,gzi) > 0)
                {
                switch (axis) {
                  case 0: x_loc = cyi-cys; y_loc = czi-czs; break;
                  case 1: x_loc = cxi-cxs; y_loc = czi-czs; break;
                  case 2: x_loc = cxi-cys; y_loc = cyi-cys; break;
                }
                for(n=0;n<n_fields;n++){
                    *(npy_float64*) PyArray_GETPTR2(c_data[n], x_loc, y_loc)
                    +=  *(npy_float64*) PyArray_GETPTR3(g_data[n], gxi, gyi, gzi) 
                        * dls[n] / refratio;
                }
                total += 1;
                }
            }
        }
    }

    Py_DECREF(g_start);
    Py_DECREF(c_start);
    Py_DECREF(g_dims);
    Py_DECREF(c_dims);
    Py_DECREF(mask);
    for(n=0;n<n_fields;n++) {
        Py_DECREF(g_data[n]);
        Py_DECREF(c_data[n]);
    }
    if(dls!=NULL)free(dls);
    if(g_data!=NULL)free(g_data);
    if(c_data!=NULL)free(c_data);
    PyObject *status = PyInt_FromLong(total);
    return status;

_fail:
    Py_XDECREF(g_start);
    Py_XDECREF(c_start);
    Py_XDECREF(g_dims);
    Py_XDECREF(c_dims);
    Py_XDECREF(mask);
    if(dls!=NULL)free(dls);
    for(n=0;n<n_fields;n++) {
        if(g_data!=NULL)if(g_data[n]!=NULL){Py_XDECREF(g_data[n]);}
        if(c_data!=NULL)if(c_data[n]!=NULL){Py_XDECREF(c_data[n]);}
    }
    if(g_data!=NULL)free(g_data);
    if(c_data!=NULL)free(c_data);
    return NULL;
}

static PyObject *_findContoursError;

int process_neighbors(PyArrayObject*, npy_int64, npy_int64, npy_int64,
                            int first);
static PyObject *
Py_FindContours(PyObject *obj, PyObject *args)
{
    PyObject *ocon_ids, *oxi, *oyi, *ozi;
    PyArrayObject *con_ids, *xi, *yi, *zi;
    xi=yi=zi=con_ids=NULL;
    npy_int64 i, j, k, n;
    int status;

    i = 0;
    if (!PyArg_ParseTuple(args, "OOOO",
        &ocon_ids, &oxi, &oyi, &ozi))
        return PyErr_Format(_findContoursError,
                    "FindContours: Invalid parameters.");
    
    con_ids   = (PyArrayObject *) PyArray_FromAny(ocon_ids,
                    PyArray_DescrFromType(NPY_INT64), 3, 3,
                    0 | NPY_UPDATEIFCOPY, NULL);
    if((con_ids==NULL) || (con_ids->nd != 3)) {
    PyErr_Format(_findContoursError,
             "FindContours: Three dimensions required for con_ids.");
    goto _fail;
    }

    xi = (PyArrayObject *) PyArray_FromAny(oxi,
                    PyArray_DescrFromType(NPY_INT64), 1, 1,
                    0, NULL);
    if(xi==NULL) {
    PyErr_Format(_findContoursError,
             "FindContours: One dimension required for xi.");
    goto _fail;
    }
    
    yi = (PyArrayObject *) PyArray_FromAny(oyi,
                    PyArray_DescrFromType(NPY_INT64), 1, 1,
                    0, NULL);
    if((yi==NULL) || (PyArray_SIZE(xi) != PyArray_SIZE(yi))) {
    PyErr_Format(_findContoursError,
             "FindContours: One dimension required for yi, same size as xi.");
    goto _fail;
    }
    
    zi = (PyArrayObject *) PyArray_FromAny(ozi,
                    PyArray_DescrFromType(NPY_INT64), 1, 1,
                    0, NULL);
    if((zi==NULL) || (PyArray_SIZE(xi) != PyArray_SIZE(zi))) {
    PyErr_Format(_findContoursError,
             "FindContours: One dimension required for zi, same size as xi.");
    goto _fail;
    }
    
    for(n=0;n<xi->dimensions[0];n++) {
      i=*(npy_int64 *)PyArray_GETPTR1(xi,n);
      j=*(npy_int64 *)PyArray_GETPTR1(yi,n);
      k=*(npy_int64 *)PyArray_GETPTR1(zi,n);
      status = process_neighbors(con_ids, i, j, k, 1);
      if(status < 0) break;
    }

    Py_DECREF(con_ids);
    Py_DECREF(xi);
    Py_DECREF(yi);
    Py_DECREF(zi);

    PyObject *retval = PyInt_FromLong(status);
    return retval;

    _fail:
        Py_XDECREF(con_ids);
        Py_XDECREF(xi);
        Py_XDECREF(yi);
        Py_XDECREF(zi);
        return NULL;
}

int process_neighbors(PyArrayObject *con_ids, npy_int64 i, npy_int64 j,
                            npy_int64 k, int first)
{
  npy_int64 off_i, off_j, off_k;
  int spawn_check, status;
  int mi, mj, mk;
  static int stack_depth;
  if (first == 1) stack_depth = 0;
  else stack_depth++;
  if (stack_depth > 10000) return -1;
  npy_int64 *fd_off, *fd_ijk;
  mi = con_ids->dimensions[0];
  mj = con_ids->dimensions[1];
  mk = con_ids->dimensions[2];
  fd_ijk = ((npy_int64*)PyArray_GETPTR3(con_ids, i, j, k));
  //if(*fd_ijk == -1){return ((npy_int64)*fd_ijk);}
  do {
    spawn_check = 0;
    for (off_i=max(i-1,0);off_i<=min(i+1,mi-1);off_i++)
      for (off_j=max(j-1,0);off_j<=min(j+1,mj-1);off_j++)
        for (off_k=max(k-1,0);off_k<=min(k+1,mk-1);off_k++) {
          if((off_i==i)&&(off_j==j)&&(off_k==k)) continue;
          fd_off = ((npy_int64*)PyArray_GETPTR3(con_ids, off_i, off_j, off_k));
          if(*fd_off == -1) continue;
          if(*fd_off > *fd_ijk){
            *fd_ijk = *fd_off;
            spawn_check += 1;
          }
          if(*fd_off < *fd_ijk){
            *fd_off = *fd_ijk;
            status = process_neighbors(con_ids, off_i, off_j, off_k, 0);
            if (*fd_off != *fd_ijk) spawn_check += 1;
            *fd_ijk = *fd_off;
            if (status < 0) return -1;
          }
        }
  } while (spawn_check > 0);
  stack_depth -= 1;
  return 1;
}

static PyObject *_interpolateError;

static void
Interpolate(long num_axis_points, npy_float64 *axis, PyArrayObject* table,
            PyArrayObject *desiredvals, long num_columns, npy_int32 *columns,
            PyArrayObject *outputvals)
{
    //int table_rows = table->dimensions[0];
    int num_desireds = desiredvals->dimensions[0];

    npy_int axis_ind, col_ind;
    npy_int32 column;
    npy_int64 desired_num;

    npy_float64 desired;

    npy_float64 logtem0 = log10(axis[0]);
    npy_float64 logtem9 = log10(axis[num_axis_points-1]);
    npy_float64 dlogtem = (logtem9-logtem0)/(num_axis_points-1);
    npy_float64 t1, t2, tdef, ki, kip;
    npy_float64 *t;

    for (desired_num = 0 ; desired_num < num_desireds ; desired_num++) {
        t = (npy_float64*)PyArray_GETPTR1(desiredvals, desired_num);
        desired = log10l(*t);
        axis_ind = min(num_axis_points-1,
                   max(0,(int)((desired-logtem0)/dlogtem)+1));
        t1 = (logtem0 + (axis_ind-1)*dlogtem);
        t2 = (logtem0 + (axis_ind+0)*dlogtem);
        tdef = t2 - t1;
        for (column = 0 ; column < num_columns ; column++) {
            col_ind = (npy_int) columns[column];
            ki  = *(npy_float64*)PyArray_GETPTR2(table, (npy_int) (axis_ind-1), col_ind);
            kip = *(npy_float64*)PyArray_GETPTR2(table, (npy_int) (axis_ind+0), col_ind);
            *(npy_float64*) PyArray_GETPTR2(outputvals, desired_num, column) =
                    ki+(desired-t1)*(kip-ki)/tdef;
        }
    }
    return;
}

static PyObject *
Py_Interpolate(PyObject *obj, PyObject *args)
{
    PyObject   *oaxis, *otable, *odesired, *ooutputvals, *ocolumns;
    PyArrayObject   *axis, *table, *desired, *outputvals, *columns;

    if (!PyArg_ParseTuple(args, "OOOOO",
        &oaxis, &otable, &odesired, &ooutputvals, &ocolumns))
        return PyErr_Format(_interpolateError,
                    "Interpolate: Invalid parameters.");

    /* Align, Byteswap, Contiguous, Typeconvert */
    axis          =  (PyArrayObject *) PyArray_FromAny(oaxis         , PyArray_DescrFromType(NPY_FLOAT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL );
    table         =  (PyArrayObject *) PyArray_FromAny(otable        , PyArray_DescrFromType(NPY_FLOAT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL );
    desired       =  (PyArrayObject *) PyArray_FromAny(odesired      , PyArray_DescrFromType(NPY_FLOAT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL );
    outputvals    =  (PyArrayObject *) PyArray_FromAny(ooutputvals   , PyArray_DescrFromType(NPY_FLOAT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL );
    columns       =  (PyArrayObject *) PyArray_FromAny(ocolumns      ,   PyArray_DescrFromType(NPY_INT32), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL );

    if (!axis || !table || !desired || !outputvals || !columns) {
        PyErr_Format(_interpolateError,
                  "Interpolate: error converting array inputs.");
        goto _fail;
    }

    if (columns->dimensions[0] != outputvals->dimensions[1]) {
        PyErr_Format(_interpolateError,
                 "Interpolate: number of columns requested must match number "
                 "of columns in output buffer. %i", (int) columns->dimensions[0]);
        goto _fail;
    }

    Interpolate(axis->dimensions[0],
              (npy_float64 *) PyArray_GETPTR1(axis, 0),
              table, desired,
              columns->dimensions[0],
              (npy_int32 *) PyArray_GETPTR1(columns, 0),
              outputvals);
    Py_DECREF(axis);
    Py_DECREF(table);
    Py_DECREF(desired);
    Py_DECREF(outputvals);
    Py_DECREF(columns);

    /* Align, Byteswap, Contiguous, Typeconvert */
    return Py_None;

  _fail:
    Py_XDECREF(axis);
    Py_XDECREF(table);
    Py_XDECREF(desired);
    Py_XDECREF(outputvals);
    Py_XDECREF(columns);


    return NULL;
}

static PyObject *_findBindingEnergyError;

static PyObject *
Py_FindBindingEnergy(PyObject *obj, PyObject *args)
{
    PyObject *omass, *ox, *oy, *oz;
    PyArrayObject *mass, *x, *y, *z;
    x=y=z=mass=NULL;
    int truncate;
    double kinetic_energy;

    if (!PyArg_ParseTuple(args, "OOOOid",
        &omass, &ox, &oy, &oz, &truncate, &kinetic_energy))
        return PyErr_Format(_findBindingEnergyError,
                    "FindBindingEnergy: Invalid parameters.");
    
    mass   = (PyArrayObject *) PyArray_FromAny(omass,
                    PyArray_DescrFromType(NPY_FLOAT64), 1, 1,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if((mass==NULL) || (mass->nd != 1)) {
    PyErr_Format(_findBindingEnergyError,
             "FindBindingEnergy: One dimension required for mass.");
    goto _fail;
    }

    x      = (PyArrayObject *) PyArray_FromAny(ox,
                    PyArray_DescrFromType(NPY_FLOAT64), 1, 1,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if((x==NULL) || (x->nd != 1) 
        || (PyArray_SIZE(x) != PyArray_SIZE(mass))) {
    PyErr_Format(_findBindingEnergyError,
             "FindBindingEnergy: x must be same size as mass.");
    goto _fail;
    }

    y      = (PyArrayObject *) PyArray_FromAny(oy,
                    PyArray_DescrFromType(NPY_FLOAT64), 1, 1,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if((y==NULL) || (y->nd != 1) 
        || (PyArray_SIZE(y) != PyArray_SIZE(mass))) {
    PyErr_Format(_findBindingEnergyError,
             "FindBindingEnergy: y must be same size as mass.");
    goto _fail;
    }

    z      = (PyArrayObject *) PyArray_FromAny(oz,
                    PyArray_DescrFromType(NPY_FLOAT64), 1, 1,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if((z==NULL) || (z->nd != 1) 
        || (PyArray_SIZE(z) != PyArray_SIZE(mass))) {
    PyErr_Format(_findBindingEnergyError,
             "FindBindingEnergy: z must be same size as mass.");
    goto _fail;
    }

    /* Do the work here. */
    int q_outer, q_inner, n_q = PyArray_SIZE(mass);
    double this_potential, total_potential;
    total_potential = 0;
    npy_float64 mass_o, x_o, y_o, z_o;
    npy_float64 mass_i, x_i, y_i, z_i;

    /* progress bar stuff */
    float totalWork = 0.5 * (pow(n_q,2.0) - n_q);
    float workDone = 0;
    int every_cells = floor(n_q / 100);
    int until_output = 1;
    for (q_outer = 0; q_outer < n_q - 1; q_outer++) {
        this_potential = 0;
        mass_o = *(npy_float64*) PyArray_GETPTR1(mass, q_outer);
        x_o = *(npy_float64*) PyArray_GETPTR1(x, q_outer);
        y_o = *(npy_float64*) PyArray_GETPTR1(y, q_outer);
        z_o = *(npy_float64*) PyArray_GETPTR1(z, q_outer);
        for (q_inner = q_outer+1; q_inner < n_q; q_inner++) {
            mass_i = *(npy_float64*) PyArray_GETPTR1(mass, q_inner);
            x_i = *(npy_float64*) PyArray_GETPTR1(x, q_inner);
            y_i = *(npy_float64*) PyArray_GETPTR1(y, q_inner);
            z_i = *(npy_float64*) PyArray_GETPTR1(z, q_inner);
            this_potential += mass_o * mass_i / 
                            sqrtl( (x_i-x_o)*(x_i-x_o)
                                 + (y_i-y_o)*(y_i-y_o)
                                 + (z_i-z_o)*(z_i-z_o) );
        }
        total_potential += this_potential;
	workDone += n_q - q_outer - 1;
    until_output -= 1;
    if(until_output == 0){
        fprintf(stderr,"Calculating Potential for %i cells: %.2f%%\t(pe/ke = %e)\r",
                n_q,((100*workDone)/totalWork),(total_potential/kinetic_energy));
        fflush(stdout); fflush(stderr);
        until_output = every_cells;
    }
        if ((truncate == 1) && (total_potential > kinetic_energy)){
            fprintf(stderr, "Truncating!\r");
            break;
        }
    }
    fprintf(stderr,"\n");
    fflush(stdout); fflush(stderr);

    Py_DECREF(mass);
    Py_DECREF(x);
    Py_DECREF(y);
    Py_DECREF(z);
    PyObject *status = PyFloat_FromDouble(total_potential);
    return status;

    _fail:
        Py_XDECREF(mass);
        Py_XDECREF(x);
        Py_XDECREF(y);
        Py_XDECREF(z);
        return NULL;
}

static PyObject *_outputFloatsToFileError;

static PyObject *
Py_OutputFloatsToFile(PyObject *obj, PyObject *args)
{
    PyObject *oarray;
    PyArrayObject *array;
    char *filename, *header = NULL;
    npy_intp i, j, imax, jmax;

    if (!PyArg_ParseTuple(args, "Os|s", &oarray, &filename, &header))
        return PyErr_Format(_outputFloatsToFileError,
                    "OutputFloatsToFile: Invalid parameters.");

    array   = (PyArrayObject *) PyArray_FromAny(oarray,
                    PyArray_DescrFromType(NPY_FLOAT64), 2, 2,
                    0, NULL);
    if(array==NULL){
    PyErr_Format(_outputFloatsToFileError,
             "OutputFloatsToFile: Failure to convert array ( nd == 2 ?)");
    goto _fail;
    }

    FILE *to_write = fopen(filename, "w");
    if(to_write == NULL){
    PyErr_Format(_outputFloatsToFileError,
             "OutputFloatsToFile: Unable to open %s for writing.", filename);
      goto _fail;
    }

    if(header!=NULL)fprintf(to_write,"%s\n", header);

    imax = PyArray_DIM(array, 0);
    jmax = PyArray_DIM(array, 1);
    for(i=0;i<imax;i++){
      for(j=0;j<jmax;j++){
        fprintf(to_write, "%0.16e",
                *((npy_float64*)PyArray_GETPTR2(array,i,j)));
        if(j<jmax-1)fprintf(to_write, "\t");
      }
      fprintf(to_write, "\n");
    }
    fclose(to_write);

    Py_DECREF(array);
    return Py_None;

   _fail:
    Py_XDECREF(array);

    return NULL;
}


static PyMethodDef _combineMethods[] = {
    {"CombineGrids", Py_CombineGrids, METH_VARARGS},
    {"Interpolate", Py_Interpolate, METH_VARARGS},
    {"DataCubeRefine", Py_DataCubeRefine, METH_VARARGS},
    {"DataCubeReplace", Py_DataCubeReplace, METH_VARARGS},
    {"FindContours", Py_FindContours, METH_VARARGS},
    {"FindBindingEnergy", Py_FindBindingEnergy, METH_VARARGS},
    {"OutputFloatsToFile", Py_OutputFloatsToFile, METH_VARARGS},
    {"FillRegion", Py_FillRegion, METH_VARARGS},
    {"FillBuffer", Py_FillBuffer, METH_VARARGS},
    {NULL, NULL} /* Sentinel */
};

/* platform independent*/
#ifdef MS_WIN32
__declspec(dllexport)
#endif

void initdata_point_utilities(void)
{
    PyObject *m, *d;
    m = Py_InitModule("data_point_utilities", _combineMethods);
    d = PyModule_GetDict(m);
    _combineGridsError = PyErr_NewException("data_point_utilities.CombineGridsError", NULL, NULL);
    PyDict_SetItemString(d, "error", _combineGridsError);
    _interpolateError = PyErr_NewException("data_point_utilities.InterpolateError", NULL, NULL);
    PyDict_SetItemString(d, "error", _interpolateError);
    _dataCubeError = PyErr_NewException("data_point_utilities.DataCubeError", NULL, NULL);
    PyDict_SetItemString(d, "error", _dataCubeError);
    _findContoursError = PyErr_NewException("data_point_utilities.FindContoursError", NULL, NULL);
    PyDict_SetItemString(d, "error", _findContoursError);
    _outputFloatsToFileError = PyErr_NewException("data_point_utilities.OutputFloatsToFileError", NULL, NULL);
    PyDict_SetItemString(d, "error", _outputFloatsToFileError);
    import_array();
}

/*
 * Local Variables:
 * mode: C
 * c-file-style: "python"
 * End:
 */
