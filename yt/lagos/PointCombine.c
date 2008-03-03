/************************************************************************
* Copyright (C) 2007 Matthew Turk.  All Rights Reserved.
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
// PointCombine
//   A module for merging points from different grids, in various ways.
//   Used for projections, interpolations, and binning profiles.
//

#include "Python.h"

#include <stdio.h>
#include <math.h>
#include <signal.h>
#include <ctype.h>

#include "numpy/ndarrayobject.h"


// Sometimes a "maximum intensity" line-integral looks better
// switch these two defs, and then fix EnzoGrid.getProjection, to switch

//#define COMB(A,B) ((A) > (B) ? (A) : (B))
#define COMB(A,B) (A + B)

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))

static PyObject *_combineGridsError;

static PyObject *
Py_CombineGrids(PyObject *obj, PyObject *args)
{
    PyObject    *ogrid_src_x, *ogrid_src_y, *ogrid_src_vals,
        *ogrid_src_mask, *ogrid_src_wgt;
    PyObject    *ogrid_dst_x, *ogrid_dst_y, *ogrid_dst_vals,
        *ogrid_dst_mask, *ogrid_dst_wgt;

    PyArrayObject    *grid_src_x, *grid_src_y, **grid_src_vals,
            *grid_src_mask, *grid_src_wgt;
    PyArrayObject    *grid_dst_x, *grid_dst_y, **grid_dst_vals,
            *grid_dst_mask, *grid_dst_wgt;

    int NumArrays, src_len, dst_len, refinement_factor;

    if (!PyArg_ParseTuple(args, "OOOOOOOOOOi",
            &ogrid_src_x, &ogrid_src_y, 
        &ogrid_src_mask, &ogrid_src_wgt, &ogrid_src_vals,
            &ogrid_dst_x, &ogrid_dst_y,
        &ogrid_dst_mask, &ogrid_dst_wgt, &ogrid_dst_vals,
        &refinement_factor))
    return PyErr_Format(_combineGridsError,
            "CombineGrids: Invalid parameters.");

    /* First the regular source arrays */

    grid_src_x    = (PyArrayObject *) PyArray_FromAny(ogrid_src_x,
                    PyArray_DescrFromType(NPY_INT64), 1, 0,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    src_len = PyArray_SIZE(grid_src_x);

    grid_src_y    = (PyArrayObject *) PyArray_FromAny(ogrid_src_y,
                    PyArray_DescrFromType(NPY_INT64), 1, 0,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if(PyArray_SIZE(grid_src_y) != src_len) {
    PyErr_Format(_combineGridsError,
             "CombineGrids: src_x and src_y must be the same shape.");
    goto _fail;
    }

    grid_src_mask = (PyArrayObject *) PyArray_FromAny(ogrid_src_mask,
                    PyArray_DescrFromType(NPY_INT64), 1, 0,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if(PyArray_SIZE(grid_src_mask) != src_len) {
    PyErr_Format(_combineGridsError,
             "CombineGrids: src_x and src_mask must be the same shape.");
    goto _fail;
    }

    grid_src_wgt  = (PyArrayObject *) PyArray_FromAny(ogrid_src_wgt,
                    PyArray_DescrFromType(NPY_FLOAT64), 1, 0,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if(PyArray_SIZE(grid_src_wgt) != src_len) {
    PyErr_Format(_combineGridsError,
             "CombineGrids: src_x and src_wgt must be the same shape.");
    goto _fail;
    }

    grid_dst_x    = (PyArrayObject *) PyArray_FromAny(ogrid_dst_x,
                    PyArray_DescrFromType(NPY_INT64), 1, 0,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    dst_len = PyArray_SIZE(grid_dst_x);

    grid_dst_y    = (PyArrayObject *) PyArray_FromAny(ogrid_dst_y,
                    PyArray_DescrFromType(NPY_INT64), 1, 0,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if(PyArray_SIZE(grid_dst_y) != dst_len) {
    PyErr_Format(_combineGridsError,
             "CombineGrids: dst_x and dst_y must be the same shape.");
    goto _fail;
    }

    grid_dst_mask = (PyArrayObject *) PyArray_FromAny(ogrid_dst_mask,
                    PyArray_DescrFromType(NPY_INT64), 1, 0,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if(PyArray_SIZE(grid_dst_mask) != dst_len) {
    PyErr_Format(_combineGridsError,
             "CombineGrids: dst_x and dst_mask must be the same shape.");
    goto _fail;
    }

    grid_dst_wgt  = (PyArrayObject *) PyArray_FromAny(ogrid_dst_wgt,
                    PyArray_DescrFromType(NPY_FLOAT64), 1, 0,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if(PyArray_SIZE(grid_dst_wgt) != dst_len) {
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
      src_vals[i] = (npy_float64 *) grid_src_vals[i]->data;
      Py_DECREF(temp_object);

      temp_object = PySequence_GetItem(ogrid_dst_vals, i);
      grid_dst_vals[i] = (PyArrayObject *) PyArray_FromAny(
          temp_object,
          PyArray_DescrFromType(NPY_FLOAT64), 1, 0,
          NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
      dst_vals[i] = (npy_float64 *) grid_dst_vals[i]->data;
      Py_DECREF(temp_object);
    }

    /* Now we're all set to call our sub-function. */

    npy_int64     *src_x    = (npy_int64 *) grid_src_x->data;
    npy_int64     *src_y    = (npy_int64 *) grid_src_y->data;
    npy_float64 *src_wgt  = (npy_float64 *) grid_src_wgt->data;
    npy_int64     *src_mask = (npy_int64 *) grid_src_mask->data;

    npy_int64     *dst_x    = (npy_int64 *) grid_dst_x->data;
    npy_int64     *dst_y    = (npy_int64 *) grid_dst_y->data;
    npy_float64 *dst_wgt  = (npy_float64 *) grid_dst_wgt->data;
    npy_int64     *dst_mask = (npy_int64 *) grid_dst_mask->data;

    int si, di, x_off, y_off;
    npy_int64  fine_x, fine_y, init_x, init_y;
    int num_found = 0;

    for (si = 0; si < src_len; si++) {
      if (src_x[si] < 0) continue;
      init_x = refinement_factor * src_x[si];
      init_y = refinement_factor * src_y[si];
      for (x_off = 0; x_off < refinement_factor; x_off++) {
        for(y_off = 0; y_off < refinement_factor; y_off++) {
          fine_x = init_x + x_off;
          fine_y = init_y + y_off;
          for (di = 0; di < dst_len; di++) {
            if (dst_x[di] < 0) continue;
            if ((fine_x == dst_x[di]) &&
                (fine_y == dst_y[di])) {
              num_found++;
              dst_wgt[di] += src_wgt[di];
              dst_mask[di] = ((src_mask[si] && dst_mask[di]) ||
                  ((refinement_factor != 1) && (dst_mask[di])));
              // So if they are on the same level, then take the logical and
              // otherwise, set it to the destination mask
              src_x[si] = -1;
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
              (npy_float64 *) axis->data,
              table, desired,
              columns->dimensions[0],
              (npy_int32 *) columns->data,
              outputvals);
    Py_XDECREF(axis);
    Py_XDECREF(table);
    Py_XDECREF(desired);
    Py_XDECREF(outputvals);
    Py_XDECREF(columns);

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

static PyMethodDef _combineMethods[] = {
    {"CombineGrids", Py_CombineGrids, METH_VARARGS},
    {"Interpolate", Py_Interpolate, METH_VARARGS},
    {NULL, NULL} /* Sentinel */
};

/* platform independent*/
#ifdef MS_WIN32
__declspec(dllexport)
#endif

void initPointCombine(void)
{
    PyObject *m, *d;
    m = Py_InitModule("PointCombine", _combineMethods);
    d = PyModule_GetDict(m);
    _combineGridsError = PyErr_NewException("PointCombine.CombineGridsError", NULL, NULL);
    PyDict_SetItemString(d, "error", _combineGridsError);
    _interpolateError = PyErr_NewException("PointCombine.InterpolateError", NULL, NULL);
    PyDict_SetItemString(d, "error", _interpolateError);
    import_array();
}

/*
 * Local Variables:
 * mode: C
 * c-file-style: "python"
 * End:
 */
