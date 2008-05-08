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
// _MPL
//   A module for making static-resolution arrays representing
//   AMR data.
//

#include "Python.h"

#include <stdio.h>
#include <math.h>
#include <signal.h>
#include <ctype.h>

#include "numpy/ndarrayobject.h"

#define min(X,Y) ((X) < (Y) ? (X) : (Y))
#define max(X,Y) ((X) > (Y) ? (X) : (Y))

static PyObject *_pixelizeError;

char _pixelizeDocstring[] =
"Returns a static-resolution pixelized version of AMR data.\n\n"
"@parameter xp: ndarray of x centers\n"
"@Parameter yp: ndarray of y centers\n"
"@parameter dxp: ndarray of x half-widths\n"
"@parameter dyp: ndarray of y half-widths\n"
"@parameter dp: ndarray of data\n"
"@parameter rows: number of pixel rows\n"
"@parameter cols: number of pixel columns\n"
"@parameter bounds: (x_min, x_max, y_min, y_max)";

static PyObject* Py_Pixelize(PyObject *obj, PyObject *args) {

  PyObject *xp, *yp, *dxp, *dyp, *dp;
  unsigned int rows, cols;
  double x_min, x_max, y_min, y_max;

    if (!PyArg_ParseTuple(args, "OOOOOII(dddd)",
        &xp, &yp, &dxp, &dyp, &dp, &cols, &rows,
        &x_min, &x_max, &y_min, &y_max))
        return PyErr_Format(_pixelizeError, "Pixelize: Invalid Parameters.");

  double width = x_max - x_min;
  double height = y_max - y_min;
  double px_dx = width / ((double) rows);
  double px_dy = height / ((double) cols);

  // Check we have something to output to
  if (rows == 0 || cols ==0)
      PyErr_Format( _pixelizeError, "Cannot scale to zero size.");

  // Get numeric arrays
  PyArrayObject *x = (PyArrayObject *) PyArray_FromAny(xp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, NPY_C_CONTIGUOUS, NULL);
  if (x == NULL) {
      PyErr_Format( _pixelizeError, "x is of incorrect type (wanted 1D float)");
      goto _fail;
  }

  PyArrayObject *y = (PyArrayObject *) PyArray_FromAny(yp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, NPY_C_CONTIGUOUS, NULL);
  if ((y == NULL) || (PyArray_SIZE(y) != PyArray_SIZE(x))) {
      PyErr_Format( _pixelizeError, "y is of incorrect type (wanted 1D float)");
      goto _fail;
  }

  PyArrayObject *d = (PyArrayObject *) PyArray_FromAny(dp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, NPY_C_CONTIGUOUS, NULL);
  if ((d == NULL) || (PyArray_SIZE(d) != PyArray_SIZE(x))) {
      PyErr_Format( _pixelizeError, "data is of incorrect type (wanted 1D float)");
      goto _fail;
  }

  PyArrayObject *dx = (PyArrayObject *) PyArray_FromAny(dxp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, NPY_C_CONTIGUOUS, NULL);
  if ((dx == NULL) || (PyArray_SIZE(dx) != PyArray_SIZE(x))) {
      PyErr_Format( _pixelizeError, "dx is of incorrect type (wanted 1D float)");
      goto _fail;
  }
  PyArrayObject *dy = (PyArrayObject *) PyArray_FromAny(dyp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, NPY_C_CONTIGUOUS, NULL);
  if ((dy == NULL) || (PyArray_SIZE(dy) != PyArray_SIZE(x))) {
      PyErr_Format( _pixelizeError, "dy is of incorrect type (wanted 1D float)");
      goto _fail;
  }

  // Check dimensions match
  int nx = x->dimensions[0];
  int ny = y->dimensions[0];
  int ndx = dx->dimensions[0];
  int ndy = dy->dimensions[0];

  // Calculate the pointer arrays to map input x to output x
  int i, j, p;
  double lc, lr, rc, rr;
  double lypx, rypx, lxpx, rxpx, overlap1, overlap2;
  npy_float64 *xs = (npy_float64 *) PyArray_GETPTR1(x, 0);
  npy_float64 *ys = (npy_float64 *) PyArray_GETPTR1(y, 0);
  npy_float64 *dxs = (npy_float64 *) PyArray_GETPTR1(dx, 0);
  npy_float64 *dys = (npy_float64 *) PyArray_GETPTR1(dy, 0);
  npy_float64 *ds = (npy_float64 *) PyArray_GETPTR1(d, 0); // We check this above

  npy_intp dims[] = {rows, cols};
  PyArrayObject *my_array =
    (PyArrayObject *) PyArray_SimpleNewFromDescr(2, dims,
              PyArray_DescrFromType(NPY_FLOAT64));
  npy_float64 *gridded = (npy_float64 *) my_array->data;

  for(p=0;p<cols*rows;p++)gridded[p]=0.0;
  for(p=0;p<nx;p++)
  {
    if(((xs[p]+dxs[p]<x_min) ||
        (xs[p]-dxs[p]>x_max)) ||
       ((ys[p]+dys[p]<y_min) ||
        (ys[p]-dys[p]>y_max))) continue;
    lc = max(((xs[p]-dxs[p]-x_min)/px_dx),0);
    lr = max(((ys[p]-dys[p]-y_min)/px_dy),0);
    rc = min(((xs[p]+dxs[p]-x_min)/px_dx), rows);
    rr = min(((ys[p]+dys[p]-y_min)/px_dy), cols);
    for (i=lr;i<rr;i++) {
      lypx = px_dy * i + y_min;
      rypx = px_dy * (i+1) + y_min;
      overlap2 = ((min(rypx, ys[p]+dys[p]) - max(lypx, (ys[p]-dys[p])))/px_dy);
      for (j=lc;j<rc;j++) {
        lxpx = px_dx * j + x_min;
        rxpx = px_dx * (j+1) + x_min;
        overlap1 = ((min(rxpx, xs[p]+dxs[p]) - max(lxpx, (xs[p]-dxs[p])))/px_dx);
        if (overlap1 < 0.0 || overlap2 < 0.0) continue;
        gridded[j*cols+i] += (ds[p]*overlap1)*(overlap2);
      }
    }
  }

  // Attatch output buffer to output buffer

  Py_DECREF(x);
  Py_DECREF(y);
  Py_DECREF(d);
  Py_DECREF(dx);
  Py_DECREF(dy);

  PyObject *return_value = Py_BuildValue("N", my_array);

  return return_value;

  _fail:

    Py_XDECREF(x);
    Py_XDECREF(y);
    Py_XDECREF(d);
    Py_XDECREF(dx);
    Py_XDECREF(dy);
    return NULL;

}

static PyObject* Py_CPixelize(PyObject *obj, PyObject *args) {

  PyObject *xp, *yp, *zp, *pxp, *pyp,
           *dxp, *dyp, *dzp, *dp,
           *centerp, *inv_matp, *indicesp;
  unsigned int rows, cols;
  double px_min, px_max, py_min, py_max;

    if (!PyArg_ParseTuple(args, "OOOOOOOOOOOOII(dddd)",
        &xp, &yp, &zp, &pxp, &pyp, &dxp, &dyp, &dzp, &centerp, &inv_matp,
        &indicesp, &dp, &cols, &rows, &px_min, &px_max, &py_min, &py_max))
        return PyErr_Format(_pixelizeError, "CPixelize: Invalid Parameters.");

  double width = px_max - px_min;
  double height = py_max - py_min;
  long double px_dx = width / ((double) rows);
  long double px_dy = height / ((double) cols);

  // Check we have something to output to
  if (rows == 0 || cols ==0)
      PyErr_Format( _pixelizeError, "Cannot scale to zero size.");

  // Get numeric arrays
  PyArrayObject *x = (PyArrayObject *) PyArray_FromAny(xp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, NPY_C_CONTIGUOUS, NULL);
  if (x == NULL) {
      PyErr_Format( _pixelizeError, "x is of incorrect type (wanted 1D float)");
      goto _fail;
  }

  PyArrayObject *y = (PyArrayObject *) PyArray_FromAny(yp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, NPY_C_CONTIGUOUS, NULL);
  if ((y == NULL) || (PyArray_SIZE(y) != PyArray_SIZE(x))) {
      PyErr_Format( _pixelizeError, "y is of incorrect type (wanted 1D float)");
      goto _fail;
  }

  PyArrayObject *z = (PyArrayObject *) PyArray_FromAny(zp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, NPY_C_CONTIGUOUS, NULL);
  if ((z == NULL) || (PyArray_SIZE(y) != PyArray_SIZE(x))) {
      PyErr_Format( _pixelizeError, "z is of incorrect type (wanted 1D float)");
      goto _fail;
  }

  PyArrayObject *px = (PyArrayObject *) PyArray_FromAny(pxp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, NPY_C_CONTIGUOUS, NULL);
  if ((px == NULL) || (PyArray_SIZE(y) != PyArray_SIZE(x))) {
      PyErr_Format( _pixelizeError, "px is of incorrect type (wanted 1D float)");
      goto _fail;
  }

  PyArrayObject *py = (PyArrayObject *) PyArray_FromAny(pyp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, NPY_C_CONTIGUOUS, NULL);
  if ((py == NULL) || (PyArray_SIZE(y) != PyArray_SIZE(x))) {
      PyErr_Format( _pixelizeError, "py is of incorrect type (wanted 1D float)");
      goto _fail;
  }

  PyArrayObject *d = (PyArrayObject *) PyArray_FromAny(dp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, NPY_C_CONTIGUOUS, NULL);
  if ((d == NULL) || (PyArray_SIZE(d) != PyArray_SIZE(x))) {
      PyErr_Format( _pixelizeError, "data is of incorrect type (wanted 1D float)");
      goto _fail;
  }

  PyArrayObject *dx = (PyArrayObject *) PyArray_FromAny(dxp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, NPY_C_CONTIGUOUS, NULL);
  if ((dx == NULL) || (PyArray_SIZE(dx) != PyArray_SIZE(x))) {
      PyErr_Format( _pixelizeError, "dx is of incorrect type (wanted 1D float)");
      goto _fail;
  }
  PyArrayObject *dy = (PyArrayObject *) PyArray_FromAny(dyp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, NPY_C_CONTIGUOUS, NULL);
  if ((dy == NULL) || (PyArray_SIZE(dy) != PyArray_SIZE(x))) {
      PyErr_Format( _pixelizeError, "dy is of incorrect type (wanted 1D float)");
      goto _fail;
  }
  PyArrayObject *dz = (PyArrayObject *) PyArray_FromAny(dzp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, NPY_C_CONTIGUOUS, NULL);
  if ((dz == NULL) || (PyArray_SIZE(dz) != PyArray_SIZE(x))) {
      PyErr_Format( _pixelizeError, "dz is of incorrect type (wanted 1D float)");
      goto _fail;
  }
  PyArrayObject *center = (PyArrayObject *) PyArray_FromAny(centerp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, NPY_C_CONTIGUOUS, NULL);
  if ((dz == NULL) || (PyArray_SIZE(center) != 3)) {
      PyErr_Format( _pixelizeError, "Center must have three points");
      goto _fail;
  }
  PyArrayObject *inv_mat = (PyArrayObject *) PyArray_FromAny(inv_matp,
            PyArray_DescrFromType(NPY_FLOAT64), 2, 2, NPY_C_CONTIGUOUS, NULL);
  if ((inv_mat == NULL) || (PyArray_SIZE(inv_mat) != 9)) {
      PyErr_Format( _pixelizeError, "inv_mat must be three by three");
      goto _fail;
  }
  PyArrayObject *indices = (PyArrayObject *) PyArray_FromAny(indicesp,
            PyArray_DescrFromType(NPY_INT64), 1, 1, NPY_C_CONTIGUOUS, NULL);
  if ((indices == NULL) || (PyArray_SIZE(indices) != PyArray_SIZE(dx))) {
      PyErr_Format( _pixelizeError, "indices must be same length as dx");
      goto _fail;
  }

  // Check dimensions match
  int nx = x->dimensions[0];

  // Calculate the pointer arrays to map input x to output x
  int i, j, p;
  int lc, lr, rc, rr;
  long double md, cxpx, cypx;
  long double cx, cy, cz;

  npy_float64 *xs = (npy_float64 *) PyArray_GETPTR1(x, 0);
  npy_float64 *ys = (npy_float64 *) PyArray_GETPTR1(y, 0);
  npy_float64 *zs = (npy_float64 *) PyArray_GETPTR1(z, 0);
  npy_float64 *pxs = (npy_float64 *) PyArray_GETPTR1(px, 0);
  npy_float64 *pys = (npy_float64 *) PyArray_GETPTR1(py, 0);
  npy_float64 *dxs = (npy_float64 *) PyArray_GETPTR1(dx, 0);
  npy_float64 *dys = (npy_float64 *) PyArray_GETPTR1(dy, 0);
  npy_float64 *dzs = (npy_float64 *) PyArray_GETPTR1(dz, 0);
  npy_float64 *ds = (npy_float64 *) PyArray_GETPTR1(d, 0); // We check this above
  npy_float64 *centers = (npy_float64 *) PyArray_GETPTR1(center,0);
  npy_int64 *indicess = (npy_int64 *) PyArray_GETPTR1(indices,0);

  npy_intp dims[] = {rows, cols};
  PyArrayObject *my_array =
    (PyArrayObject *) PyArray_SimpleNewFromDescr(2, dims,
              PyArray_DescrFromType(NPY_FLOAT64));
  npy_float64 *gridded = (npy_float64 *) my_array->data;

  npy_float64 inv_mats[3][3];
  for(i=0;i<3;i++)for(j=0;j<3;j++)
      inv_mats[i][j]=*(npy_float64*)PyArray_GETPTR2(inv_mat,i,j);

  int pp;
  for(p=0;p<cols*rows;p++)gridded[p]=0.0;
  for(pp=0; pp<nx; pp++)
  {
    p = indicess[pp];
    // Any point we want to plot is at most this far from the center
    md = 2.0*sqrtl(dxs[p]*dxs[p] + dys[p]*dys[p] + dzs[p]*dzs[p]);
    if(((pxs[p]+md<px_min) ||
        (pxs[p]-md>px_max)) ||
       ((pys[p]+md<py_min) ||
        (pys[p]-md>py_max))) continue;
    lc = max(floorl((pxs[p]-md-px_min)/px_dx),0);
    lr = max(floorl((pys[p]-md-py_min)/px_dy),0);
    rc = min(ceill((pxs[p]+md-px_min)/px_dx),rows);
    rr = min(ceill((pys[p]+md-py_min)/px_dy),cols);
    for (i=lr;i<rr;i++) {
      cypx = px_dy * (i+0.5) + py_min;
      for (j=lc;j<rc;j++) {
        cxpx = px_dx * (j+0.5) + px_min;
        cx = inv_mats[0][0]*cxpx + inv_mats[0][1]*cypx + centers[0];
        cy = inv_mats[1][0]*cxpx + inv_mats[1][1]*cypx + centers[1];
        cz = inv_mats[2][0]*cxpx + inv_mats[2][1]*cypx + centers[2];
        if( (fabs(xs[p]-cx)>dxs[p]) || 
            (fabs(ys[p]-cy)>dys[p]) ||
            (fabs(zs[p]-cz)>dzs[p])) continue;
        gridded[j*cols+i] += ds[p];
      }
    }
  }

  // Attatch output buffer to output buffer

  Py_DECREF(x);
  Py_DECREF(y);
  Py_DECREF(z);
  Py_DECREF(px);
  Py_DECREF(py);
  Py_DECREF(d);
  Py_DECREF(dx);
  Py_DECREF(dy);
  Py_DECREF(dz);
  Py_DECREF(center);
  Py_DECREF(indices);
  Py_DECREF(inv_mat);

  PyObject *return_value = Py_BuildValue("N", my_array);

  return return_value;

  _fail:

    Py_XDECREF(x);
    Py_XDECREF(y);
    Py_XDECREF(z);
    Py_XDECREF(px);
    Py_XDECREF(py);
    Py_XDECREF(d);
    Py_XDECREF(dx);
    Py_XDECREF(dy);
    Py_XDECREF(dz);
    Py_XDECREF(center);
    Py_XDECREF(indices);
    Py_XDECREF(inv_mat);

    return NULL;

}

static PyMethodDef __MPLMethods[] = {
    {"Pixelize", Py_Pixelize, METH_VARARGS, _pixelizeDocstring},
    {"CPixelize", Py_CPixelize, METH_VARARGS, NULL},
    {NULL, NULL} /* Sentinel */
};

/* platform independent*/
#ifdef MS_WIN32
__declspec(dllexport)
#endif

void init_MPL(void)
{
    PyObject *m, *d;
    m = Py_InitModule("_MPL", __MPLMethods);
    d = PyModule_GetDict(m);
    _pixelizeError = PyErr_NewException("_MPL.error", NULL, NULL);
    PyDict_SetItemString(d, "error", _pixelizeError);
    import_array();
}
