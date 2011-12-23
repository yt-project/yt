/************************************************************************
* Copyright (C) 2007-2011 Matthew Turk.  All Rights Reserved.
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
  PyArrayObject *x, *y, *dx, *dy, *d;
  xp = yp = dxp = dyp = dp = NULL;
  x = y = dx = dy = d = NULL;
  unsigned int rows, cols;
  int antialias = 1;
  double x_min, x_max, y_min, y_max;
  double period_x, period_y;
  period_x = period_y = 0;
  int check_period = 1;

  if (!PyArg_ParseTuple(args, "OOOOOII(dddd)|i(dd)i",
      &xp, &yp, &dxp, &dyp, &dp, &cols, &rows,
      &x_min, &x_max, &y_min, &y_max,
      &antialias, &period_x, &period_y, &check_period))
      return PyErr_Format(_pixelizeError, "Pixelize: Invalid Parameters.");

  double width = x_max - x_min;
  double height = y_max - y_min;
  double px_dx = width / ((double) rows);
  double px_dy = height / ((double) cols);
  double ipx_dx = 1.0 / px_dx;
  double ipx_dy = 1.0 / px_dy;

  // Check we have something to output to
  if (rows == 0 || cols ==0)
      PyErr_Format( _pixelizeError, "Cannot scale to zero size.");

  // Get numeric arrays
  x = (PyArrayObject *) PyArray_FromAny(xp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, 0, NULL);
  if (x == NULL) {
      PyErr_Format( _pixelizeError, "x is of incorrect type (wanted 1D float)");
      goto _fail;
  }

  y = (PyArrayObject *) PyArray_FromAny(yp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, 0, NULL);
  if ((y == NULL) || (PyArray_SIZE(y) != PyArray_SIZE(x))) {
      PyErr_Format( _pixelizeError, "y is of incorrect type (wanted 1D float)");
      goto _fail;
  }

  d = (PyArrayObject *) PyArray_FromAny(dp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, 0, NULL);
  if ((d == NULL) || (PyArray_SIZE(d) != PyArray_SIZE(x))) {
      PyErr_Format( _pixelizeError, "data is of incorrect type (wanted 1D float)");
      goto _fail;
  }

  dx = (PyArrayObject *) PyArray_FromAny(dxp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, 0, NULL);
  if ((dx == NULL) || (PyArray_SIZE(dx) != PyArray_SIZE(x))) {
      PyErr_Format( _pixelizeError, "dx is of incorrect type (wanted 1D float)");
      goto _fail;
  }
  dy = (PyArrayObject *) PyArray_FromAny(dyp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, 0, NULL);
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
  int i, j, p, xi, yi;
  double lc, lr, rc, rr;
  double lypx, rypx, lxpx, rxpx, overlap1, overlap2;
  npy_float64 oxsp, oysp, xsp, ysp, dxsp, dysp, dsp;
  int xiter[2], yiter[2];
  double xiterv[2], yiterv[2];

  

  npy_intp dims[] = {rows, cols};
  PyArrayObject *my_array =
    (PyArrayObject *) PyArray_SimpleNewFromDescr(2, dims,
              PyArray_DescrFromType(NPY_FLOAT64));
  //npy_float64 *gridded = (npy_float64 *) my_array->data;

  xiter[0] = yiter[0] = 0;
  xiterv[0] = yiterv[0] = 0.0;

  Py_BEGIN_ALLOW_THREADS
  for(i=0;i<rows;i++)for(j=0;j<cols;j++)
      *(npy_float64*) PyArray_GETPTR2(my_array, i, j) = 0.0;
  for(p=0;p<nx;p++)
  {
    // these are cell-centered
    oxsp = *((npy_float64 *)PyArray_GETPTR1(x, p));
    oysp = *((npy_float64 *)PyArray_GETPTR1(y, p));
    // half-width
    dxsp = *((npy_float64 *)PyArray_GETPTR1(dx, p));
    dysp = *((npy_float64 *)PyArray_GETPTR1(dy, p));
    dsp = *((npy_float64 *)PyArray_GETPTR1(d, p));
    xiter[1] = yiter[1] = 999;
    if(check_period == 1) {
      if (oxsp < x_min) {xiter[1] = +1; xiterv[1] = period_x;}
      else if (oxsp > x_max) {xiter[1] = -1; xiterv[1] = -period_x;}
      if (oysp < y_min) {yiter[1] = +1; yiterv[1] = period_y;}
      else if (oysp > y_max) {yiter[1] = -1; yiterv[1] = -period_y;}
    }
    for(xi = 0; xi < 2; xi++) {
      if(xiter[xi] == 999)continue;
      xsp = oxsp + xiterv[xi];
      if((xsp+dxsp<x_min) || (xsp-dxsp>x_max)) continue;
      for(yi = 0; yi < 2; yi++) {
        if(yiter[yi] == 999)continue;
        ysp = oysp + yiterv[yi];
        if((ysp+dysp<y_min) || (ysp-dysp>y_max)) continue;
        lc = max(((xsp-dxsp-x_min)*ipx_dx),0);
        lr = max(((ysp-dysp-y_min)*ipx_dy),0);
        rc = min(((xsp+dxsp-x_min)*ipx_dx), rows);
        rr = min(((ysp+dysp-y_min)*ipx_dy), cols);
        for (i=lr;i<rr;i++) {
          lypx = px_dy * i + y_min;
          rypx = px_dy * (i+1) + y_min;
          overlap2 = ((min(rypx, ysp+dysp) - max(lypx, (ysp-dysp)))*ipx_dy);
          for (j=lc;j<rc;j++) {
            lxpx = px_dx * j + x_min;
            rxpx = px_dx * (j+1) + x_min;
            overlap1 = ((min(rxpx, xsp+dxsp) - max(lxpx, (xsp-dxsp)))*ipx_dx);
            if (overlap1 < 0.0 || overlap2 < 0.0) continue;
            if (antialias == 1)
              *(npy_float64*) PyArray_GETPTR2(my_array, j, i) +=
                    (dsp*overlap1)*overlap2;
            else *(npy_float64*) PyArray_GETPTR2(my_array, j, i) = dsp;
          }
        }
      }
    }
  }
  Py_END_ALLOW_THREADS

  // Attatch output buffer to output buffer

  Py_DECREF(x);
  Py_DECREF(y);
  Py_DECREF(d);
  Py_DECREF(dx);
  Py_DECREF(dy);

  PyObject *return_value = Py_BuildValue("N", my_array);

  return return_value;

  _fail:

    if(x!=NULL)Py_XDECREF(x);
    if(y!=NULL)Py_XDECREF(y);
    if(d!=NULL)Py_XDECREF(d);
    if(dx!=NULL)Py_XDECREF(dx);
    if(dy!=NULL)Py_XDECREF(dy);
    return NULL;

}

static PyObject* Py_CPixelize(PyObject *obj, PyObject *args) {

  PyObject *xp, *yp, *zp, *pxp, *pyp,
           *dxp, *dyp, *dzp, *dp,
           *centerp, *inv_matp, *indicesp;

  xp = yp = zp = pxp = pyp = dxp = dyp = dzp = dp = NULL;
  centerp = inv_matp = indicesp = NULL;

  PyArrayObject *x, *y, *z, *px, *py, *d,
                *dx, *dy, *dz, *center, *inv_mat, *indices;

  x = y = z = px = py = dx = dy = dz = d = NULL;
  center = inv_mat = indices = NULL;

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
  x = (PyArrayObject *) PyArray_FromAny(xp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, 0, NULL);
  if (x == NULL) {
      PyErr_Format( _pixelizeError, "x is of incorrect type (wanted 1D float)");
      goto _fail;
  }

  y = (PyArrayObject *) PyArray_FromAny(yp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, 0, NULL);
  if ((y == NULL) || (PyArray_SIZE(y) != PyArray_SIZE(x))) {
      PyErr_Format( _pixelizeError, "y is of incorrect type (wanted 1D float)");
      goto _fail;
  }

  z = (PyArrayObject *) PyArray_FromAny(zp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, 0, NULL);
  if ((z == NULL) || (PyArray_SIZE(y) != PyArray_SIZE(x))) {
      PyErr_Format( _pixelizeError, "z is of incorrect type (wanted 1D float)");
      goto _fail;
  }

  px = (PyArrayObject *) PyArray_FromAny(pxp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, 0, NULL);
  if ((px == NULL) || (PyArray_SIZE(y) != PyArray_SIZE(x))) {
      PyErr_Format( _pixelizeError, "px is of incorrect type (wanted 1D float)");
      goto _fail;
  }

  py = (PyArrayObject *) PyArray_FromAny(pyp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, 0, NULL);
  if ((py == NULL) || (PyArray_SIZE(y) != PyArray_SIZE(x))) {
      PyErr_Format( _pixelizeError, "py is of incorrect type (wanted 1D float)");
      goto _fail;
  }

  d = (PyArrayObject *) PyArray_FromAny(dp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, 0, NULL);
  if ((d == NULL) || (PyArray_SIZE(d) != PyArray_SIZE(x))) {
      PyErr_Format( _pixelizeError, "data is of incorrect type (wanted 1D float)");
      goto _fail;
  }

  dx = (PyArrayObject *) PyArray_FromAny(dxp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, 0, NULL);
  if ((dx == NULL) || (PyArray_SIZE(dx) != PyArray_SIZE(x))) {
      PyErr_Format( _pixelizeError, "dx is of incorrect type (wanted 1D float)");
      goto _fail;
  }
  dy = (PyArrayObject *) PyArray_FromAny(dyp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, 0, NULL);
  if ((dy == NULL) || (PyArray_SIZE(dy) != PyArray_SIZE(x))) {
      PyErr_Format( _pixelizeError, "dy is of incorrect type (wanted 1D float)");
      goto _fail;
  }
  dz = (PyArrayObject *) PyArray_FromAny(dzp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, 0, NULL);
  if ((dz == NULL) || (PyArray_SIZE(dz) != PyArray_SIZE(x))) {
      PyErr_Format( _pixelizeError, "dz is of incorrect type (wanted 1D float)");
      goto _fail;
  }
  center = (PyArrayObject *) PyArray_FromAny(centerp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, NPY_C_CONTIGUOUS, NULL);
  if ((dz == NULL) || (PyArray_SIZE(center) != 3)) {
      PyErr_Format( _pixelizeError, "Center must have three points");
      goto _fail;
  }
  inv_mat = (PyArrayObject *) PyArray_FromAny(inv_matp,
            PyArray_DescrFromType(NPY_FLOAT64), 2, 2, 0, NULL);
  if ((inv_mat == NULL) || (PyArray_SIZE(inv_mat) != 9)) {
      PyErr_Format( _pixelizeError, "inv_mat must be three by three");
      goto _fail;
  }
  indices = (PyArrayObject *) PyArray_FromAny(indicesp,
            PyArray_DescrFromType(NPY_INT64), 1, 1, 0, NULL);
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

  npy_float64 xsp, ysp, zsp, pxsp, pysp, dxsp, dysp, dzsp, dsp;
  npy_float64 *centers = (npy_float64 *) PyArray_GETPTR1(center,0);

  npy_intp dims[] = {rows, cols};
  PyArrayObject *my_array =
    (PyArrayObject *) PyArray_SimpleNewFromDescr(2, dims,
              PyArray_DescrFromType(NPY_FLOAT64));
  npy_float64 *gridded = (npy_float64 *) my_array->data;
  npy_float64 *mask = malloc(sizeof(npy_float64)*rows*cols);

  npy_float64 inv_mats[3][3];
  for(i=0;i<3;i++)for(j=0;j<3;j++)
      inv_mats[i][j]=*(npy_float64*)PyArray_GETPTR2(inv_mat,i,j);

  int pp;
  npy_float64 radius;
  for(p=0;p<cols*rows;p++)gridded[p]=mask[p]=0.0;
  for(pp=0; pp<nx; pp++)
  {
    p = *((npy_int64 *) PyArray_GETPTR1(indices, pp));
    npy_float64 xsp = *((npy_float64 *) PyArray_GETPTR1(x, p));
    npy_float64 ysp = *((npy_float64 *) PyArray_GETPTR1(y, p));
    npy_float64 zsp = *((npy_float64 *) PyArray_GETPTR1(z, p));
    npy_float64 pxsp = *((npy_float64 *) PyArray_GETPTR1(px, p));
    npy_float64 pysp = *((npy_float64 *) PyArray_GETPTR1(py, p));
    npy_float64 dxsp = *((npy_float64 *) PyArray_GETPTR1(dx, p));
    npy_float64 dysp = *((npy_float64 *) PyArray_GETPTR1(dy, p));
    npy_float64 dzsp = *((npy_float64 *) PyArray_GETPTR1(dz, p));
    npy_float64 dsp = *((npy_float64 *) PyArray_GETPTR1(d, p)); // We check this above
    // Any point we want to plot is at most this far from the center
    md = 2.0*sqrtl(dxsp*dxsp + dysp*dysp + dzsp*dzsp);
    if(((pxsp+md<px_min) ||
        (pxsp-md>px_max)) ||
       ((pysp+md<py_min) ||
        (pysp-md>py_max))) continue;
    lc = max(floorl((pxsp-md-px_min)/px_dx),0);
    lr = max(floorl((pysp-md-py_min)/px_dy),0);
    rc = min(ceill((pxsp+md-px_min)/px_dx),rows);
    rr = min(ceill((pysp+md-py_min)/px_dy),cols);
    for (i=lr;i<rr;i++) {
      cypx = px_dy * (i+0.5) + py_min;
      for (j=lc;j<rc;j++) {
        cxpx = px_dx * (j+0.5) + px_min;
        cx = inv_mats[0][0]*cxpx + inv_mats[0][1]*cypx + centers[0];
        cy = inv_mats[1][0]*cxpx + inv_mats[1][1]*cypx + centers[1];
        cz = inv_mats[2][0]*cxpx + inv_mats[2][1]*cypx + centers[2];
        if( (fabs(xsp-cx)*0.95>dxsp) || 
            (fabs(ysp-cy)*0.95>dysp) ||
            (fabs(zsp-cz)*0.95>dzsp)) continue;
        mask[j*cols+i] += 1;
        gridded[j*cols+i] += dsp;
      }
    }
  }
  for(p=0;p<cols*rows;p++)gridded[p]=gridded[p]/mask[p];

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
  free(mask);

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
