//
// AMRPixelize
//   A module for making static-resolution arrays representing
//   AMR data.
//
// Written by: Matthew Turk (mturk@stanford.edu) May 2007
// Modified:
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

static PyObject* Py_Pixelize(PyObject *obj, PyObject *args) {

  PyObject *xp, *yp, *dxp, *dyp, *dp;
  unsigned int rows, cols;
  double x_min, x_max, y_min, y_max;

    if (!PyArg_ParseTuple(args, "OOOOOIIdddd",
        &xp, &yp, &dxp, &dyp, &dp, &rows, &cols, 
        &x_min, &x_max, &y_min, &y_max))
        return PyErr_Format(_pixelizeError, "Pixelize: Invalid Parameters.");

  double width = x_max - x_min;
  double height = y_max - y_min;
  double px_dx = width / ((double) cols);
  double px_dy = height / ((double) rows);

  // Check we have something to output to
  if (rows == 0 || cols ==0)
      PyErr_Format( _pixelizeError, "Cannot scale to zero size.");

  // Get numeric arrays
  PyArrayObject *x = (PyArrayObject *) PyArray_FromAny(xp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, NPY_C_CONTIGUOUS, NULL);
  if (x == NULL)
      PyErr_Format( _pixelizeError, "x is of incorrect type (wanted 1D float)");
  PyArrayObject *y = (PyArrayObject *) PyArray_FromAny(yp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, NPY_C_CONTIGUOUS, NULL);
  if (y == NULL) {
      Py_XDECREF(x);
      PyErr_Format( _pixelizeError, "y is of incorrect type (wanted 1D float)");
  }
  PyArrayObject *d = (PyArrayObject *) PyArray_FromAny(dp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, NPY_C_CONTIGUOUS, NULL);
  if (d == NULL) {
      Py_XDECREF(x);
      Py_XDECREF(y);
      PyErr_Format( _pixelizeError, "data is of incorrect type (wanted 1D float)");
  }

  PyArrayObject *dx = (PyArrayObject *) PyArray_FromAny(dxp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, NPY_C_CONTIGUOUS, NULL);
  if (x == NULL) {
      Py_XDECREF(x);
      Py_XDECREF(y);
      Py_XDECREF(d);
      PyErr_Format( _pixelizeError, "dx is of incorrect type (wanted 1D float)");
  }
  PyArrayObject *dy = (PyArrayObject *) PyArray_FromAny(dyp,
            PyArray_DescrFromType(NPY_FLOAT64), 1, 1, NPY_C_CONTIGUOUS, NULL);
  if (x == NULL) {
      Py_XDECREF(x);
      Py_XDECREF(y);
      Py_XDECREF(d);
      Py_XDECREF(dx);
      PyErr_Format( _pixelizeError, "dy is of incorrect type (wanted 1D float)");
  }

  // Check dimensions match
  int nx = x->dimensions[0];
  int ny = y->dimensions[0];
  int ndx = dx->dimensions[0];
  int ndy = dy->dimensions[0];
  if (nx != d->dimensions[0] || ny != d->dimensions[0] ||
      ndx != d->dimensions[0] || ndy != d->dimensions[0]) {
      Py_XDECREF(x);
      Py_XDECREF(y);
      Py_XDECREF(d);
      Py_XDECREF(dx);
      Py_XDECREF(dy);
      PyErr_Format( _pixelizeError, "data and axis dimensions do not match");
  }

  npy_float64 *gridded;
    gridded = malloc(sizeof(npy_float64)*rows*cols);
  if (gridded == NULL) {
      Py_XDECREF(x);
      Py_XDECREF(y);
      Py_XDECREF(d);
      Py_XDECREF(dx);
      Py_XDECREF(dy);
      PyErr_Format(_pixelizeError, "Could not allocate memory for image");
  }

  // Calculate the pointer arrays to map input x to output x
  int i, j, p;
  double lc, lr, rc, rr;
  double lypx, rypx, lxpx, rxpx, overlap, overlap1, overlap2;
  npy_float64 *xs = (npy_float64 *) PyArray_GETPTR1(x, 0);
  npy_float64 *ys = (npy_float64 *) PyArray_GETPTR1(y, 0);
  npy_float64 *dxs = (npy_float64 *) PyArray_GETPTR1(dx, 0);
  npy_float64 *dys = (npy_float64 *) PyArray_GETPTR1(dy, 0);
  npy_float64 *ds = (npy_float64 *) PyArray_GETPTR1(d, 0); // We check this above

  // Upper left is (center - (dx, dy))/dx_per_pixel
  // Lower right is (center + (dx, dy))/dx_per_pixel
  // X width in pixels is int(dx / (dx_per_pixel))
  // Y width in pixels is int(dy / (dy_per_pixel))
  // Center is at (center - LE)/dx_per_pixel

  for (p=0;p<cols*rows;p++) gridded[p] = 0;
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
    for (i=lr;i<rr;i++)
      for (j=lc;j<rc;j++) {
        lypx = px_dy * i + y_min;
        rypx = px_dy * (i+1) + y_min;
        lxpx = px_dx * j + x_min;
        rxpx = px_dx * (j+1) + x_min;
        overlap1 = ((min(rxpx, xs[p]+dxs[p]) - max(lxpx, (xs[p]-dxs[p])))/px_dx);
        overlap2 = ((min(rypx, ys[p]+dys[p]) - max(lypx, (ys[p]-dys[p])))/px_dy);
        if (overlap1 < 0.0 || overlap2 < 0.0) continue;
        gridded[j*cols+i] += (ds[p]*overlap1)*(overlap2);
      }
  }

  // Attatch output buffer to output buffer

  Py_XDECREF(x);
  Py_XDECREF(y);
  Py_XDECREF(d);
  Py_XDECREF(dx);
  Py_XDECREF(dy);

  int dims[] = {rows, cols};

  PyObject* gridret = PyArray_FromDimsAndData(
        2, dims, NPY_FLOAT64, (char *) gridded);

  return gridret;
}

static PyMethodDef _AMRPixelizeMethods[] = {
    {"Pixelize", Py_Pixelize, METH_VARARGS},
    {NULL, NULL} /* Sentinel */
};

/* platform independent*/
#ifdef MS_WIN32
__declspec(dllexport)
#endif

void initAMRPixelize(void)
{
    PyObject *m, *d;
    m = Py_InitModule("AMRPixelize", _AMRPixelizeMethods);
    d = PyModule_GetDict(m);
    _pixelizeError = PyErr_NewException("AMRPixelize.error", NULL, NULL);
    PyDict_SetItemString(d, "error", _pixelizeError);
    import_array();
}

