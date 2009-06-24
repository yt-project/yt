#ifndef _NUMPY_HOP_H
#include "Python.h"
#include "numpy/ndarrayobject.h"

#define NP_DENS(kd, in) \
    (*(npy_float64*)PyArray_GETPTR1(kd->np_densities, kd->p[in].np_index))
#define NP_POS(kd, in, dim) \
    (*(npy_float64*)PyArray_GETPTR1(kd->np_pos[dim], kd->p[in].np_index))
#define NP_MASS(kd, in) \
    (*(npy_float64*)PyArray_GETPTR1(kd->np_masses, kd->p[in].np_index))/kd->totalmass

#endif
