#ifndef _NUMPY_HOP_H
#include "Python.h"
#include "numpy/ndarrayobject.h"

#define NP_DENS(kd, in) \
    kd->np_densities[kd->p[in].np_index]
#define NP_POS(kd, in, dim) \
    kd->np_pos[dim][kd->p[in].np_index]
#define NP_MASS(kd, in) \
    (kd->np_masses[kd->p[in].np_index]/kd->totalmass)

#endif
