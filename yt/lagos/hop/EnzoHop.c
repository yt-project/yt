/************************************************************************
* Copyright (C) 2008-2009 Matthew Turk.  All Rights Reserved.
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
// EnzoHop
//   A module for running HOP halo finding on a set of particles
//

#include "Python.h"
#include <stdio.h>
#include <math.h>
#include <signal.h>
#include <ctype.h>
#include "kd.h"
#include "hop.h"

#include "numpy/ndarrayobject.h"

void initgrouplist(Grouplist *g);
void hop_main(KD kd, HC *my_comm, float densthres);
void regroup_main(float dens_outer, HC *my_comm,
        int mingroupsize, float dsaddle);
static PyObject *_HOPerror;

static PyObject *
Py_EnzoHop(PyObject *obj, PyObject *args)
{
    PyObject    *oxpos, *oypos, *ozpos,
                *omass;

    PyArrayObject    *xpos, *ypos, *zpos,
                     *mass;
    xpos=ypos=zpos=mass=NULL;
    npy_float64 totalmass = 0.0;
    float thresh = 160.0;
    int mingroupsize = -1;
    float dsaddle  = -1;

    int i;

    if (!PyArg_ParseTuple(args, "OOOO|fif",
        &oxpos, &oypos, &ozpos, &omass, &thresh,
        &mingroupsize, &dsaddle))
    return PyErr_Format(_HOPerror,
            "EnzoHop: Invalid parameters.");
    if (dsaddle == -1) dsaddle = 2.5 * thresh;

    /* First the regular source arrays */

    xpos    = (PyArrayObject *) PyArray_FromAny(oxpos,
                    PyArray_DescrFromType(NPY_FLOAT64), 1, 1,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if(!xpos){
    PyErr_Format(_HOPerror,
             "EnzoHop: xpos didn't work.");
    goto _fail;
    }
    int num_particles = PyArray_SIZE(xpos);

    ypos    = (PyArrayObject *) PyArray_FromAny(oypos,
                    PyArray_DescrFromType(NPY_FLOAT64), 1, 1,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if((!ypos)||(PyArray_SIZE(ypos) != num_particles)) {
    PyErr_Format(_HOPerror,
             "EnzoHop: xpos and ypos must be the same length.");
    goto _fail;
    }

    zpos    = (PyArrayObject *) PyArray_FromAny(ozpos,
                    PyArray_DescrFromType(NPY_FLOAT64), 1, 1,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if((!zpos)||(PyArray_SIZE(zpos) != num_particles)) {
    PyErr_Format(_HOPerror,
             "EnzoHop: xpos and zpos must be the same length.");
    goto _fail;
    }

    mass    = (PyArrayObject *) PyArray_FromAny(omass,
                    PyArray_DescrFromType(NPY_FLOAT64), 1, 1,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if((!mass)||(PyArray_SIZE(mass) != num_particles)) {
    PyErr_Format(_HOPerror,
             "EnzoHop: xpos and mass must be the same length.");
    goto _fail;
    }

    for(i = 0; i < num_particles; i++)
        totalmass+=*(npy_float64*)PyArray_GETPTR1(mass,i);

  /* initialize the kd hop structure */

  KD kd;
  int nBucket = 16, kdcount = 0;
  kdInit(&kd, nBucket);
  kd->nActive = num_particles;
  kd->p = malloc(sizeof(PARTICLE)*num_particles);
  if (kd->p == NULL) {
    fprintf(stderr, "failed allocating particles.\n");
    goto _fail;
  }
  
 	/* Copy positions into kd structure. */

    fprintf(stdout, "Filling in %d particles\n", num_particles);
	for (i = 0; i < num_particles; i++) {
	  kd->p[i].r[0] = (float)(*(npy_float64*) PyArray_GETPTR1(xpos, i));
	  kd->p[i].r[1] = (float)(*(npy_float64*) PyArray_GETPTR1(ypos, i));
	  kd->p[i].r[2] = (float)(*(npy_float64*) PyArray_GETPTR1(zpos, i));
	  kd->p[i].fMass = (float)(*(npy_float64*) PyArray_GETPTR1(mass, i)/totalmass);
	}

    HC my_comm;
    my_comm.s = newslice();
    my_comm.gl = (Grouplist*)malloc(sizeof(Grouplist));
    if(my_comm.gl == NULL) {
        fprintf(stderr, "failed allocating Grouplist\n");
        goto _fail;
    }
    initgrouplist(my_comm.gl);

    fprintf(stderr, "Calling hop... %d %0.3e\n",num_particles,thresh);
    hop_main(kd, &my_comm, thresh);

    fprintf(stderr, "Calling regroup...\n");

    regroup_main(thresh, &my_comm, mingroupsize, dsaddle);

    // Now we need to get the groupID, realID and the density.
    // This will give us the index into the original array.
    // Additionally, note that we don't really need to tie the index
    // back to the ID in this code, as we can do that back in the python code.
    // All we need to do is provide density and group information.
    
    // Tags (as per writetagsf77) are in gl.s->ntag+1 and there are gl.s->numlist of them.
    PyArrayObject *particle_density = (PyArrayObject *)
            PyArray_SimpleNewFromDescr(1, PyArray_DIMS(xpos),
                    PyArray_DescrFromType(NPY_FLOAT64));
    PyArrayObject *particle_group_id = (PyArrayObject *)
            PyArray_SimpleNewFromDescr(1, PyArray_DIMS(xpos),
                    PyArray_DescrFromType(NPY_INT32));
    
    for (i = 0; i < num_particles; i++) {
      // Density is in kd->p[i].fDensity
      *(npy_float64*)(PyArray_GETPTR1(particle_density, i)) =
            (npy_float64) kd->p[i].fDensity;
      // tag is in gl.s->ntag[i+1]
      *(npy_int32*)(PyArray_GETPTR1(particle_group_id, i)) =
            (npy_int32) my_comm.s->ntag[i+1];
    }

	kdFinish(kd);
    free(my_comm.gl);
    free_slice(my_comm.s);

    PyArray_UpdateFlags(particle_density, NPY_OWNDATA | particle_density->flags);
    PyArray_UpdateFlags(particle_group_id, NPY_OWNDATA | particle_group_id->flags);
    PyObject *return_value = Py_BuildValue("NN", particle_density, particle_group_id);

    Py_DECREF(xpos);
    Py_DECREF(ypos);
    Py_DECREF(zpos);
    Py_DECREF(mass);

    /* We don't need this, as it's done in kdFinish
    if(kd->p!=NULL)free(kd->p);
    */

    return return_value;

_fail:
    Py_XDECREF(xpos);
    Py_XDECREF(ypos);
    Py_XDECREF(zpos);
    Py_XDECREF(mass);

    if(kd->p!=NULL)free(kd->p);

    return NULL;

}

static PyMethodDef _HOPMethods[] = {
    {"RunHOP", Py_EnzoHop, METH_VARARGS},
    {NULL, NULL} /* Sentinel */
};

/* platform independent*/
#ifdef MS_WIN32
__declspec(dllexport)
#endif

void initEnzoHop(void)
{
    PyObject *m, *d;
    m = Py_InitModule("EnzoHop", _HOPMethods);
    d = PyModule_GetDict(m);
    _HOPerror = PyErr_NewException("EnzoHop.HOPerror", NULL, NULL);
    PyDict_SetItemString(d, "error", _HOPerror);
    import_array();
}

/*
 * Local Variables:
 * mode: C
 * c-file-style: "python"
 * End:
 */
