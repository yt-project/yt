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
#include "structmember.h"
#include <stdio.h>
#include <math.h>
#include <signal.h>
#include <ctype.h>
#include "kd.h"
#include "hop.h"
#include "hop_numpy.h"

#include "numpy/ndarrayobject.h"

void initgrouplist(Grouplist *g);
void hop_main(KD kd, HC *my_comm, float densthres);
void regroup_main(float dens_outer, HC *my_comm);
static PyObject *_HOPerror;

int convert_particle_arrays(
    PyObject *oxpos, PyObject *oypos, PyObject *ozpos, PyObject *omass,
    PyArrayObject **xpos, PyArrayObject **ypos, PyArrayObject **zpos,
      PyArrayObject **mass)
{

    /* First the regular source arrays */

    *xpos    = (PyArrayObject *) PyArray_FromAny(oxpos,
                    PyArray_DescrFromType(NPY_FLOAT64), 1, 1,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if(!*xpos){
    PyErr_Format(_HOPerror,
             "EnzoHop: xpos didn't work.");
    return -1;
    }
    int num_particles = PyArray_SIZE(*xpos);

    *ypos    = (PyArrayObject *) PyArray_FromAny(oypos,
                    PyArray_DescrFromType(NPY_FLOAT64), 1, 1,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if((!*ypos)||(PyArray_SIZE(*ypos) != num_particles)) {
    PyErr_Format(_HOPerror,
             "EnzoHop: xpos and ypos must be the same length.");
    return -1;
    }

    *zpos    = (PyArrayObject *) PyArray_FromAny(ozpos,
                    PyArray_DescrFromType(NPY_FLOAT64), 1, 1,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if((!*zpos)||(PyArray_SIZE(*zpos) != num_particles)) {
    PyErr_Format(_HOPerror,
             "EnzoHop: xpos and zpos must be the same length.");
    return -1;
    }

    *mass    = (PyArrayObject *) PyArray_FromAny(omass,
                    PyArray_DescrFromType(NPY_FLOAT64), 1, 1,
                    NPY_INOUT_ARRAY | NPY_UPDATEIFCOPY, NULL);
    if((!*mass)||(PyArray_SIZE(*mass) != num_particles)) {
    PyErr_Format(_HOPerror,
             "EnzoHop: xpos and mass must be the same length.");
    return -1;
    }

    return num_particles;

}
    

static PyObject *
Py_EnzoHop(PyObject *obj, PyObject *args)
{
    PyObject    *oxpos, *oypos, *ozpos,
                *omass;

    PyArrayObject    *xpos, *ypos, *zpos,
                     *mass;
    xpos=ypos=zpos=mass=NULL;
    npy_float64 totalmass = 0.0;
    float normalize_to = 1.0;
    float thresh = 160.0;

    int i;

    if (!PyArg_ParseTuple(args, "OOOO|ff",
        &oxpos, &oypos, &ozpos, &omass, &thresh, &normalize_to))
    return PyErr_Format(_HOPerror,
            "EnzoHop: Invalid parameters.");

    int num_particles = convert_particle_arrays(
            oxpos, oypos, ozpos, omass,
            &xpos, &ypos, &zpos, &mass);
    if (num_particles < 0) goto _fail;

    for(i = 0; i < num_particles; i++)
        totalmass+=*(npy_float64*)PyArray_GETPTR1(mass,i);
    totalmass /= normalize_to;

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
    PyArrayObject *particle_density = (PyArrayObject *)
            PyArray_SimpleNewFromDescr(1, PyArray_DIMS(xpos),
                    PyArray_DescrFromType(NPY_FLOAT64));

    fprintf(stdout, "Copying arrays for %d particles\n", num_particles);
    kd->np_masses = mass;
    kd->np_pos[0] = xpos;
    kd->np_pos[1] = ypos;
    kd->np_pos[2] = zpos;
    kd->np_densities = particle_density;
    kd->totalmass = totalmass;
	for (i = 0; i < num_particles; i++) kd->p[i].np_index = i;

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
    regroup_main(thresh, &my_comm);

    // Now we need to get the groupID, realID and the density.
    // This will give us the index into the original array.
    // Additionally, note that we don't really need to tie the index
    // back to the ID in this code, as we can do that back in the python code.
    // All we need to do is provide density and group information.
    
    // Tags (as per writetagsf77) are in gl.s->ntag+1 and there are gl.s->numlist of them.
    PyArrayObject *particle_group_id = (PyArrayObject *)
            PyArray_SimpleNewFromDescr(1, PyArray_DIMS(xpos),
                    PyArray_DescrFromType(NPY_INT32));
    
    for (i = 0; i < num_particles; i++) {
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

//
// Now a fun wrapper class for the kD-tree
//

typedef struct {
    PyObject_HEAD
    KD kd;
    PyArrayObject *xpos, *ypos, *zpos;
    PyArrayObject *mass;
    int num_particles;
} kDTreeType;

static int
kDTreeType_init(kDTreeType *self, PyObject *args, PyObject *kwds)
{
    int nBuckets = 16;
    static char *kwlist[] = {"xpos", "ypos", "zpos", "mass",
                             "nbuckets", NULL};
    PyObject    *oxpos, *oypos, *ozpos,
                *omass;
    self->xpos=self->ypos=self->zpos=self->mass=NULL;

    /* 16 buckets, I s'pose */
    kdInit(&self->kd, 16);

    fprintf(stderr, "1\n");
    if (! PyArg_ParseTupleAndKeywords(args, kwds, "OOOO|i", kwlist, 
                           &oxpos, &oypos, &ozpos, &omass,
                           &nBuckets))
        return -1;  /* Should this give an error? */

    fprintf(stderr, "2\n");
    self->num_particles = convert_particle_arrays(
            oxpos, oypos, ozpos, omass,
            &self->xpos, &self->ypos, &self->zpos, &self->mass);
    self->kd->nActive = self->num_particles;
    self->kd->p = malloc(sizeof(PARTICLE)*self->num_particles);
    if (self->kd->p == NULL) {
      fprintf(stderr, "failed allocating particles.\n");
      goto _fail;
    }

    return 0;

    _fail:
        Py_XDECREF(self->xpos);
        Py_XDECREF(self->ypos);
        Py_XDECREF(self->zpos);
        Py_XDECREF(self->mass);

        if(self->kd->p!=NULL)free(self->kd->p);

        return -1;
}

static void
kDTreeType_dealloc(kDTreeType *self)
{
   kdFinish(self->kd);
   Py_XDECREF(self->xpos);
   Py_XDECREF(self->ypos);
   Py_XDECREF(self->zpos);
   Py_XDECREF(self->mass);

   self->ob_type->tp_free((PyObject*)self);
}

static PyMemberDef kDTreeType_members[] = {
   { "xpos",  T_OBJECT,    offsetof(kDTreeType, xpos), 0,
               "The xposition array."},
   { "ypos",  T_OBJECT,    offsetof(kDTreeType, ypos), 0,
               "The yposition array."},
   { "zpos",  T_OBJECT,    offsetof(kDTreeType, zpos), 0,
               "The zposition array."},
   { "mass",  T_OBJECT,    offsetof(kDTreeType, mass), 0,
               "The mass array."},
   { "nparticles", T_INT,  offsetof(kDTreeType, num_particles), 0,
               "The number of particles"},
   { NULL }
};

static PyMethodDef
kDTreeType_methods[] = {
   /*{ "set",    (PyCFunction) CountDict_set, METH_VARARGS,
               "Set a key and increment the count." },*/
   // typically there would be more here...
   
   { NULL }
};

static PyTypeObject
kDTreeTypeDict = {
   PyObject_HEAD_INIT(NULL)
   0,                         /* ob_size */
   "kDTree",               /* tp_name */
   sizeof(kDTreeType),         /* tp_basicsize */
   0,                         /* tp_itemsize */
   (destructor)kDTreeType_dealloc, /* tp_dealloc */
   0,                         /* tp_print */
   0,                         /* tp_getattr */
   0,                         /* tp_setattr */
   0,                         /* tp_compare */
   0,                         /* tp_repr */
   0,                         /* tp_as_number */
   0,                         /* tp_as_sequence */
   0,                         /* tp_as_mapping */
   0,                         /* tp_hash */
   0,                         /* tp_call */
   0,                         /* tp_str */
   0,                         /* tp_getattro */
   0,                         /* tp_setattro */
   0,                         /* tp_as_buffer */
   Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /* tp_flags*/
   "kDTree object",           /* tp_doc */
   0,                         /* tp_traverse */
   0,                         /* tp_clear */
   0,                         /* tp_richcompare */
   0,                         /* tp_weaklistoffset */
   0,                         /* tp_iter */
   0,                         /* tp_iternext */
   kDTreeType_methods,         /* tp_methods */
   kDTreeType_members,         /* tp_members */
   0,                         /* tp_getset */
   0,                         /* tp_base */
   0,                         /* tp_dict */
   0,                         /* tp_descr_get */
   0,                         /* tp_descr_set */
   0,                         /* tp_dictoffset */
   (initproc)kDTreeType_init,     /* tp_init */
   0,                         /* tp_alloc */
   0,                         /* tp_new */
};

void initEnzoHop(void)
{
    PyObject *m, *d;
    m = Py_InitModule("EnzoHop", _HOPMethods);
    d = PyModule_GetDict(m);
    _HOPerror = PyErr_NewException("EnzoHop.HOPerror", NULL, NULL);
    PyDict_SetItemString(d, "error", _HOPerror);

    kDTreeTypeDict.tp_new = PyType_GenericNew;
    if (PyType_Ready(&kDTreeTypeDict) < 0) {
       return;
    }

   Py_INCREF(&kDTreeTypeDict);
   PyModule_AddObject(m, "kDTree", (PyObject*)&kDTreeTypeDict);

   import_array();
}

/*
 * Local Variables:
 * mode: C
 * c-file-style: "python"
 * End:
 */
