/*******************************************************************************
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
*******************************************************************************/

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

#ifndef Py_TYPE
    #define Py_TYPE(ob) (((PyObject*)(ob))->ob_type)
#endif

void initgrouplist(Grouplist *g);
void hop_main(KD kd, HC *my_comm, float densthres);
void regroup_main(float dens_outer, HC *my_comm);
static PyObject *_HOPerror;

int convert_particle_arrays(
    PyObject *oxpos, PyObject *oypos, PyObject *ozpos, PyObject *omass,
    PyArrayObject **xpos, PyArrayObject **ypos, PyArrayObject **zpos,
      PyArrayObject **mass)
{
    int num_particles;

    /* First the regular source arrays */

    *xpos    = (PyArrayObject *) PyArray_FromAny(oxpos,
                    PyArray_DescrFromType(NPY_FLOAT64), 1, 1,
                    NPY_ARRAY_INOUT_ARRAY | NPY_ARRAY_UPDATEIFCOPY, NULL);
    if(!*xpos){
    PyErr_Format(_HOPerror,
             "EnzoHop: xpos didn't work.");
    return -1;
    }
    num_particles = PyArray_SIZE(*xpos);

    *ypos    = (PyArrayObject *) PyArray_FromAny(oypos,
                    PyArray_DescrFromType(NPY_FLOAT64), 1, 1,
                    NPY_ARRAY_INOUT_ARRAY | NPY_ARRAY_UPDATEIFCOPY, NULL);
    if((!*ypos)||(PyArray_SIZE(*ypos) != num_particles)) {
    PyErr_Format(_HOPerror,
             "EnzoHop: xpos and ypos must be the same length.");
    return -1;
    }

    *zpos    = (PyArrayObject *) PyArray_FromAny(ozpos,
                    PyArray_DescrFromType(NPY_FLOAT64), 1, 1,
                    NPY_ARRAY_INOUT_ARRAY | NPY_ARRAY_UPDATEIFCOPY, NULL);
    if((!*zpos)||(PyArray_SIZE(*zpos) != num_particles)) {
    PyErr_Format(_HOPerror,
             "EnzoHop: xpos and zpos must be the same length.");
    return -1;
    }

    *mass    = (PyArrayObject *) PyArray_FromAny(omass,
                    PyArray_DescrFromType(NPY_FLOAT64), 1, 1,
                    NPY_ARRAY_INOUT_ARRAY | NPY_ARRAY_UPDATEIFCOPY, NULL);
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
    npy_float64 totalmass = 0.0;
    float normalize_to = 1.0;
    float thresh = 160.0;
    int i, num_particles;
    KD kd;
    int nBucket = 16, kdcount = 0;
    PyArrayObject *particle_density;
    HC my_comm;
    PyArrayObject *particle_group_id;
    PyObject *return_value;

    xpos=ypos=zpos=mass=NULL;

    if (!PyArg_ParseTuple(args, "OOOO|ff",
        &oxpos, &oypos, &ozpos, &omass, &thresh, &normalize_to))
    return PyErr_Format(_HOPerror,
            "EnzoHop: Invalid parameters.");

    num_particles = convert_particle_arrays(
            oxpos, oypos, ozpos, omass,
            &xpos, &ypos, &zpos, &mass);
    if (num_particles < 0) goto _fail;

    for(i = 0; i < num_particles; i++)
        totalmass+=*(npy_float64*)PyArray_GETPTR1(mass,i);
    totalmass /= normalize_to;

  /* initialize the kd hop structure */

  kdInit(&kd, nBucket);
  kd->nActive = num_particles;
  kd->p = malloc(sizeof(PARTICLE)*num_particles);
  if (kd->p == NULL) {
    fprintf(stderr, "failed allocating particles.\n");
    goto _fail;
  }
  
 	/* Copy positions into kd structure. */
    particle_density = (PyArrayObject *)
            PyArray_SimpleNewFromDescr(1, PyArray_DIMS(xpos),
                    PyArray_DescrFromType(NPY_FLOAT64));

    fprintf(stdout, "Copying arrays for %d particles\n", num_particles);
    kd->np_masses = (npy_float64*) PyArray_DATA(mass);
    kd->np_pos[0] = (npy_float64*) PyArray_DATA(xpos);
    kd->np_pos[1] = (npy_float64*) PyArray_DATA(ypos);
    kd->np_pos[2] = (npy_float64*) PyArray_DATA(zpos);
    kd->np_densities = (npy_float64*) PyArray_DATA(particle_density);
    kd->totalmass = totalmass;
	for (i = 0; i < num_particles; i++) kd->p[i].np_index = i;

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
    particle_group_id = (PyArrayObject *)
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

    PyArray_UpdateFlags(particle_density, NPY_ARRAY_OWNDATA | PyArray_FLAGS(particle_density));
    PyArray_UpdateFlags(particle_group_id, NPY_ARRAY_OWNDATA | PyArray_FLAGS(particle_group_id));
    return_value = Py_BuildValue("NN", particle_density, particle_group_id);

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
    PyArrayObject *mass, *densities;
    int num_particles;
} kDTreeType;

static int
kDTreeType_init(kDTreeType *self, PyObject *args, PyObject *kwds)
{
    int nBuckets = 16, i;
    float normalize_to = 1.0;
    static char *kwlist[] = {"xpos", "ypos", "zpos", "mass",
                             "nbuckets", "norm", NULL};
    PyObject    *oxpos, *oypos, *ozpos,
                *omass;
    npy_float64 totalmass = 0.0;

    self->xpos=self->ypos=self->zpos=self->mass=NULL;


    if (! PyArg_ParseTupleAndKeywords(args, kwds, "OOOO|if", kwlist, 
                           &oxpos, &oypos, &ozpos, &omass,
                           &nBuckets, &normalize_to))
        return -1;  /* Should this give an error? */

    kdInit(&self->kd, nBuckets);

    self->num_particles = convert_particle_arrays(
            oxpos, oypos, ozpos, omass,
            &self->xpos, &self->ypos, &self->zpos, &self->mass);

    self->kd->nActive = self->num_particles;
    self->kd->p = malloc(sizeof(PARTICLE)*self->num_particles);
    if (self->kd->p == NULL) {
      fprintf(stderr, "failed allocating particles.\n");
      goto _fail;
    }

    /* Now we set up our Density array */
    self->densities = (PyArrayObject *)
            PyArray_SimpleNewFromDescr(1, PyArray_DIMS(self->xpos),
                    PyArray_DescrFromType(NPY_FLOAT64));

    for(i= 0; i < self->num_particles; i++) {
        self->kd->p[i].np_index = i;
        *(npy_float64*)(PyArray_GETPTR1(self->densities, i)) = 0.0;
        totalmass+=*(npy_float64*)PyArray_GETPTR1(self->mass,i);
    }
    totalmass /= normalize_to;


    self->kd->np_masses = (npy_float64 *)PyArray_DATA(self->mass);
    self->kd->np_pos[0] = (npy_float64 *)PyArray_DATA(self->xpos);
    self->kd->np_pos[1] = (npy_float64 *)PyArray_DATA(self->ypos);
    self->kd->np_pos[2] = (npy_float64 *)PyArray_DATA(self->zpos);
    self->kd->np_densities = (npy_float64 *)PyArray_DATA(self->densities);
    self->kd->totalmass = totalmass;

    PrepareKD(self->kd);
    kdBuildTree(self->kd);

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

   Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject *
kDTreeType_up_pass(kDTreeType *self, PyObject *args) {
    int iCell;

    if (!PyArg_ParseTuple(args, "i", &iCell))
        return PyErr_Format(_HOPerror,
            "kDTree.up_pass: invalid parameters.");

    if(iCell >= self->num_particles)
        return PyErr_Format(_HOPerror,
            "kDTree.up_pass: iCell cannot be >= num_particles!");

    kdUpPass(self->kd, iCell);
    return Py_None;
}

static PyObject *
kDTreeType_median_jst(kDTreeType *self, PyObject *args) {
    int d, l, u, median;
    PyObject *omedian;

    if (!PyArg_ParseTuple(args, "iii", &d, &l, &u))
        return PyErr_Format(_HOPerror,
            "kDTree.median_jst: invalid parameters.");

    if(d >= 3)
        return PyErr_Format(_HOPerror,
            "kDTree.median_jst: d cannot be >= 3!");

    if(l >= self->num_particles)
        return PyErr_Format(_HOPerror,
            "kDTree.median_jst: l cannot be >= num_particles!");

    if(u >= self->num_particles)
        return PyErr_Format(_HOPerror,
            "kDTree.median_jst: u cannot be >= num_particles!");

    median = kdMedianJst(self->kd, d, l, u);

    omedian = PyLong_FromLong((long)median);
    return omedian;
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
   { "densities",  T_OBJECT,    offsetof(kDTreeType, densities), 0,
               "The density array."},
   { "num_particles", T_INT,  offsetof(kDTreeType, num_particles), 0,
               "The number of particles"},
   { NULL }
};

static PyMethodDef
kDTreeType_methods[] = {
   { "up_pass",    (PyCFunction) kDTreeType_up_pass, METH_VARARGS,
               "Pass up something or another, I'm not really sure."},
   { "median_jst",    (PyCFunction) kDTreeType_median_jst, METH_VARARGS,
               "Use the JST Median algorithm on two points along a dimension."},
   // typically there would be more here...
   
   { NULL }
};

static PyTypeObject
kDTreeTypeDict = {
   PyVarObject_HEAD_INIT(NULL, 0)
                            /* ob_size */
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

PyMODINIT_FUNC
#if PY_MAJOR_VERSION >= 3
#define _RETVAL m
PyInit_EnzoHop(void)
#else
#define _RETVAL 
initEnzoHop(void)
#endif
{
    PyObject *m, *d;
#if PY_MAJOR_VERSION >= 3
    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "EnzoHop",           /* m_name */
        "EnzoHop Module",    /* m_doc */
        -1,                  /* m_size */
        _HOPMethods,          /* m_methods */
        NULL,                /* m_reload */
        NULL,                /* m_traverse */
        NULL,                /* m_clear */
        NULL,                /* m_free */
    };
    m = PyModule_Create(&moduledef); 
#else
    m = Py_InitModule("EnzoHop", _HOPMethods);
#endif
    d = PyModule_GetDict(m);
    _HOPerror = PyErr_NewException("EnzoHop.HOPerror", NULL, NULL);
    PyDict_SetItemString(d, "error", _HOPerror);

    kDTreeTypeDict.tp_new = PyType_GenericNew;
    if (PyType_Ready(&kDTreeTypeDict) < 0) {
       return _RETVAL;
    }

   Py_INCREF(&kDTreeTypeDict);
   PyModule_AddObject(m, "kDTree", (PyObject*)&kDTreeTypeDict);

   import_array();
   return _RETVAL;
}

/*
 * Local Variables:
 * mode: C
 * c-file-style: "python"
 * End:
 */
