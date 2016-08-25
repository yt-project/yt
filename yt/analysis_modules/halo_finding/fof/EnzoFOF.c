/*******************************************************************************
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
*******************************************************************************/

//
// EnzoFOF
//   A module for running friends-of-friends halo finding on a set of particles
//

#include "Python.h"
#include <stdio.h>
#include <math.h>
#include <signal.h>
#include <ctype.h>
#include "kd.h"
#include "tipsydefs.h"

#include "numpy/ndarrayobject.h"


static PyObject *_FOFerror;

static PyObject *
Py_EnzoFOF(PyObject *obj, PyObject *args)
{
    PyObject    *oxpos, *oypos, *ozpos;
    PyArrayObject    *xpos, *ypos, *zpos;
    float link = 0.2;
    float fPeriod[3] = {1.0, 1.0, 1.0};
	int nMembers = 8;
    int i, num_particles;
	KDFOF kd;
	int nBucket,j;
	float fEps;
	int nGroup,bVerbose=1;
	int sec,usec;
	PyArrayObject *particle_group_id;
    PyObject *return_value;

    xpos=ypos=zpos=NULL;

    if (!PyArg_ParseTuple(args, "OOO|f(fff)i",
        &oxpos, &oypos, &ozpos, &link,
        &fPeriod[0], &fPeriod[1], &fPeriod[2],
        &nMembers))
    return PyErr_Format(_FOFerror,
            "EnzoFOF: Invalid parameters.");

    /* First the regular source arrays */

    xpos    = (PyArrayObject *) PyArray_FromAny(oxpos,
                    PyArray_DescrFromType(NPY_FLOAT64), 1, 1,
                    NPY_ARRAY_INOUT_ARRAY | NPY_ARRAY_UPDATEIFCOPY, NULL);
    if(!xpos){
    PyErr_Format(_FOFerror,
             "EnzoFOF: xpos didn't work.");
    goto _fail;
    }
    num_particles = PyArray_SIZE(xpos);

    ypos    = (PyArrayObject *) PyArray_FromAny(oypos,
                    PyArray_DescrFromType(NPY_FLOAT64), 1, 1,
                    NPY_ARRAY_INOUT_ARRAY | NPY_ARRAY_UPDATEIFCOPY, NULL);
    if((!ypos)||(PyArray_SIZE(ypos) != num_particles)) {
    PyErr_Format(_FOFerror,
             "EnzoFOF: xpos and ypos must be the same length.");
    goto _fail;
    }

    zpos    = (PyArrayObject *) PyArray_FromAny(ozpos,
                    PyArray_DescrFromType(NPY_FLOAT64), 1, 1,
                    NPY_ARRAY_INOUT_ARRAY | NPY_ARRAY_UPDATEIFCOPY, NULL);
    if((!zpos)||(PyArray_SIZE(zpos) != num_particles)) {
    PyErr_Format(_FOFerror,
             "EnzoFOF: xpos and zpos must be the same length.");
    goto _fail;
    }

    /* let's get started with the FOF stuff */
	
	/* linking length */
	fprintf(stdout, "Link length is %f\n", link);
	fEps = link;
	
	nBucket = 16;

    /* initialize the kd FOF structure */

	kdInitFoF(&kd,nBucket,fPeriod);
	
	/* kdReadTipsyFoF(kd,stdin,bDark,bGas,bStar); */

 	/* Copy positions into kd structure. */

    fprintf(stdout, "Filling in %d particles\n", num_particles);
    kd->nActive = num_particles;
	kd->p = (PARTICLEFOF *)malloc(kd->nActive*sizeof(PARTICLEFOF));
	assert(kd->p != NULL);
	for (i = 0; i < num_particles; i++) {
	  kd->p[i].iOrder = i;
	  kd->p[i].r[0] = (float)(*(npy_float64*) PyArray_GETPTR1(xpos, i));
	  kd->p[i].r[1] = (float)(*(npy_float64*) PyArray_GETPTR1(ypos, i));
	  kd->p[i].r[2] = (float)(*(npy_float64*) PyArray_GETPTR1(zpos, i));
	}
	
	kdBuildTreeFoF(kd);
	kdTimeFoF(kd,&sec,&usec);
	nGroup = kdFoF(kd,fEps);
	kdTimeFoF(kd,&sec,&usec);
	if (bVerbose) printf("Number of initial groups:%d\n",nGroup);
	nGroup = kdTooSmallFoF(kd,nMembers);
	if (bVerbose) {
		printf("Number of groups:%d\n",nGroup);
		printf("FOF CPU TIME: %d.%06d secs\n",sec,usec);
		}
	kdOrderFoF(kd);

	/* kdOutGroupFoF(kd,ach); */
	
    // Now we need to get the groupID, realID.
    // This will give us the index into the original array.
    // Additionally, note that we don't really need to tie the index
    // back to the ID in this code, as we can do that back in the python code.
    // All we need to do is group information.
    
    // Tags are in kd->p[i].iGroup
    particle_group_id = (PyArrayObject *)
            PyArray_SimpleNewFromDescr(1, PyArray_DIMS(xpos),
                    PyArray_DescrFromType(NPY_INT32));
    
    for (i = 0; i < num_particles; i++) {
      // group tag is in kd->p[i].iGroup
      *(npy_int32*)(PyArray_GETPTR1(particle_group_id, i)) =
            (npy_int32) kd->p[i].iGroup;
    }

	kdFinishFoF(kd);

    PyArray_UpdateFlags(particle_group_id,
        NPY_ARRAY_OWNDATA | PyArray_FLAGS(particle_group_id));
    return_value = Py_BuildValue("N", particle_group_id);

    Py_DECREF(xpos);
    Py_DECREF(ypos);
    Py_DECREF(zpos);

    /* We don't need this, as it's done in kdFinish
    if(kd->p!=NULL)free(kd->p);
    */

    return return_value;

_fail:
    Py_XDECREF(xpos);
    Py_XDECREF(ypos);
    Py_XDECREF(zpos);

    if(kd->p!=NULL)free(kd->p);

    return NULL;

}

static PyMethodDef _FOFMethods[] = {
    {"RunFOF", Py_EnzoFOF, METH_VARARGS},
    {NULL, NULL} /* Sentinel */
};

/* platform independent*/
#ifdef MS_WIN32
__declspec(dllexport)
#endif

PyMODINIT_FUNC
#if PY_MAJOR_VERSION >= 3
#define _RETVAL m
PyInit_EnzoFOF(void)
#else
#define _RETVAL 
initEnzoFOF(void)
#endif
{
    PyObject *m, *d;
#if PY_MAJOR_VERSION >= 3
    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "EnzoFOF",           /* m_name */
        "EnzoFOF Module",    /* m_doc */
        -1,                  /* m_size */
        _FOFMethods,          /* m_methods */
        NULL,                /* m_reload */
        NULL,                /* m_traverse */
        NULL,                /* m_clear */
        NULL,                /* m_free */
    };
    m = PyModule_Create(&moduledef); 
#else
    m = Py_InitModule("EnzoFOF", _FOFMethods);
#endif
    d = PyModule_GetDict(m);
    _FOFerror = PyErr_NewException("EnzoFOF.FOFerror", NULL, NULL);
    PyDict_SetItemString(d, "error", _FOFerror);
    import_array();
    return _RETVAL;
}

/*
 * Local Variables:
 * mode: C
 * c-file-style: "python"
 * End:
 */
