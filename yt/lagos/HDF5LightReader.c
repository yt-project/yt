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
// HDF5_LightReader
//   A module for light-weight reading of HDF5 files
//

#include "Python.h"

#include <stdio.h>
#include <math.h>
#include <signal.h>
#include <ctype.h>
#include "H5LT.h"
#include "hdf5.h"

#include "numpy/ndarrayobject.h"


static PyObject *_hdf5ReadError;

static PyObject *
Py_ReadHDF5DataSet(PyObject *obj, PyObject *args)
{
    char *filename, *nodename;

    hsize_t *my_dims;
    hid_t file_id;
    herr_t my_error;
    H5T_class_t class_id;
    size_t type_size;
    int my_typenum, my_rank, i, my_size;
    void *my_data;
    PyArrayObject *my_array;

    if (!PyArg_ParseTuple(args, "ss",
            &filename, &nodename))
        return PyErr_Format(_hdf5ReadError,
               "ReadHDF5DataSet: Invalid parameters.");

    //char* filename = *ofilename;
    //char* nodename = *onodename;

    file_id = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT); 
    
    my_error = H5LTget_dataset_ndims ( file_id, nodename, &my_rank );
    if(my_error) {
        PyErr_Format(_hdf5ReadError,
                 "ReadHDF5DataSet: Problem getting dataset info (%s, %s).",
                                    filename, nodename);
        goto _fail;
    }

    my_dims = malloc(sizeof(hsize_t) * my_rank);
    my_error = H5LTget_dataset_info ( file_id, nodename,
                my_dims, &class_id, &type_size );
    if(my_error) {
        PyErr_Format(_hdf5ReadError,
                 "ReadHDF5DataSet: Problem getting dataset info (%s, %s).",
                                    filename, nodename);
        goto _fail;
    }

    my_size = 1;
    npy_intp *dims = malloc(my_rank * sizeof(npy_intp));
    for (i = 0; i < my_rank; i++) {
      dims[i] = (npy_intp) my_dims[i];
      my_size *= my_dims[i];
    }

    if (!(class_id == H5T_FLOAT)){
      PyErr_Format(_hdf5ReadError,
          "ReadHDF5DataSet: Unrecognized datatype, size %i.", type_size);
      goto _fail;
    }

    /*
    switch (type_size) {
      case 4:
        fprintf(stderr, "Reading (%i) %i\n", my_size, type_size);
        my_typenum = NPY_FLOAT32;
        H5LTread_dataset_float(file_id, nodename, my_data);
        break;
      case 8:
        fprintf(stderr, "Reading (%i) %i\n", my_size, type_size);
        my_typenum = NPY_FLOAT64;
        H5LTread_dataset_double(file_id, nodename, my_data);
        break;
      default:
        PyErr_Format(_hdf5ReadError,
            "ReadHDF5DataSet: Unrecognized datatype, size %i.", type_size);
        goto _fail;
        break; //haha goto!
    }
    */
    my_array = (PyArrayObject *) PyArray_SimpleNewFromDescr(my_rank, dims,
                PyArray_DescrFromType(NPY_DOUBLE));

    H5LTread_dataset_double(file_id, nodename, (void *) my_array->data);
    H5Fclose(file_id);
    /*
    my_array = (PyArrayObject *) PyArray_SimpleNewFromData(my_rank, dims,
                                    NPY_FLOAT64, (void *)my_data);
    */

    PyArray_UpdateFlags(my_array, NPY_OWNDATA | my_array->flags);
    // 'N' does not increase the reference count
    PyObject *return_value = Py_BuildValue("N", my_array);

    free(dims);
    free(my_dims);

    return return_value;

    _fail:
      Py_XDECREF(my_array);
      if(file_id) H5Fclose(file_id);
      if(my_data) free(my_data);
      if(my_dims) free(my_dims);
      if(dims) free(dims);
      return NULL;
}

static PyMethodDef _hdf5LightReaderMethods[] = {
    {"ReadData", Py_ReadHDF5DataSet, METH_VARARGS},
    {NULL, NULL} /* Sentinel */
};

/* platform independent*/
#ifdef MS_WIN32
__declspec(dllexport)
#endif

void initHDF5LightReader(void)
{
    PyObject *m, *d;
    m = Py_InitModule("HDF5LightReader", _hdf5LightReaderMethods);
    d = PyModule_GetDict(m);
    _hdf5ReadError = PyErr_NewException("HDF5LightReader.ReadingError", NULL, NULL);
    PyDict_SetItemString(d, "error", _hdf5ReadError);
    import_array();
}

/*
 * Local Variables:
 * mode: C
 * c-file-style: "python"
 * End:
 */
