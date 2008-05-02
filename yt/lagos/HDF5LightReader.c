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
herr_t iterate_dataset(hid_t loc_id, const char *name, void *nodelist);

int get_my_desc_type(hid_t native_type_id){

   /*
   http://terra.rice.edu/comp.res/apps/h/hdf5/docs/RM_H5T.html#Datatype-GetNativeType

   hid_t H5Tget_native_type(hid_t type_id, H5T_direction_t direction  ) returns
   from the following list:
        H5T_NATIVE_CHAR         NPY_??
        H5T_NATIVE_SHORT        NPY_SHORT
        H5T_NATIVE_INT          NPY_INT
        H5T_NATIVE_LONG         NPY_LONG
        H5T_NATIVE_LLONG        NPY_LONGLONG

        H5T_NATIVE_UCHAR        NPY_??
        H5T_NATIVE_USHORT       NPY_USHORT
        H5T_NATIVE_UINT         NPY_UINT
        H5T_NATIVE_ULONG        NPY_ULONG
        H5T_NATIVE_ULLONG       NPY_ULONGLONG

        H5T_NATIVE_FLOAT        NPY_FLOAT
        H5T_NATIVE_DOUBLE       NPY_DOUBLE
        H5T_NATIVE_LDOUBLE      NPY_LONGDOUBLE
    */

       if(H5Tequal(native_type_id, H5T_NATIVE_SHORT   ) > 0){return NPY_SHORT;}
  else if(H5Tequal(native_type_id, H5T_NATIVE_INT     ) > 0){return NPY_INT;}
  else if(H5Tequal(native_type_id, H5T_NATIVE_LONG    ) > 0){return NPY_LONG;}
  else if(H5Tequal(native_type_id, H5T_NATIVE_LLONG   ) > 0){return NPY_LONGLONG;}
  else if(H5Tequal(native_type_id, H5T_NATIVE_USHORT  ) > 0){return NPY_USHORT;}
  else if(H5Tequal(native_type_id, H5T_NATIVE_UINT    ) > 0){return NPY_UINT;}
  else if(H5Tequal(native_type_id, H5T_NATIVE_ULONG   ) > 0){return NPY_ULONG;}
  else if(H5Tequal(native_type_id, H5T_NATIVE_ULLONG  ) > 0){return NPY_ULONGLONG;}
  else if(H5Tequal(native_type_id, H5T_NATIVE_FLOAT   ) > 0){return NPY_FLOAT;}
  else if(H5Tequal(native_type_id, H5T_NATIVE_DOUBLE  ) > 0){return NPY_DOUBLE;}
  else if(H5Tequal(native_type_id, H5T_NATIVE_LDOUBLE ) > 0){return NPY_LONGDOUBLE;}
  else {return -1;}

}

static PyObject *
Py_ReadHDF5DataSet(PyObject *obj, PyObject *args)
{
    char *filename, *nodename;

    hsize_t *my_dims = NULL;
    npy_intp *dims = NULL;
    hid_t file_id, datatype_id, native_type_id, dataset;
    herr_t my_error;
    htri_t file_exists;
    H5T_class_t class_id;
    size_t type_size;
    int my_typenum, my_rank, i;
    H5E_auto_t err_func;
    void *err_datastream;
    PyArrayObject *my_array = NULL;
    file_id = datatype_id = native_type_id = dataset = 0;

    if (!PyArg_ParseTuple(args, "ss",
            &filename, &nodename))
        return PyErr_Format(_hdf5ReadError,
               "ReadHDF5DataSet: Invalid parameters.");

    /* How portable is this? */
    if (access(filename, R_OK) < 0) {
        PyErr_Format(_hdf5ReadError,
                 "ReadHDF5DataSet: %s does not exist, or no read permissions\n",
                     filename);
        goto _fail;
    }

    file_exists = H5Fis_hdf5(filename);
    if (file_exists == 0) {
        PyErr_Format(_hdf5ReadError,
                 "ReadHDF5DataSet: %s is not an HDF5 file", filename);
        goto _fail;
    }

    file_id = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT); 

    if (file_id < 0) {
        PyErr_Format(_hdf5ReadError,
                 "ReadHDF5DataSet: Unable to open %s", nodename);
        goto _fail;
    }

    /* We turn off error reporting briefly, because it turns out that
       reading datasets with group names is more forgiving than finding
       datasets with group names using the high-level interface. */

    H5Eget_auto(&err_func, &err_datastream);
    H5Eset_auto(NULL, NULL);
    dataset = H5Dopen(file_id, nodename);
    H5Eset_auto(err_func, err_datastream);

    if(dataset < 0) {
        PyErr_Format(_hdf5ReadError,
                 "ReadHDF5DataSet: Unable to open dataset (%s, %s).",
                                    filename, nodename);
        goto _fail;
    }

    my_error = H5LTget_dataset_ndims ( file_id, nodename, &my_rank );
    if(my_error) {
        PyErr_Format(_hdf5ReadError,
                 "ReadHDF5DataSet: Problem getting dataset rank (%s, %s).",
                                    filename, nodename);
        goto _fail;
    }

    /* How do we keep this from leaking in failures? */
    my_dims = malloc(sizeof(hsize_t) * my_rank);
    my_error = H5LTget_dataset_info ( file_id, nodename,
                my_dims, &class_id, &type_size );

    if(my_error < 0) {
        PyErr_Format(_hdf5ReadError,
                 "ReadHDF5DataSet: Problem getting dataset info (%s, %s).",
                                    filename, nodename);
        goto _fail;
    }

    dims = malloc(my_rank * sizeof(npy_intp));
    for (i = 0; i < my_rank; i++) dims[i] = (npy_intp) my_dims[i];

    datatype_id = H5Dget_type(dataset);
    native_type_id = H5Tget_native_type(datatype_id, H5T_DIR_ASCEND);

    /* Behavior here is intentionally undefined for non-native types */

    int my_desc_type = get_my_desc_type(native_type_id);
    if (my_desc_type == -1) {
          PyErr_Format(_hdf5ReadError,
                       "ReadHDF5DataSet: Unrecognized datatype.  Use a more advanced reader.");
          goto _fail;
    }

    my_array = (PyArrayObject *) PyArray_SimpleNewFromDescr(my_rank, dims,
                PyArray_DescrFromType(my_desc_type));
    if (!my_array) {
        PyErr_Format(_hdf5ReadError,
                 "ReadHDF5DataSet: Unable to create NumPy array.");
        goto _fail;
    }

    H5LTread_dataset(file_id, nodename, native_type_id, (void *) my_array->data);

    PyArray_UpdateFlags(my_array, NPY_OWNDATA | my_array->flags);
    PyObject *return_value = Py_BuildValue("N", my_array);

    H5Fclose(file_id);
    H5Dclose(dataset);
    H5Tclose(native_type_id);
    H5Tclose(datatype_id);
    free(my_dims);
    free(dims);

    return return_value;

    _fail:
      Py_XDECREF(my_array);
      if(!(file_id <= 0)&&(H5Iget_ref(file_id))) H5Fclose(file_id);
      if(!(dataset <= 0)&&(H5Iget_ref(dataset))) H5Dclose(dataset);
      if(!(native_type_id <= 0)&&(H5Iget_ref(native_type_id))) H5Tclose(native_type_id);
      if(!(datatype_id <= 0)&&(H5Iget_ref(datatype_id))) H5Tclose(datatype_id);
      if(my_dims != NULL) free(my_dims);
      if(dims != NULL) free(dims);
      return NULL;
}

static PyObject *
Py_ReadHDF5DataSetSlice(PyObject *obj, PyObject *args)
{
    char *filename, *nodename;

    hsize_t *my_dims = NULL;
    npy_intp *dims = NULL;
    hid_t file_id, datatype_id, native_type_id, dataset, dataspace, memspace;
    herr_t my_error;
    htri_t file_exists;
    H5T_class_t class_id;
    size_t type_size;
    int my_typenum, my_rank, i, axis, coord;
    H5E_auto_t err_func;
    void *err_datastream;
    PyArrayObject *my_array = NULL;
    file_id = datatype_id = native_type_id = dataset = dataspace = memspace = 0;

    if (!PyArg_ParseTuple(args, "ssII",
            &filename, &nodename, &axis, &coord))
        return PyErr_Format(_hdf5ReadError,
               "ReadHDF5DataSetSlice: Invalid parameters.");

    /* How portable is this? */
    if (access(filename, R_OK) < 0) {
        PyErr_Format(_hdf5ReadError,
                 "ReadHDF5DataSetSlice: %s does not exist, or no read permissions\n",
                     filename);
        goto _fail;
    }

    file_exists = H5Fis_hdf5(filename);
    if (file_exists == 0) {
        PyErr_Format(_hdf5ReadError,
                 "ReadHDF5DataSetSlice: %s is not an HDF5 file", filename);
        goto _fail;
    }

    file_id = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT); 

    if (file_id < 0) {
        PyErr_Format(_hdf5ReadError,
                 "ReadHDF5DataSetSlice: Unable to open %s", nodename);
        goto _fail;
    }

    /* We turn off error reporting briefly, because it turns out that
       reading datasets with group names is more forgiving than finding
       datasets with group names using the high-level interface. */

    H5Eget_auto(&err_func, &err_datastream);
    H5Eset_auto(NULL, NULL);
    dataset = H5Dopen(file_id, nodename);
    H5Eset_auto(err_func, err_datastream);

    if(dataset < 0) {
        PyErr_Format(_hdf5ReadError,
                 "ReadHDF5DataSetSlice: Unable to open dataset (%s, %s).",
                                    filename, nodename);
        goto _fail;
    }

    my_error = H5LTget_dataset_ndims ( file_id, nodename, &my_rank );
    if(my_error) {
        PyErr_Format(_hdf5ReadError,
                 "ReadHDF5DataSetSlice: Problem getting dataset rank (%s, %s).",
                                    filename, nodename);
        goto _fail;
    }
    if(my_rank!=3) {
        PyErr_Format(_hdf5ReadError,
                 "ReadHDF5DataSetSlice: Sorry, I only know how to slice 3D into 2D.");
        goto _fail;
    }

    /* How do we keep this from leaking in failures? */
    my_dims = malloc(sizeof(hsize_t) * my_rank);
    my_error = H5LTget_dataset_info ( file_id, nodename,
                my_dims, &class_id, &type_size );

    if(my_error < 0) {
        PyErr_Format(_hdf5ReadError,
                 "ReadHDF5DataSetSlice: Problem getting dataset info (%s, %s).",
                                    filename, nodename);
        goto _fail;
    }

    dims = malloc(my_rank * sizeof(npy_intp));
    for (i = 0; i < my_rank; i++) dims[i] = (npy_intp) my_dims[i];

    datatype_id = H5Dget_type(dataset);
    native_type_id = H5Tget_native_type(datatype_id, H5T_DIR_ASCEND);

    /* Behavior here is intentionally undefined for non-native types */

    int my_desc_type = get_my_desc_type(native_type_id);
    if (my_desc_type == -1) {
          PyErr_Format(_hdf5ReadError,
                       "ReadHDF5DataSetSlice: Unrecognized datatype.  Use a more advanced reader.");
          goto _fail;
    }

    /* Okay, now let's figure out what the dataspaces will look like */

    hsize_t slice_coords[3] = {0, 0, 0};
    slice_coords[axis] = coord;
    hsize_t slice_blocks[3] = {dims[0], dims[1], dims[2]};
    slice_blocks[axis] = 1;

    dataspace = H5Dget_space(dataset);
    if(my_error) {
        PyErr_Format(_hdf5ReadError,
                 "ReadHDF5DataSetSlice: Problem getting datasspace.");
        goto _fail;
    }

    my_error = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, slice_coords, NULL, 
     slice_blocks, NULL);
    if(my_error) {
        PyErr_Format(_hdf5ReadError,
                 "ReadHDF5DataSetSlice: Problem selecting hyperslab.");
        goto _fail;
    }

    hsize_t slice_dims[2];
    int j = 0;
    for (i=0;i<3;i++)if(i!=axis)slice_dims[j++]=dims[i];
    memspace = H5Screate_simple(2,slice_dims,NULL); 

    npy_intp slice_dims_i[2] = {slice_dims[0], slice_dims[1]};
    my_array = (PyArrayObject *) PyArray_SimpleNewFromDescr(2, slice_dims_i,
                PyArray_DescrFromType(my_desc_type));
    if (!my_array) {
        PyErr_Format(_hdf5ReadError,
                 "ReadHDF5DataSetSlice: Unable to create NumPy array.");
        goto _fail;
    }

    my_error = H5Dread(dataset, native_type_id, memspace, dataspace, H5P_DEFAULT,    
                       my_array->data);

    PyArray_UpdateFlags(my_array, NPY_OWNDATA | my_array->flags);
    PyObject *return_value = Py_BuildValue("N", my_array);

    H5Fclose(file_id);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Sclose(memspace);
    H5Tclose(native_type_id);
    H5Tclose(datatype_id);
    free(my_dims);
    free(dims);

    return return_value;

    _fail:
      Py_XDECREF(my_array);
      if(!(file_id <= 0)&&(H5Iget_ref(file_id))) H5Fclose(file_id);
      if(!(dataset <= 0)&&(H5Iget_ref(dataset))) H5Dclose(dataset);
      if(!(dataspace <= 0)&&(H5Iget_ref(dataspace))) H5Sclose(dataspace);
      if(!(memspace <= 0)&&(H5Iget_ref(memspace))) H5Sclose(memspace);
      if(!(native_type_id <= 0)&&(H5Iget_ref(native_type_id))) H5Tclose(native_type_id);
      if(!(datatype_id <= 0)&&(H5Iget_ref(datatype_id))) H5Tclose(datatype_id);
      if(my_dims != NULL) free(my_dims);
      if(dims != NULL) free(dims);
      return NULL;

}

static PyObject *
Py_ReadListOfDatasets(PyObject *obj, PyObject *args)
{
    char *filename, *nodename;

    hid_t file_id;
    herr_t my_error;
    htri_t file_exists;
    H5T_class_t class_id;
    H5E_auto_t err_func;
    file_id = 0;

    if (!PyArg_ParseTuple(args, "ss",
            &filename, &nodename))
        return PyErr_Format(_hdf5ReadError,
               "ReadListOfDatasets: Invalid parameters.");

    /* How portable is this? */
    if (access(filename, R_OK) < 0) {
        PyErr_Format(_hdf5ReadError,
                 "ReadListOfDatasets: %s does not exist, or no read permissions\n",
                     filename);
        goto _fail;
    }

    file_exists = H5Fis_hdf5(filename);
    if (file_exists == 0) {
        PyErr_Format(_hdf5ReadError,
                 "ReadListOfDatasets: %s is not an HDF5 file", filename);
        goto _fail;
    }

    file_id = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT); 
    PyObject *nodelist = PyList_New(0);
    if (nodelist == NULL) {
        PyErr_Format(_hdf5ReadError,
                 "ReadListOfDatasets: List couldn't be made!");
        goto _fail;
    }

    my_error = H5Giterate(file_id, nodename, NULL, iterate_dataset, (void *) nodelist);
    H5Fclose(file_id);
    if (my_error) {
        PyErr_Format(_hdf5ReadError,
                 "ReadListOfDatasets: Problem iterating over HDF5 set.");
        goto _fail;
    }
    
    PyObject *return_value = Py_BuildValue("N", nodelist);
    return return_value;

    _fail:
      Py_XDECREF(nodelist);
      if(!(file_id <= 0)&&(H5Iget_ref(file_id))) H5Fclose(file_id);
      return NULL;
    
}

herr_t iterate_dataset(hid_t loc_id, const char *name, void *nodelist)
{
    H5G_stat_t statbuf;
    PyObject* node_name, node_list;

    H5Gget_objinfo(loc_id, name, 0, &statbuf);
    if (statbuf.type == H5G_DATASET) {
        node_name = PyString_FromString(name);
        if (node_name == NULL) {return -1;}
        if (PyList_Append((PyObject *)nodelist, node_name)) {return -1;}
    }
    return 0;
};


static PyMethodDef _hdf5LightReaderMethods[] = {
    {"ReadData", Py_ReadHDF5DataSet, METH_VARARGS},
    {"ReadDataSlice", Py_ReadHDF5DataSetSlice, METH_VARARGS},
    {"ReadListOfDatasets", Py_ReadListOfDatasets, METH_VARARGS},
    {NULL, NULL} 
};

#ifdef MS_WIN32
__declspec(dllexport)
#endif

void initHDF5LightReader(void)
{
    PyObject *m, *d;
    m = Py_InitModule("HDF5LightReader", _hdf5LightReaderMethods);
    d = PyModule_GetDict(m);
    _hdf5ReadError = PyErr_NewException("HDF5LightReader.ReadingError", NULL, NULL);
    PyDict_SetItemString(d, "ReadingError", _hdf5ReadError);
    import_array();
}
