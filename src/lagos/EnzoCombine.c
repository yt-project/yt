//
// EnzoCombine
//   A module for merging points from different grids, in various ways.
//   Used for projections, interpolations, and binning profiles.
//
// Written by: Matthew Turk (mturk@stanford.edu) Nov 2006
// Modified:
//

#include "Python.h"

#include <stdio.h>
#include <math.h>
#include <signal.h>
#include <ctype.h>

#include "numpy/ndarrayobject.h"


// Sometimes a "maximum intensity" line-integral looks better
// switch these two defs, and then fix EnzoGrid.getProjection, to switch

//#define COMB(A,B) ((A) > (B) ? (A) : (B))
#define COMB(A,B) (A + B)

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))

static PyObject *_combineError;

static void 
RefineCoarseData(long fpoints, npy_int64 *finedata_x, npy_int64 *finedata_y, npy_float64 *finedata_vals, npy_float64 *finedata_wgt,
                 long cpoints, npy_int64 *coarsedata_x, npy_int64 *coarsedata_y, npy_float64 *coarsedata_vals, npy_float64 *coarsedata_wgt,
                 int refinementFactor, int *totalRefined)
{
    long fi, ci;

    int flagged[cpoints], tr = 0;

    long double rf = refinementFactor;

    npy_int64 coarseCell_x, coarseCell_y;

    for (ci=0; ci<cpoints; ci++) flagged[ci]=0;

    for (fi = 0; fi < fpoints; fi++) {
        coarseCell_x = floorl(finedata_x[fi]/rf);
        coarseCell_y = floorl(finedata_y[fi]/rf);
        for (ci = 0; ci < cpoints; ci++) {
            if ((coarseCell_x == coarsedata_x[ci]) &&
                (coarseCell_y == coarsedata_y[ci])) {
                    tr += 1;
                    finedata_vals[fi] = COMB(coarsedata_vals[ci], finedata_vals[fi]);
                    finedata_wgt[fi]  = COMB(coarsedata_wgt[ci],  finedata_wgt[fi]);
                    flagged[ci] = 1;
                    break;  // Each fine cell maps to one and only one coarse
                            // cell
            }
        }
    }
    *totalRefined = tr;

    for (ci=0; ci<cpoints; ci++) {
        if (flagged[ci]==1) {
            coarsedata_x[ci]=-1;
            coarsedata_y[ci]=-1;
            coarsedata_vals[ci]=0.0;
            coarsedata_wgt[ci]=0.0;
        }
    }
}

static PyObject *
Py_RefineCoarseData(PyObject *obj, PyObject *args)
{
    PyObject   *ofinedata_x, *ofinedata_y, *ofinedata_vals, *ofinedata_wgt,
               *ocoarsedata_x, *ocoarsedata_y, *ocoarsedata_vals, *ocoarsedata_wgt,
               *oreturn_index;
    PyArrayObject   *finedata_x, *finedata_y, *finedata_vals, *finedata_wgt,
	    *coarsedata_x, *coarsedata_y, *coarsedata_vals, *coarsedata_wgt;
    int refinementFactor;

    if (!PyArg_ParseTuple(args, "OOOOOOOOi", 
        &ofinedata_x, &ofinedata_y, &ofinedata_vals, &ofinedata_wgt,
        &ocoarsedata_x, &ocoarsedata_y, &ocoarsedata_vals, &ocoarsedata_wgt,
        &refinementFactor))
        return PyErr_Format(_combineError, 
                    "CombineData: Invalid parameters.");

    /* Align, Byteswap, Contiguous, Typeconvert */

    finedata_x    = (PyArrayObject *) PyArray_FromAny(ofinedata_x   , PyArray_DescrFromType(NPY_INT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL);
    finedata_y    = (PyArrayObject *) PyArray_FromAny(ofinedata_y   , PyArray_DescrFromType(NPY_INT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL);
    finedata_vals = (PyArrayObject *) PyArray_FromAny(ofinedata_vals, PyArray_DescrFromType(NPY_FLOAT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL);
    finedata_wgt  = (PyArrayObject *) PyArray_FromAny(ofinedata_wgt , PyArray_DescrFromType(NPY_FLOAT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL);

    coarsedata_x    = (PyArrayObject *) PyArray_FromAny(ocoarsedata_x   , PyArray_DescrFromType(NPY_INT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL);
    coarsedata_y    = (PyArrayObject *) PyArray_FromAny(ocoarsedata_y   , PyArray_DescrFromType(NPY_INT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL);
    coarsedata_vals = (PyArrayObject *) PyArray_FromAny(ocoarsedata_vals, PyArray_DescrFromType(NPY_FLOAT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL);
    coarsedata_wgt  = (PyArrayObject *) PyArray_FromAny(ocoarsedata_wgt , PyArray_DescrFromType(NPY_FLOAT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL);

    if (!finedata_x || !finedata_y || !finedata_vals || !finedata_wgt
     || !coarsedata_x || !coarsedata_y || !coarsedata_vals || !coarsedata_wgt) {
        /*PyErr_Format( _combineError, 
                  "CombineData: error converting array inputs.");
        */
        if (!coarsedata_x)
          PyErr_Format( _combineError, 
              "CombineData: error converting coarsedata_x");
        if (!coarsedata_y)
          PyErr_Format( _combineError, 
              "CombineData: error converting coarsedata_y");
        if (!coarsedata_vals)
          PyErr_Format( _combineError, 
              "CombineData: error converting coarsedata_vals");
        if (!coarsedata_wgt)
          PyErr_Format( _combineError, 
              "CombineData: error converting coarsedata_wgt");
        if (!finedata_x)
          PyErr_Format( _combineError, 
              "CombineData: error converting finedata_x");
        if (!finedata_y)
          PyErr_Format( _combineError, 
              "CombineData: error converting finedata_y");
        if (!finedata_vals)
          PyErr_Format( _combineError, 
              "CombineData: error converting finedata_vals");
        if (!finedata_wgt)
          PyErr_Format( _combineError, 
              "CombineData: error converting finedata_wgt");
        goto _fail;
    }

/*
    if (!NA_ShapeEqual(finedata_x, finedata_y)) {
        PyErr_Format(_combineError,
                 "CombineData: x and y finedata numarray need identitcal shapes.");
        goto _fail;
    }

    if (!NA_ShapeEqual(finedata_x, finedata_vals)) {
        PyErr_Format(_combineError,
                 "CombineData: x and vals finedata numarray need identitcal shapes.");
        goto _fail;
    }

    if (!NA_ShapeEqual(finedata_x, finedata_wgt)) {
        PyErr_Format(_combineError,
                 "CombineData: x and weight finedata numarray need identitcal shapes.");
        goto _fail;
    }

    if (!NA_ShapeEqual(coarsedata_x, coarsedata_y)) {
        PyErr_Format(_combineError,
                 "CombineData: x and y coarsedata numarray need identitcal shapes.");
        goto _fail;
    }

    if (!NA_ShapeEqual(coarsedata_x, coarsedata_vals)) {
        PyErr_Format(_combineError,
                 "CombineData: x and vals coarsedata numarray need identitcal shapes.");
        goto _fail;
    }

    if (!NA_ShapeEqual(coarsedata_x, coarsedata_wgt)) {
        PyErr_Format(_combineError,
                 "CombineData: x and weight coarsedata numarray need identitcal shapes.");
        goto _fail;
    }
*/

    int totalRefined;

    RefineCoarseData(finedata_vals->dimensions[0], 
                     (npy_int64 *) finedata_x->data,
                     (npy_int64 *) finedata_y->data,
                     (npy_float64 *) finedata_vals->data,
                     (npy_float64 *) finedata_wgt->data,
                     coarsedata_vals->dimensions[0],   
                     (npy_int64 *) coarsedata_x->data,
                     (npy_int64 *) coarsedata_y->data,
                     (npy_float64 *) coarsedata_vals->data,
                     (npy_float64 *) coarsedata_wgt->data,
                     refinementFactor, &totalRefined);

    Py_XDECREF(finedata_x);
    Py_XDECREF(finedata_y);
    Py_XDECREF(finedata_vals);
    Py_XDECREF(finedata_wgt);
    Py_XDECREF(coarsedata_x);
    Py_XDECREF(coarsedata_y);
    Py_XDECREF(coarsedata_vals);
    Py_XDECREF(coarsedata_wgt);

    /* Align, Byteswap, Contiguous, Typeconvert */
    oreturn_index = PyInt_FromLong((long)totalRefined);
    return oreturn_index;

  _fail:
    Py_XDECREF(finedata_x);
    Py_XDECREF(finedata_y);
    Py_XDECREF(finedata_vals);
    Py_XDECREF(finedata_wgt);
    Py_XDECREF(coarsedata_x);
    Py_XDECREF(coarsedata_y);
    Py_XDECREF(coarsedata_vals);
    Py_XDECREF(coarsedata_wgt);

    return NULL;      
}

static void 
CombineData(long dpoints, npy_int64 *alldata_x, npy_int64 *alldata_y, npy_float64 *alldata_vals, npy_int64 *alldata_mask, npy_float64 *alldata_wgt,
	    long gpoints, npy_int64 *griddata_x, npy_int64 *griddata_y, npy_float64 *griddata_vals, npy_int64 *griddata_mask, npy_float64 *griddata_wgt,
                 int *lastIndex)
{
    long di, gi;

    int appendData;
    int myIndex = *lastIndex;
    myIndex = 0;

    //for (gi=0; gi<gpoints; gi++) flagged[gi]=0;

    for (di = 0; di < dpoints; di++) {
        appendData = 0;
        if (alldata_x[di] < 0) continue;
        for (gi = 0; gi < gpoints; gi++) {
            if ((alldata_x[di] == griddata_x[gi]) &&
                (alldata_y[di] == griddata_y[gi])) {
                    alldata_vals[di] = COMB(alldata_vals[di], griddata_vals[gi]);
                    alldata_wgt[di]  = COMB(alldata_wgt[di], griddata_wgt[gi]);
                    alldata_mask[di] = (alldata_mask[di] && griddata_mask[gi]);
                    griddata_x[gi] = -1;
                    appendData = 0;
                    myIndex++;
                    break; // We map to one alldata point AT MOST
            }
        }
        if (appendData) {
            alldata_x[myIndex] = griddata_x[gi];
            alldata_y[myIndex] = griddata_y[gi];
            alldata_mask[myIndex] = griddata_mask[gi];
            alldata_vals[myIndex] = griddata_vals[gi];
            alldata_wgt[myIndex] = griddata_wgt[gi];
            myIndex++;
        }
    }
    *lastIndex = myIndex;
}

static PyObject *
Py_CombineData(PyObject *obj, PyObject *args)
{
    PyObject   *oalldata_x, *oalldata_y, *oalldata_vals, *oalldata_mask, *oalldata_wgt,
    	       *ogriddata_x, *ogriddata_y, *ogriddata_vals, *ogriddata_mask, *ogriddata_wgt,
               *oreturn_index;
    PyArrayObject   *alldata_x, *alldata_y, *alldata_vals, *alldata_mask, *alldata_wgt,
	    *griddata_x, *griddata_y, *griddata_vals, *griddata_mask, *griddata_wgt;
    int lastIndex;

    if (!PyArg_ParseTuple(args, "OOOOOOOOOOi", 
        &oalldata_x, &oalldata_y, &oalldata_vals, &oalldata_mask, &oalldata_wgt,
        &ogriddata_x, &ogriddata_y, &ogriddata_vals, &ogriddata_mask, &ogriddata_wgt,
        &lastIndex))
        return PyErr_Format(_combineError, 
                    "CombineData: Invalid parameters.");

    /* Align, Byteswap, Contiguous, Typeconvert */
    alldata_x    = (PyArrayObject *) PyArray_FromAny(oalldata_x   , PyArray_DescrFromType(NPY_INT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL);
    alldata_y    = (PyArrayObject *) PyArray_FromAny(oalldata_y   , PyArray_DescrFromType(NPY_INT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL);
    alldata_vals = (PyArrayObject *) PyArray_FromAny(oalldata_vals, PyArray_DescrFromType(NPY_FLOAT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL);
    alldata_wgt  = (PyArrayObject *) PyArray_FromAny(oalldata_wgt , PyArray_DescrFromType(NPY_FLOAT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL);
    alldata_mask = (PyArrayObject *) PyArray_FromAny(oalldata_mask, PyArray_DescrFromType(NPY_INT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL);

    griddata_x    = (PyArrayObject *) PyArray_FromAny(ogriddata_x   , PyArray_DescrFromType(NPY_INT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL);
    griddata_y    = (PyArrayObject *) PyArray_FromAny(ogriddata_y   , PyArray_DescrFromType(NPY_INT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL);
    griddata_vals = (PyArrayObject *) PyArray_FromAny(ogriddata_vals, PyArray_DescrFromType(NPY_FLOAT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL);
    griddata_wgt  = (PyArrayObject *) PyArray_FromAny(ogriddata_wgt , PyArray_DescrFromType(NPY_FLOAT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL);
    griddata_mask = (PyArrayObject *) PyArray_FromAny(ogriddata_mask, PyArray_DescrFromType(NPY_INT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL);

    if (!alldata_x || !alldata_y || !alldata_vals || !alldata_wgt || !alldata_mask
     || !griddata_x || !griddata_y || !griddata_vals || !griddata_wgt || !griddata_mask) {
        if (!alldata_x)
          PyErr_Format( _combineError, 
              "CombineData: error converting alldata_x");
        if (!alldata_y)
          PyErr_Format( _combineError, 
              "CombineData: error converting alldata_y");
        if (!alldata_vals)
          PyErr_Format( _combineError, 
              "CombineData: error converting alldata_vals");
        if (!alldata_wgt)
          PyErr_Format( _combineError, 
              "CombineData: error converting alldata_wgt");
        if (!alldata_mask)
          PyErr_Format( _combineError, 
              "CombineData: error converting alldata_mask");
        if (!griddata_x)
          PyErr_Format( _combineError, 
              "CombineData: error converting griddata_x");
        if (!griddata_y)
          PyErr_Format( _combineError, 
              "CombineData: error converting griddata_y");
        if (!griddata_vals)
          PyErr_Format( _combineError, 
              "CombineData: error converting griddata_vals");
        if (!griddata_wgt)
          PyErr_Format( _combineError, 
              "CombineData: error converting griddata_wgt");
        if (!griddata_mask)
          PyErr_Format( _combineError, 
              "CombineData: error converting griddata_mask");
        goto _fail;
    }


/*
    if (!NA_ShapeEqual(alldata_x, alldata_y)) {
        PyErr_Format(_combineError,
                 "CombineData: x and y alldata numarray need identitcal shapes.");
        goto _fail;
    }

    if (!NA_ShapeEqual(alldata_x, alldata_vals)) {
        PyErr_Format(_combineError,
                 "CombineData: x and vals alldata numarray need identitcal shapes.");
        goto _fail;
    }

    if (!NA_ShapeEqual(alldata_x, alldata_wgt)) {
        PyErr_Format(_combineError,
                 "CombineData: x and weight alldata numarray need identitcal shapes.");
        goto _fail;
    }

    if (!NA_ShapeEqual(griddata_x, griddata_y)) {
        PyErr_Format(_combineError,
                 "CombineData: x and y griddata numarray need identitcal shapes.");
        goto _fail;
    }

    if (!NA_ShapeEqual(griddata_x, griddata_vals)) {
        PyErr_Format(_combineError,
                 "CombineData: x and vals griddata numarray need identitcal shapes.");
        goto _fail;
    }

    if (!NA_ShapeEqual(griddata_wgt, griddata_vals)) {
        PyErr_Format(_combineError,
                 "CombineData: weight and vals griddata numarray need identitcal shapes.");
        goto _fail;
    }*/

    CombineData(alldata_vals->dimensions[0], 
                (npy_int64 *) alldata_x->data,
                (npy_int64 *) alldata_y->data,
                (npy_float64 *) alldata_vals->data,
                (npy_int64 *) alldata_mask->data,
                (npy_float64 *) alldata_wgt->data,		
                griddata_vals->dimensions[0],   
                (npy_int64 *) griddata_x->data,
                (npy_int64 *) griddata_y->data,
                (npy_float64 *) griddata_vals->data,
                (npy_int64 *) griddata_mask->data,
                (npy_float64 *) griddata_wgt->data,
                &lastIndex);

    Py_XDECREF(alldata_x);
    Py_XDECREF(alldata_y);
    Py_XDECREF(alldata_vals);
    Py_XDECREF(alldata_mask);
    Py_XDECREF(alldata_wgt);
    Py_XDECREF(griddata_x);
    Py_XDECREF(griddata_y);
    Py_XDECREF(griddata_vals);
    Py_XDECREF(griddata_mask);
    Py_XDECREF(griddata_wgt);

    /* Align, Byteswap, Contiguous, Typeconvert */
    oreturn_index = PyInt_FromLong((long)lastIndex);
    return oreturn_index;

  _fail:
    Py_XDECREF(alldata_x);
    Py_XDECREF(alldata_y);
    Py_XDECREF(alldata_vals);
    Py_XDECREF(alldata_wgt);
    Py_XDECREF(griddata_x);
    Py_XDECREF(griddata_y);
    Py_XDECREF(griddata_vals);
    Py_XDECREF(griddata_wgt);

    return NULL;      
}

static void 
FindUpper(long ipoints, npy_float64 *input_axis, long vpoints, npy_float64 *wanted_vals, 
          npy_int64 *upper_inds)
{
    npy_int64 i;
    npy_int64 v;
    for (v = 0 ; v < vpoints ; v++)
        for (i = 0 ; i < ipoints ; i++)
            if (input_axis[i] >= wanted_vals[v]) {
                upper_inds[v] = i;
                break;
            }
    return;
}

static PyObject *
Py_FindUpper(PyObject *obj, PyObject *args)
{
    PyObject   *oinput_axis, *odesired_vals, *ooutput_ind;
    PyArrayObject   *input_axis, *desired_vals, *output_ind;

    if (!PyArg_ParseTuple(args, "OOO", 
        &oinput_axis, &odesired_vals, &ooutput_ind))
        return PyErr_Format(_combineError, 
                    "FindUpper: Invalid parameters.");

    /* Align, Byteswap, Contiguous, Typeconvert */
    input_axis    =  (PyArrayObject *) PyArray_FromAny(oinput_axis   , PyArray_DescrFromType(NPY_FLOAT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL );
    desired_vals  =  (PyArrayObject *) PyArray_FromAny(odesired_vals , PyArray_DescrFromType(NPY_FLOAT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL );
    output_ind    =  (PyArrayObject *) PyArray_FromAny(ooutput_ind   ,   PyArray_DescrFromType(NPY_INT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL);
    
    if (!input_axis || !desired_vals) {
        PyErr_Format( _combineError, 
                  "FindUpper: error converting array inputs.");
        goto _fail;
    }

/*
    if (!NA_ShapeEqual(output_ind, desired_vals)) {
        PyErr_Format(_combineError,
                 "FindUpper: vals and output need identical shapes.");
        goto _fail;
    }
*/
    FindUpper(input_axis->dimensions[0], 
              (npy_float64 *) input_axis->data,
              desired_vals->dimensions[0],   
              (npy_float64 *) desired_vals->data,
              (npy_int64 *) output_ind->data);

    Py_XDECREF(input_axis);
    Py_XDECREF(desired_vals);
    Py_XDECREF(output_ind);

    /* Align, Byteswap, Contiguous, Typeconvert */
    return Py_None;

  _fail:
    Py_XDECREF(input_axis);
    Py_XDECREF(desired_vals);
    Py_XDECREF(output_ind);

    return NULL;      
}

static void 
Interpolate(long num_axis_points, npy_float64 *axis, PyArrayObject* table,
            PyArrayObject *desiredvals, long num_columns, npy_int32 *columns,
            PyArrayObject *outputvals)
{
    //int table_rows = table->dimensions[0];
    int num_desireds = desiredvals->dimensions[0];

    npy_int axis_ind, col_ind;
    npy_int32 column;
    npy_int64 desired_num;
    
    npy_float64 desired;

    npy_float64 logtem0 = log10(axis[0]);
    npy_float64 logtem9 = log10(axis[num_axis_points-1]);
    npy_float64 dlogtem = (logtem9-logtem0)/(num_axis_points-1);
    npy_float64 t1, t2, tdef, ki, kip;
    npy_float64 *t;

    for (desired_num = 0 ; desired_num < num_desireds ; desired_num++) {
        t = (npy_float64*)PyArray_GETPTR1(desiredvals, desired_num);
        desired = log10l(*t);
        axis_ind = min(num_axis_points-1,
                   max(0,(int)((desired-logtem0)/dlogtem)+1));
        t1 = (logtem0 + (axis_ind-1)*dlogtem);
        t2 = (logtem0 + (axis_ind+0)*dlogtem);
        tdef = t2 - t1;
        for (column = 0 ; column < num_columns ; column++) {
            col_ind = (npy_int) columns[column];
            ki  = *(npy_float64*)PyArray_GETPTR2(table, (npy_int) (axis_ind-1), col_ind);
            kip = *(npy_float64*)PyArray_GETPTR2(table, (npy_int) (axis_ind+0), col_ind);
            *(npy_float64*) PyArray_GETPTR2(outputvals, desired_num, column) =
                    ki+(desired-t1)*(kip-ki)/tdef;
        }
    }
    return;
}

static PyObject *
Py_Interpolate(PyObject *obj, PyObject *args)
{
    PyObject   *oaxis, *otable, *odesired, *ooutputvals, *ocolumns;
    PyArrayObject   *axis, *table, *desired, *outputvals, *columns;

    if (!PyArg_ParseTuple(args, "OOOOO", 
        &oaxis, &otable, &odesired, &ooutputvals, &ocolumns))
        return PyErr_Format(_combineError, 
                    "Interpolate: Invalid parameters.");

    /* Align, Byteswap, Contiguous, Typeconvert */
    axis          =  (PyArrayObject *) PyArray_FromAny(oaxis         , PyArray_DescrFromType(NPY_FLOAT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL );
    table         =  (PyArrayObject *) PyArray_FromAny(otable        , PyArray_DescrFromType(NPY_FLOAT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL );
    desired       =  (PyArrayObject *) PyArray_FromAny(odesired      , PyArray_DescrFromType(NPY_FLOAT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL );
    outputvals    =  (PyArrayObject *) PyArray_FromAny(ooutputvals   , PyArray_DescrFromType(NPY_FLOAT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL );
    columns       =  (PyArrayObject *) PyArray_FromAny(ocolumns      ,   PyArray_DescrFromType(NPY_INT32), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL );

    if (!axis || !table || !desired || !outputvals || !columns) {
        PyErr_Format( _combineError, 
                  "Interpolate: error converting array inputs.");
        goto _fail;
    }

    if (columns->dimensions[0] != outputvals->dimensions[1]) {
        PyErr_Format(_combineError,
                 "Interpolate: number of columns requested must match number "
                 "of columns in output buffer. %i", (int) columns->dimensions[0]);
        goto _fail;
    }

    Interpolate(axis->dimensions[0], 
              (npy_float64 *) axis->data,
              table, desired,
              columns->dimensions[0],
              (npy_int32 *) columns->data,
              outputvals);
    Py_XDECREF(axis);
    Py_XDECREF(table);
    Py_XDECREF(desired);
    Py_XDECREF(outputvals);
    Py_XDECREF(columns);

    /* Align, Byteswap, Contiguous, Typeconvert */
    return Py_None;

  _fail:
    Py_XDECREF(axis);
    Py_XDECREF(table);
    Py_XDECREF(desired);
    Py_XDECREF(outputvals);
    Py_XDECREF(columns);


    return NULL;      
}

static void 
BinProfile(long num_bins, npy_int64 *binindices, npy_float64 *profilevalues,
           long num_elements, npy_float64 *fieldvalues, npy_float64 *weightvalues)
{
    // We are going to assume -- making an ass out of you and me -- that the
    // values of the profile are already zero.  

    npy_float64 weights[num_bins];
    npy_float64 myweight;
    npy_int64 mybin = -1;
    int bin;
    int element;

    for (bin = 0; bin < num_bins ; bin++) {
        weights[bin] = 0;
    }

    if (weightvalues[0] != -999) {
        for (element = 0; element < num_elements ; element++) {
            mybin = binindices[element];
            myweight = weightvalues[mybin];
            profilevalues[mybin] += myweight*fieldvalues[element];
            weights[mybin] += myweight;
        }
        for (bin = 0; bin < num_bins ; bin++) {
            profilevalues[bin] = profilevalues[bin] / weights[bin];
        }
    } else {
        for (element = 0; element < num_elements ; element++) {
            mybin = binindices[element];
            myweight = weightvalues[mybin];
            profilevalues[mybin] += fieldvalues[element];
        }
        for (bin = 1; bin < num_bins ; bin++) {
            profilevalues[bin] += profilevalues[bin-1];
            //printf("ProfileValues: %i\t%0.4e\n", bin, profilevalues[bin]);
        }
    }

    // Okay, and we're done.

}

static PyObject *
Py_BinProfile(PyObject *obj, PyObject *args)
{
    PyObject   *ofield, *obinindices, *oprofile, *oweightfield;
    PyArrayObject   *field, *binindices, *profile, *weightfield;

    if (!PyArg_ParseTuple(args, "OOOO",
                &ofield, &obinindices, &oprofile, &oweightfield))
        return PyErr_Format(_combineError, 
                    "Interpolate: Invalid parameters.");

    /* Align, Byteswap, Contiguous, Typeconvert */
    field       =  (PyArrayObject *) PyArray_FromAny(ofield         , PyArray_DescrFromType(NPY_FLOAT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL);
    binindices  =  (PyArrayObject *) PyArray_FromAny(obinindices    , PyArray_DescrFromType(NPY_INT64)  , 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL);
    profile     =  (PyArrayObject *) PyArray_FromAny(oprofile       , PyArray_DescrFromType(NPY_FLOAT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL);
    weightfield =  (PyArrayObject *) PyArray_FromAny(oweightfield   , PyArray_DescrFromType(NPY_FLOAT64), 1, 0, NPY_ENSURECOPY | NPY_UPDATEIFCOPY, NULL);

    if (!field){
        PyErr_Format( _combineError, 
                  "BinProfile: error converting array inputs. (field)");
        goto _fail;
    }

    if (!binindices){
        PyErr_Format( _combineError, 
                  "BinProfile: error converting array inputs. (binindices)");
        goto _fail;
    }

    if (!profile ){
        PyErr_Format( _combineError, 
                  "BinProfile: error converting array inputs. (profile)");
        goto _fail;
    }

    if (!weightfield) {
        PyErr_Format( _combineError, 
                  "BinProfile: error converting array inputs. (weightfield)");
        goto _fail;
    }

    if (field->dimensions[0] != binindices->dimensions[0]) {
        PyErr_Format(_combineError,
                 "BinProfile: number of field values must match number of "
                  "indices.  %i, %i", (int) field->dimensions[0],
                                      (int) binindices->dimensions[0] );
        goto _fail;
    }

    BinProfile(profile->dimensions[0],
               (npy_int64 *) binindices->data,
               (npy_float64 *) profile->data,
               field->dimensions[0],
               (npy_float64 *) field->data,
               (npy_float64 *) weightfield->data);

    Py_XDECREF(field);
    Py_XDECREF(binindices);
    Py_XDECREF(profile);
    Py_XDECREF(weightfield);

    return Py_None;

  _fail:
    Py_XDECREF(field);
    Py_XDECREF(binindices);
    Py_XDECREF(profile);
    Py_XDECREF(weightfield);

    return NULL;      
}

static PyMethodDef _combineMethods[] = {
    {"RefineCoarseData", Py_RefineCoarseData, METH_VARARGS},
    {"CombineData", Py_CombineData, METH_VARARGS},
    {"FindUpper", Py_FindUpper, METH_VARARGS},
    {"Interpolate", Py_Interpolate, METH_VARARGS},
    {"BinProfile", Py_BinProfile, METH_VARARGS},
    {NULL, NULL} /* Sentinel */
};

/* platform independent*/
#ifdef MS_WIN32
__declspec(dllexport)
#endif

void initEnzoCombine(void)
{
    PyObject *m, *d;
    m = Py_InitModule("EnzoCombine", _combineMethods);
    d = PyModule_GetDict(m);
    _combineError = PyErr_NewException("EnzoCombine.error", NULL, NULL);
    PyDict_SetItemString(d, "error", _combineError);
    import_array();
}

/*
 * Local Variables:
 * mode: C
 * c-file-style: "python"
 * End:
 */
