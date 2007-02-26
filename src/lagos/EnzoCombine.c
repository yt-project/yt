//
// point_combine:
//   A module for merging points from different grids
//
// Written by: Matthew Turk (mturk@stanford.edu) Nov 2006
// Modified:
//

#include "Python.h"

#include <stdio.h>
#include <math.h>
#include <signal.h>
#include <ctype.h>

#include "libnumarray.h"


// Sometimes a "maximum intensity" line-integral looks better
// switch these two defs, and then fix EnzoGrid.getProjection, to switch

//#define COMB(A,B) ((A) > (B) ? (A) : (B))
#define COMB(A,B) (A + B)

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))

static PyObject *_combineError;

static void 
RefineCoarseData(long fpoints, Int64 *finedata_x, Int64 *finedata_y, Float64 *finedata_vals, Float64 *finedata_wgt,
                 long cpoints, Int64 *coarsedata_x, Int64 *coarsedata_y, Float64 *coarsedata_vals, Float64 *coarsedata_wgt,
                 int refinementFactor, int *totalRefined)
{
    long fi, ci;

    int flagged[cpoints], tr = 0;

    long double rf = refinementFactor;

    Int64 coarseCell_x, coarseCell_y;

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
    finedata_x    = NA_IoArray(ofinedata_x   , tInt64, NUM_C_ARRAY | NUM_WRITABLE);
    finedata_y    = NA_IoArray(ofinedata_y   , tInt64, NUM_C_ARRAY | NUM_WRITABLE);
    finedata_vals = NA_IoArray(ofinedata_vals, tFloat64, NUM_C_ARRAY | NUM_WRITABLE);
    finedata_wgt  = NA_IoArray(ofinedata_wgt , tFloat64, NUM_C_ARRAY | NUM_WRITABLE);

    coarsedata_x    = NA_IoArray(ocoarsedata_x   , tInt64, NUM_C_ARRAY | NUM_WRITABLE);
    coarsedata_y    = NA_IoArray(ocoarsedata_y   , tInt64, NUM_C_ARRAY | NUM_WRITABLE);
    coarsedata_vals = NA_IoArray(ocoarsedata_vals, tFloat64, NUM_C_ARRAY | NUM_WRITABLE);
    coarsedata_wgt  = NA_IoArray(ocoarsedata_wgt , tFloat64, NUM_C_ARRAY | NUM_WRITABLE);

    //convolved = NA_OptionalOutputArray(oconvolved, tFloat64, NUM_C_ARRAY, data);

    if (!finedata_x || !finedata_y || !finedata_vals || !finedata_wgt
     || !coarsedata_x || !coarsedata_y || !coarsedata_vals || !coarsedata_wgt) {
        PyErr_Format( _combineError, 
                  "CombineData: error converting array inputs.");
        goto _fail;
    }

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

    /*if ((kernel->nd != 1) || (data->nd != 1)) {
        PyErr_Format(_combineError,
                 "CombineData: numarray must have 1 dimensions.");
        goto _fail;
    }*/

    int totalRefined;

    RefineCoarseData(finedata_vals->dimensions[0], 
                     NA_OFFSETDATA(finedata_x),
                     NA_OFFSETDATA(finedata_y),
                     NA_OFFSETDATA(finedata_vals),
                     NA_OFFSETDATA(finedata_wgt),
                     coarsedata_vals->dimensions[0],   
                     NA_OFFSETDATA(coarsedata_x),
                     NA_OFFSETDATA(coarsedata_y),
                     NA_OFFSETDATA(coarsedata_vals),
                     NA_OFFSETDATA(coarsedata_wgt),
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
CombineData(long dpoints, Int64 *alldata_x, Int64 *alldata_y, Float64 *alldata_vals, Int64 *alldata_mask, Float64 *alldata_wgt,
	    long gpoints, Int64 *griddata_x, Int64 *griddata_y, Float64 *griddata_vals, Int64 *griddata_mask, Float64 *griddata_wgt,
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
    alldata_x    = NA_IoArray(oalldata_x   , tInt64, NUM_C_ARRAY | NUM_WRITABLE);
    alldata_y    = NA_IoArray(oalldata_y   , tInt64, NUM_C_ARRAY | NUM_WRITABLE);
    alldata_vals = NA_IoArray(oalldata_vals, tFloat64, NUM_C_ARRAY | NUM_WRITABLE);
    alldata_wgt  = NA_IoArray(oalldata_wgt , tFloat64, NUM_C_ARRAY | NUM_WRITABLE);
    alldata_mask = NA_IoArray(oalldata_mask, tInt64, NUM_C_ARRAY | NUM_WRITABLE);

    griddata_x    = NA_IoArray(ogriddata_x   , tInt64, NUM_C_ARRAY | NUM_WRITABLE);
    griddata_y    = NA_IoArray(ogriddata_y   , tInt64, NUM_C_ARRAY | NUM_WRITABLE);
    griddata_vals = NA_IoArray(ogriddata_vals, tFloat64, NUM_C_ARRAY | NUM_WRITABLE);
    griddata_wgt  = NA_IoArray(ogriddata_wgt , tFloat64, NUM_C_ARRAY | NUM_WRITABLE);
    griddata_mask = NA_IoArray(ogriddata_mask, tInt64, NUM_C_ARRAY | NUM_WRITABLE);

    //convolved = NA_OptionalOutputArray(oconvolved, tFloat64, NUM_C_ARRAY, data);

    if (!alldata_x || !alldata_y || !alldata_vals || !alldata_mask || !alldata_wgt
     || !griddata_x || !griddata_y || !griddata_vals || !griddata_mask || !griddata_wgt) {
        PyErr_Format( _combineError, 
                  "CombineData: error converting array inputs.");
        goto _fail;
    }

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
    }

    /*if ((kernel->nd != 1) || (data->nd != 1)) {
        PyErr_Format(_combineError,
                 "CombineData: numarray must have 1 dimensions.");
        goto _fail;
    }*/

    CombineData(alldata_vals->dimensions[0], 
                NA_OFFSETDATA(alldata_x),
                NA_OFFSETDATA(alldata_y),
                NA_OFFSETDATA(alldata_vals),
                NA_OFFSETDATA(alldata_mask),
                NA_OFFSETDATA(alldata_wgt),		
                griddata_vals->dimensions[0],   
                NA_OFFSETDATA(griddata_x),
                NA_OFFSETDATA(griddata_y),
                NA_OFFSETDATA(griddata_vals),
                NA_OFFSETDATA(griddata_mask),
                NA_OFFSETDATA(griddata_wgt),
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
FindUpper(long ipoints, Float64 *input_axis, long vpoints, Float64 *wanted_vals, 
          Int64 *upper_inds)
{
    Int64 i;
    Int64 v;
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
    input_axis    =  NA_InputArray(oinput_axis   , tFloat64, NUM_C_ARRAY );
    desired_vals  =  NA_InputArray(odesired_vals , tFloat64, NUM_C_ARRAY );
    output_ind    = NA_OutputArray(ooutput_ind   ,   tInt64, NUM_C_ARRAY | NUM_WRITABLE);
    
    if (!input_axis || !desired_vals) {
        PyErr_Format( _combineError, 
                  "FindUpper: error converting array inputs.");
        goto _fail;
    }

    if (!NA_ShapeEqual(output_ind, desired_vals)) {
        PyErr_Format(_combineError,
                 "FindUpper: vals and output need identical shapes.");
        goto _fail;
    }
    FindUpper(input_axis->dimensions[0], 
              NA_OFFSETDATA(input_axis),
              desired_vals->dimensions[0],   
              NA_OFFSETDATA(desired_vals),
              NA_OFFSETDATA(output_ind));

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
Interpolate(long num_axis_points, Float64 *axis, PyArrayObject* table,
            PyArrayObject *desiredvals, long num_columns, Int32 *columns,
            PyArrayObject *outputvals)
{
    //int table_rows = table->dimensions[0];
    int num_desireds = desiredvals->dimensions[0];

    Int64 axis_ind;
    Int32 column;
    Int64 desired_num;
    
    Float64 desired;

    Float64 logtem0 = log10(axis[0]);
    Float64 logtem9 = log10(axis[num_axis_points-1]);
    Float64 dlogtem = (logtem9-logtem0)/(num_axis_points-1);
    Float64 t1, t2, tdef, ki, kip;

    for (desired_num = 0 ; desired_num < num_desireds ; desired_num++) {
        desired = log10l(NA_GET1(desiredvals, Float64, desired_num));
        axis_ind = min(num_axis_points-1,
                   max(0,(int)((desired-logtem0)/dlogtem)+1));
        t1 = (logtem0 + (axis_ind-1)*dlogtem);
        t2 = (logtem0 + (axis_ind+0)*dlogtem);
        tdef = t2 - t1;
        for (column = 0 ; column < num_columns ; column++) {
            ki  = NA_GET2(table, Float64, axis_ind-1, columns[column]);
            kip = NA_GET2(table, Float64, axis_ind+0, columns[column]);
            NA_SET2(outputvals, Float64, desired_num, column, 
                    ki+(desired-t1)*(kip-ki)/tdef );
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
    axis          =  NA_InputArray(oaxis         , tFloat64, NUM_C_ARRAY );
    table         =  NA_InputArray(otable        , tFloat64, 0           );
    desired       =  NA_InputArray(odesired      , tFloat64, 0           );
    outputvals    = NA_OutputArray(ooutputvals   , tFloat64, NUM_WRITABLE);
    columns       =  NA_InputArray(ocolumns      ,   tInt32, NUM_C_ARRAY );

    if (!axis || !table || !desired || !outputvals || !columns) {
        PyErr_Format( _combineError, 
                  "Interpolate: error converting array inputs.");
        goto _fail;
    }

    if (columns->dimensions[0] != outputvals->dimensions[1]) {
        PyErr_Format(_combineError,
                 "Interpolate: number of columns requested must match number "
                 "of columns in output buffer. %i", columns->dimensions[0]);
        goto _fail;
    }

    Interpolate(axis->dimensions[0], 
              NA_OFFSETDATA(axis),
              table, desired,
              columns->dimensions[0],
              NA_OFFSETDATA(columns),
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

static PyMethodDef _combineMethods[] = {
    {"RefineCoarseData", Py_RefineCoarseData, METH_VARARGS},
    {"CombineData", Py_CombineData, METH_VARARGS},
    {"FindUpper", Py_FindUpper, METH_VARARGS},
    {"Interpolate", Py_Interpolate, METH_VARARGS},
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
    import_libnumarray();
}

/*
 * Local Variables:
 * mode: C
 * c-file-style: "python"
 * End:
 */
