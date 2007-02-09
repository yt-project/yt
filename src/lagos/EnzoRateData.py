#
# chemistry:
#   Some functions to include chemistry stuff,
#   although it should work for generalized
#   tables, if fed the appropriate key.
#
# Written by: Matthew Turk (mturk@stanford.edu) Nov 2006
# Modified:
#

from numarray import *
import numarray.nd_image as nd
import RavenCombine
import types, exceptions

class EnzoTable:
    # This is a class for storing tables of data from enzo.
    # Specifically, we will be storing chemical rates and cooling rates
    def __init__(self, filename, key):
        # We get fed the filename and a table that is the key to the column names
        self.cols = {}
        self.invcols = {}
        for i in range(len(key)):
            self.cols[key[i]] = i
            self.invcols[i] = key[i]
        self.filename = filename
        self.readTable()
    def clearAll(self):
        try:
            del self.columns
        except:
            pass
        try:
            del self.params
        except:
            pass
    def readTable(self):
        self.clearAll()
        lines = open(self.filename).readlines()
        # We may have comment lines, which we need to strip out
        # There may be a more elegant way to do this
        # We assume that all unit lines have only one value
        # This is slow, due to the stripping and splitting, but it shouldn't be
        # a huge portion of the whole runtime
        i = 0
        self.params = {}
        toArray = []
        for line in lines:
            if line.strip()[0] == "#":
                # It's a comment/unit line
                p, v = line[1:].split("=")
                self.params[p.replace(" ","")[1:]] = float(v)
            else:
                toArray.append(map(float, line.split()))
        self.columns = array(toArray, Float32)
        print "Found %s bins of values" % (self.columns.shape[0])
    def __getitem__(self, item):
        ## This WILL get cleaned up, but it does what I want for now
        x_vals = []
        toReshape = None
        if isinstance(item, types.TupleType):
            if isinstance(item[0], types.FloatType) or \
               isinstance(item[0], types.IntType):
                x_vals.append(float(item[0]))
            elif isinstance(item[0], ArrayType):
                toReshape = item[0].shape
                x_vals = (item[0].flat)
            colsToReturn = []
            for col in item[1:]:
                if isinstance(col, types.StringType):
                    colsToReturn.append(self.cols[col])
                else:
                    colsToReturn.append(int(col))
        elif isinstance(item, types.FloatType) or \
             isinstance(item, types.IntType):
            x_vals.append(float(item))
            colsToReturn=arange(1,len(self.cols))
        elif isinstance(item, ArrayType):
            toReshape = item.shape
            x_vals = item.flat
            colsToReturn=arange(1,len(self.cols))
        elif isinstance(item, types.ListType):
            colsToReturn=arange(1,len(self.cols))
        elif isinstance(item, types.StringType):
            return self.columns[:,self.cols[item]]
        else:
            raise exceptions.TypeError()
        colsToReturn = array(colsToReturn,Int32)
        valsToReturn = zeros((len(x_vals),len(colsToReturn)),Float32)
        ind = zeros(len(x_vals),Int64)
        x_axis = self.columns[:,0]
        x_vals_arr = array(x_vals, Float32)
        RavenCombine.Interpolate(x_axis, self.columns, x_vals_arr, valsToReturn, colsToReturn)
        if toReshape != None:
            if len(colsToReturn == 1):
                valsToReturn = reshape(valsToReturn, toReshape)
            else:
                newShape = list(toReshape)
                newShape.append(len(colsToReturn))
                valsToReturn = reshape(valsToReturn, newShape)
        return valsToReturn
