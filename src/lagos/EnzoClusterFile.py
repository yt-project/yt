"""
Some functions to include the output from enzo_anyl, and generate it if need
be.

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
@todo: Create an AnalyzeCluster file from a given set of parameters, and a
hierarchy file.
@todo: Run enzo_anyl and parse the resultant data.
@todo: Derived fields, somehow?  This could be done with an addColumn.  This
should definitely be separate from EnzoDerivedFields, I think.
@todo: Merge back into EnzoDataTypes.py -- separate for now while still
testing.  (Does it make sense to put it in there?)
"""

from yt.lagos import *
import types, exceptions

class AnalyzeClusterOutput:
    # This is a class for storing the results of enzo_anyl runs
    def __init__(self, filename, id=None):
        """
        We receive a filename, and then parse it to get the key of columns, and
        then parse the columns to get the array of values
        """
        self.filename = filename
        if isinstance(filename, types.StringType):  # For the __add__ method
            self.lines = open(self.filename).readlines()
            self.parseKey()
            self.parseData()
        elif isinstance(filename, types.DictType):
            self.filename = id
            ks = filename.keys()
            ks.sort()
            numBins = filename[ks[0]].shape[0]
            self.data = na.zeros((len(ks), numBins), nT.Float64)
            self.columns = {}
            self.rcolumns = {}
            i = 1
            for column in ks:
                #print column, self.data.shape, self.data[i,:].shape, filename[column].shape
                self.data[i-1,:] = filename[column]
                self.columns[i] = column
                self.rcolumns[column] = i
                i += 1

    def __getitem__(self, item):
        if isinstance(item, types.StringType):
            i = self.rcolumns[item]
        else:
            i = int(item)
        return self.data[i]

    def parseKey(self):
        """
        Scans, grabs the column names
        """
        self.columns = {}
        self.rcolumns = {}
        index = self.lines.index("# COLUMN   DESCRIPTION\n")
        for line in self.lines[index+1:]:
            if line.rstrip() == "#":
                break
            col, name = line[1:].strip().split(" ",1)
            col = int(col) - 1
            self.columns[col] = name.strip().rstrip()
            self.rcolumns[name.strip().rstrip()] = col

    def parseData(self):
        """
        Parse the data into self.data

        @todo: Convert to use a regular expression and counting the number of
        matches
        """
        numBins = 0
        for line in self.lines:
            if line[0] != "#" and len(line.strip()) > 0:
                numBins += 1
        self.data = na.zeros((len(self.columns), numBins), nT.Float64)
        i = 0
        nc = len(self.columns)
        for line in self.lines:
            if line[0] != "#" and len(line.strip()) > 0:
                self.data[:,i] = map(float, line.split()[:nc])
                i += 1

    def outputHippo(self, filename):
        """
        This outputs the columns in a format that HippoDraw can read easily
        """
        fs = "\t".join(["%0.15e"] * len(self.columns)) + "\n"
        f = open(filename,"w")
        # Dicts are unsorted, so we do this.  Could probably be better.
        fields = [] 
        keys = self.columns.keys()
        keys.sort()
        for key in keys:
            fields.append(self.columns[key])
        header = "\t".join(fields) + "\n"
        f.write(header)
        for i in xrange(self.data.shape[1]):
            f.write(fs % tuple(self.data[:,i]))
        f.close()

    def __add__(self, other):
        """
        This will add all non-duplicated columns.

        @note: We kind of have to do this the slow way, to ensure that no
        double-columning goes on.  I couldn't come up with a better way than
        this.  Maybe I will, some day.  But for now, it works!  Hooray!
        """
        comb = AnalyzeClusterOutput(None)
        comb.columns = self.columns.copy()
        comb.rcolumns = self.rcolumns.copy()
        cols = other.columns.keys()
        offset = len(comb.columns) # 0-indexed
        i = offset
        for col in cols:
            if other.columns[col] not in self.columns.values():
                comb.columns[i] = other.columns[col]
                comb.rcolumns[col] = i
                i += 1
        comb.data = na.zeros((len(comb.columns), self.data.shape[1]), nT.Float64)
        comb.data[:len(self.columns),:] = self.data
        i = offset
        for col in cols:
            if other.columns[col] not in self.columns.values():
                comb.data[i,:] = other.data[col,:] # data is zero-indexed
                i += 1
        #print comb.data[offset:,:].shape, other.data[1:,:].shape
        #comb.data[offset:,:] = other.data[1:,:]
        return comb
