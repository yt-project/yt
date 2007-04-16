"""
Some functions to include the output from enzo_anyl, and generate it if need
be.

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
@todo: Create an AnalyzeCluster file from a given set of parameters, and a
hierarchy file.
@todo: Run enzo_anyl and parse the resultant data.
@todo: Implement the + operator so that we can add two different ClusterFiles
together (for instance, the Species and the plain) and get a resultant,
enormous one.
@todo: Derived fields, somehow?  This could be done with an addColumn.  This
should definitely be separate from EnzoDerivedFields, I think.
@todo: Merge back into EnzoDataTypes.py -- separate for now while still
testing.  (Does it make sense to put it in there?)
"""

from yt.lagos import *
#from numarray import *
#import numarray.nd_image as nd
import types, exceptions

class AnalyzeClusterOutput:
    # This is a class for storing the results of enzo_anyl runs
    def __init__(self, filename):
        """
        We receive a filename, and then parse it to get the key of columns, and
        then parse the columns to get the array of values
        """
        self.filename = filename
        self.lines = open(self.filename).readlines()
        self.parseKey()
        self.parseData()

    def parseKey(self):
        """
        Scans, grabs the column names
        """
        self.columns = {}
        index = self.lines.index("# COLUMN   DESCRIPTION\n")
        for line in self.lines[index+1:]:
            if line.rstrip() == "#":
                break
            col, name = line[1:].strip().split(" ",1)
            col = int(col)
            self.columns[col] = name.strip().rstrip()

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
        self.data = na.zeros((len(self.columns), numBins), na.Float64)
        i = 0
        for line in self.lines:
            if line[0] != "#" and len(line.strip()) > 0:
                self.data[:,i] = map(float, line.split()[:-3])
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

    def giveHippo(self, myEnzoHippo, name):
        """
        This will add a datastore to an existing HippoDraw instance.

        @param myEnzoHippo: the EnzoHippo instance to add it to as a datastore.  Note that
        this EnzoHippo instance can later go away, since we're only going to be using
        it to access the appropriate HippoApp.
        @type myEnzoHippo: L{yt.raven.EnzoHippo<yt.raven.EnzoHippo>}
        @param name: The name to pass in for the datastore
        @type name: string
        """
        pass
