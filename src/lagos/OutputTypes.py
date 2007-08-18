"""
Generalized Enzo output objects, both static and time-series.

Presumably at some point EnzoRun will be absorbed into here.
"""

__author__ = "U{Matthew Turk<http://www.stanford.edu/~mturk>}"
__organization__ = "U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}"
__contact__ = "U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}"
__license__ = "GPL-3"


from yt.lagos import *
from yt.funcs import *
from collections import defaultdict
import string, re, gc, time, os, os.path

# We want to support the movie format in the future.
# When such a thing comes to pass, I'll move all the stuff that is contant up
# to here, and then have it instantiate EnzoStaticOutputs as appropriate.
class EnzoOutput:
    pass

class EnzoStaticOutput(EnzoOutput):
    """
    This class is a stripped down class that simply reads and parses, without
    looking at the hierarchy.

    @todo: Move some of the parameters to the EnzoRun?
           Maybe it is just more appropriate to think of time series data and
           single-time data?
    
    @param filename: The filename of the parameterfile we want to load
    @type filename: String
    """
    __hierarchy = None
    def __init__(self, filename, data_style=4):
        self.data_style = data_style
        if filename.endswith(".hierarchy"):
            filename = filename[:-10]
        self.parameterFilename = "%s" % (filename)
        self.basename = os.path.basename(filename)
        self.directory = os.path.dirname(filename)
        self.fullpath = os.path.abspath(self.directory)
        if len(self.directory) == 0:
            self.directory = "."
        self.conversionFactors = defaultdict(lambda: 1.0)
        self.parameters = {}
        self.parameters["CurrentTimeIdentifier"] = \
            int(os.stat(self.parameterFilename)[ST_CTIME])
        self.parseParameterFile()
        self.setUnits()
        # These can be taken out if you so desire
        rp = os.path.join(self.directory, "rates.out")
        if os.path.exists(rp):
            try:
                self.rates = EnzoTable(rp, rates_out_key)
            except:
                pass
        cp = os.path.join(self.directory, "cool_rates.out")
        if os.path.exists(cp):
            try:
                self.cool = EnzoTable(cp, cool_out_key)
            except:
                pass

    def getTimeID(self):
        return time.ctime(float(self["CurrentTimeIdentifier"]))

    def __xattrs__(self, mode="default"):
        return ("basename", "getTimeID()")
        
    def __getitem__(self, key):
        """
        Returns units, parameters, or conversionFactors in that order
        """
        if self.units.has_key(key):
            return self.units[key]
        elif self.timeUnits.has_key(key):
            return self.timeUnits[key]
        elif self.parameters.has_key(key):
            return self.parameters[key]
        return self.conversionFactors[key]

    def __repr__(self):
        return self.basename

    def keys(self):
        """
        Returns a list of possible keys, from units, parameters and
        conversionFactors
        """
        return self.units.keys() \
             + self.timeUnits.keys() \
             + self.parameters.keys() \
             + self.conversionFactors.keys()

    def has_key(self, key):
        """
        Returns true or false
        """
        return (self.units.has_key(key) or \
                self.timeUnits.has_key(key) or \
                self.parameters.has_key(key) or \
                self.conversionFactors.has_key(key))

    def parseParameterFile(self):
        """
        Parses the parameter file and establishes the various
        dictionaries.
        """
        # Let's read the file
        lines = open(self.parameterFilename).readlines()
        for lineI, line in enumerate(lines):
            if len(line) < 2:
                continue
            param, vals = map(strip,map(rstrip,line.split("=")))
            if parameterDict.has_key(param):
                t = map(parameterDict[param], vals.split())
                if len(t) == 1:
                    self.parameters[param] = t[0]
                else:
                    self.parameters[param] = t
                if param.endswith("Units") and not param.startswith("Temperature"):
                    dataType = param[:-5]
                    self.conversionFactors[dataType] = self.parameters[param]
            elif param.startswith("#DataCGS"):
                # Assume of the form: #DataCGSConversionFactor[7] = 2.38599e-26 g/cm^3
                if lines[lineI-1].find("Label") >= 0:
                    kk = lineI-1
                elif lines[lineI-2].find("Label") >= 0:
                    kk = lineI-2
                dataType = lines[kk].split("=")[-1].rstrip().strip()
                convFactor = float(line.split("=")[-1].split()[0])
                self.conversionFactors[dataType] = convFactor
            elif param.startswith("#CGSConversionFactor"):
                dataType = param[20:].rstrip()
                convFactor = float(line.split("=")[-1])
                self.conversionFactors[dataType] = convFactor

    def setUnits(self):
        """
        Generates the conversion to various physical units based on the parameter file
        """
        self.units = {}
        self.timeUnits = {}
        if len(self.parameters) == 0:
            self.parseParameterFile()
        if self["ComovingCoordinates"]:
            z = self["CosmologyCurrentRedshift"]
            boxh = self["CosmologyComovingBoxSize"]
            self.units['aye']  = (1.0 + self["CosmologyInitialRedshift"])/(z - 1.0)
            if not self.has_key("Time"):
                LengthUnit = 3.086e24 * boxh / self["CosmologyHubbleConstantNow"] \
                             / (1+self["CosmologyInitialRedshift"])
                self.conversionFactors["Time"] = LengthUnit / self["x-velocity"]
        elif self.has_key("LengthUnits"):
            # We are given LengthUnits, which is number of cm per box length
            # So we convert that to box-size in Mpc
            z = 0
            boxh = 3.24077e-25 * self["LengthUnits"]
            self.units['aye']  = 1.0
        else:
            z = 0
            boxh = 1.0
            self.units['aye'] = 1.0
        seconds = self["Time"]
        box = boxh/(1+z)
        for unit in unitList.keys():
            self.units[unit] = unitList[unit] * box
            self.units[unit+'h'] = unitList[unit] * boxh
        self.timeUnits['1']     = 1
        self.units['1']     = 1
        self.timeUnits['years'] = seconds / (365*3600*24.0)
        self.timeUnits['days']  = seconds / (3600*24.0)

    def _get_hierarchy(self):
        if self.__hierarchy == None:
            self.__hierarchy = EnzoHierarchy(self, data_style=self.data_style)
        return self.__hierarchy

    def _set_hierarchy(self, newh):
        if self.__hierarchy != None:
            mylog.warning("Overriding hierarchy attribute!  This is probably unwise!")
        self.__hierarchy = newh

    hierarchy = property(_get_hierarchy, _set_hierarchy)
    h = property(_get_hierarchy, _set_hierarchy)
