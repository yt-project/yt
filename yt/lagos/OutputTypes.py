"""
Generalized Enzo output objects, both static and time-series.

Presumably at some point EnzoRun will be absorbed into here.
@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
@license:
  Copyright (C) 2007 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

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
    def __init__(self, filename, data_style=None):
        self.data_style = data_style
        if filename.endswith(".hierarchy"):
            filename = filename[:-10]
        self.parameter_filename = "%s" % (filename)
        self.basename = os.path.basename(filename)
        self.directory = os.path.dirname(filename)
        self.fullpath = os.path.abspath(self.directory)
        if len(self.directory) == 0:
            self.directory = "."
        #self.conversion_factors = defaultdict(lambda: 1.0)
        self.conversion_factors = {}
        self.parameters = {}
        self.parameters["CurrentTimeIdentifier"] = \
            int(os.stat(self.parameter_filename)[ST_CTIME])
        self.__parse_parameter_file()
        self.__set_units()
        # These can be taken out if you so desire
        rp = os.path.join(self.directory, "rates.out")
        if os.path.exists(rp):
            self.rates = EnzoTable(rp, rates_out_key)
        cp = os.path.join(self.directory, "cool_rates.out")
        if os.path.exists(cp):
            self.cool = EnzoTable(cp, cool_out_key)

    def getTimeID(self):
        return time.ctime(float(self["CurrentTimeIdentifier"]))

    def __xattrs__(self, mode="default"):
        return ("basename", "getTimeID()")

    def __getitem__(self, key):
        """
        Returns _units, parameters, or _conversion_factors in that order
        """
        if self.units.has_key(key):
            return self.units[key]
        elif self.time_units.has_key(key):
            return self.time_units[key]
        elif self.parameters.has_key(key):
            return self.parameters[key]
        return self.conversion_factors[key]

    def __repr__(self):
        return self.basename

    def keys(self):
        """
        Returns a list of possible keys, from _units, parameters and
        _conversion_factors
        """
        return self.units.keys() \
             + self.time_units.keys() \
             + self.parameters.keys() \
             + self.conversion_factors.keys()

    def has_key(self, key):
        """
        Returns true or false
        """
        return (self.units.has_key(key) or \
                self.time_units.has_key(key) or \
                self.parameters.has_key(key) or \
                self.conversion_factors.has_key(key))

    def __parse_parameter_file(self):
        """
        Parses the parameter file and establishes the various
        dictionaries.
        """
        # Let's read the file
        lines = open(self.parameter_filename).readlines()
        for lineI, line in enumerate(lines):
            if line.find("#") >= 1: # Keep the commented lines
                line=line[:line.find("#")]
            line=line.strip().rstrip()
            if len(line) < 2:
                continue
            try:
                param, vals = map(strip,map(rstrip,line.split("=")))
            except ValueError:
                mylog.error("ValueError: '%s'", line)
            if parameterDict.has_key(param):
                t = map(parameterDict[param], vals.split())
                if len(t) == 1:
                    self.parameters[param] = t[0]
                else:
                    self.parameters[param] = t
                if param.endswith("Units") and not param.startswith("Temperature"):
                    dataType = param[:-5]
                    self.conversion_factors[dataType] = self.parameters[param]
            elif param.startswith("#DataCGS"):
                # Assume of the form: #DataCGSConversionFactor[7] = 2.38599e-26 g/cm^3
                if lines[lineI-1].find("Label") >= 0:
                    kk = lineI-1
                elif lines[lineI-2].find("Label") >= 0:
                    kk = lineI-2
                dataType = lines[kk].split("=")[-1].rstrip().strip()
                convFactor = float(line.split("=")[-1].split()[0])
                self.conversion_factors[dataType] = convFactor
            elif param.startswith("#CGSConversionFactor"):
                dataType = param[20:].rstrip()
                convFactor = float(line.split("=")[-1])
                self.conversion_factors[dataType] = convFactor
            elif param.startswith("DomainRightEdge"):
                self.parameters["DomainRightEdge"] = \
                    na.array([float(i) for i in vals.split()])

    def __set_units(self):
        """
        Generates the conversion to various physical _units based on the parameter file
        """
        self.units = {}
        self.time_units = {}
        if len(self.parameters) == 0:
            self.__parse_parameter_file()
        if self["ComovingCoordinates"]:
            self.__setup_comoving_units()
        elif self.has_key("LengthUnits"):
            self.__setup_getunits_units()
        else:
            self.__setup_nounits_units()
        self.time_units['1'] = 1
        self.units['1'] = 1
        seconds = self["Time"]
        self.time_units['years'] = seconds / (365*3600*24.0)
        self.time_units['days']  = seconds / (3600*24.0)

    def __setup_comoving_units(self):
        z = self["CosmologyCurrentRedshift"]
        h = self["CosmologyHubbleConstantNow"]
        boxcm_cal = self["CosmologyComovingBoxSize"]
        boxcm_uncal = boxcm_cal / h
        box_proper = boxcm_uncal/(1+z)
        self.units['aye']  = (1.0 + self["CosmologyInitialRedshift"])/(z - 1.0)
        if not self.has_key("Time"):
            LengthUnit = 3.086e24 * box_proper
            self.conversion_factors["Time"] = LengthUnit / self["x-velocity"]
        for unit in unitList.keys():
            self.units[unit] = unitList[unit] * box_proper
            self.units[unit+'h'] = unitList[unit] * box_proper * h
            self.units[unit+'hcm'] = unitList[unit] * boxcm_cal

    def __setup_getunits_units(self):
        # We are given LengthUnits, which is number of cm per box length
        # So we convert that to box-size in Mpc
        box_proper = 3.24077e-25 * self["LengthUnits"]
        self.units['aye']  = 1.0
        for unit in unitList.keys():
            self.units[unit] = unitList[unit] * box_proper

    def __setup_nounits_units(self):
        z = 0
        box_proper = 1.0
        self.units['aye'] = 1.0
        mylog.warning("No length units.  Setting 1.0 = 1 proper Mpc.")
        if not self.has_key("TimeUnits"):
            mylog.warning("No time units.  Setting 1.0 = 1 second.")
            self.conversion_factors["Time"] = 1.0
        for unit in unitList.keys():
            self.units[unit] = unitList[unit] * box_proper
        return box, boxh

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

    def cosmology_get_units(self):
        k = {}
        k["utim"] = 2.52e17/na.sqrt(self.parameters["CosmologyOmegaMatterNow"])\
                       / self.parameters["CosmologyHubbleConstantNow"] \
                       / (1+self.parameters["CosmologyInitialRedshift"])**1.5
        k["urho"] = 1.88e-29 * self.parameters["CosmologyOmegaMatterNow"] \
                        * self.parameters["CosmologyHubbleConstantNow"]**2 \
                        * (1.0 + self.parameters["CosmologyCurrentRedshift"])**3
        k["uxyz"] = 3.086e24 * \
               self.parameters["CosmologyComovingBoxSize"] / \
               self.parameters["CosmologyHubbleConstantNow"] / \
               (1.0 + self.parameters["CosmologyCurrentRedshift"])
        k["uaye"] = 1.0/(1.0 + self.parameters["CosmologyInitialRedshift"])
        k["uvel"] = 1.225e7*self.parameters["CosmologyComovingBoxSize"] \
                      *na.sqrt(self.parameters["CosmologyOmegaMatterNow"]) \
                      *na.sqrt(1+ self.parameters["CosmologyInitialRedshift"])
        k["utem"] = 1.88e6 * (self.parameters["CosmologyComovingBoxSize"]**2) \
                      * self.parameters["CosmologyOmegaMatterNow"] \
                      * (1.0 + self.parameters["CosmologyInitialRedshift"])
        k["aye"]  = (1.0 + self.parameters["CosmologyInitialRedshift"]) / \
               (1.0 + self.parameters["CosmologyCurrentRedshift"])
        return k
