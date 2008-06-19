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
import string, re, gc, time, os, os.path

# We want to support the movie format in the future.
# When such a thing comes to pass, I'll move all the stuff that is contant up
# to here, and then have it instantiate EnzoStaticOutputs as appropriate.

class StaticOutput(object):
    def __init__(self, filename, data_style=None):
        """
        Base class for generating new output types.  Principally consists of
        a *filename* and a *data_style* which will be passed on to children.
        """
        self.data_style = data_style
        self.parameter_filename = str(filename)
        self.basename = os.path.basename(filename)
        self.directory = os.path.expanduser(os.path.dirname(filename))
        self.fullpath = os.path.abspath(self.directory)
        if len(self.directory) == 0:
            self.directory = "."
        self.conversion_factors = {}
        self.parameters = {}
        self._parse_parameter_file()
        self._set_units()
        # These can be taken out if you so desire

    def __repr__(self):
        return self.basename

    def __getitem__(self, key):
        """
        Returns _units, parameters, or _conversion_factors in that order
        """
        for d in [self.units, self.time_units, self.parameters, \
                  self.conversion_factors]:
            if key in d: return d[key]
        raise KeyError(key)

    def keys(self):
        """
        Returns a list of possible keys, from _units, parameters and
        _conversion_factors
        """
        return self.units.keys() \
             + self.time_units.keys() \
             + self.parameters.keys() \
             + self.conversion_factors.keys()

    def get_smallest_appropriate_unit(self, v):
        max_nu = 1e30
        good_u = None
        for unit in ['mpc','kpc','pc','au','rsun','cm']:
            vv = v*self[unit]
            if vv < max_nu and vv > 1.0:
                good_u = unit
                max_nu = v*self[unit]
        return good_u

    def has_key(self, key):
        """
        Returns true or false
        """
        return key in self.units or \
               key in self.time_units or \
               key in self.parameters or \
               key in self.conversion_factors

    def _get_hierarchy(self):
        if self.__hierarchy == None:
            if self._hierarchy_class == None:
                raise RuntimeError("You should not instantiate StaticOutput.")
            self.__hierarchy = self._hierarchy_class(self, data_style=self.data_style)
        return self.__hierarchy

    def _set_hierarchy(self, newh):
        if self.__hierarchy != None:
            mylog.warning("Overriding hierarchy attribute!  This is probably unwise!")
        self.__hierarchy = newh

    __hierarchy = None
    hierarchy = property(_get_hierarchy, _set_hierarchy)
    h = property(_get_hierarchy, _set_hierarchy)


class EnzoStaticOutput(StaticOutput):
    """
    Enzo-specific output, set at a fixed time.
    """
    _hierarchy_class = EnzoHierarchy
    def __init__(self, filename, data_style=None,
                 parameter_override = None,
                 conversion_override = None):
        """
        This class is a stripped down class that simply reads and parses
        *filename* without looking at the hierarchy.  *data_style* gets passed
        to the hierarchy to pre-determine the style of data-output.  However,
        it is not strictly necessary.  Optionally you may specify a
        *parameter_override* dictionary that will override anything in the
        paarmeter file and a *conversion_override* dictionary that consists
        of {fieldname : conversion_to_cgs} that will override the #DataCGS.
        """
        if parameter_override is None: parameter_override = {}
        self.__parameter_override = parameter_override
        if conversion_override is None: conversion_override = {}
        self.__conversion_override = conversion_override

        StaticOutput.__init__(self, filename, data_style)
        rp = os.path.join(self.directory, "rates.out")
        if os.path.exists(rp):
            self.rates = EnzoTable(rp, rates_out_key)
        cp = os.path.join(self.directory, "cool_rates.out")
        if os.path.exists(cp):
            self.cool = EnzoTable(cp, cool_out_key)

    def _parse_parameter_file(self):
        """
        Parses the parameter file and establishes the various
        dictionaries.
        """
        # Let's read the file
        self.parameters["CurrentTimeIdentifier"] = \
            int(os.stat(self.parameter_filename)[ST_CTIME])
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
        for p, v in self.__parameter_override.items():
            self.parameters[p] = v
        for p, v in self.__conversion_override.items():
            self.conversion_factors[p] = v

    def _set_units(self):
        """
        Generates the conversion to various physical _units based on the parameter file
        """
        self.units = {}
        self.time_units = {}
        if len(self.parameters) == 0:
            self._parse_parameter_file()
        if self["ComovingCoordinates"]:
            self._setup_comoving_units()
        elif self.has_key("LengthUnits"):
            self._setup_getunits_units()
        else:
            self._setup_nounits_units()
        self.time_units['1'] = 1
        self.units['1'] = 1
        seconds = self["Time"]
        self.time_units['years'] = seconds / (365*3600*24.0)
        self.time_units['days']  = seconds / (3600*24.0)

    def _setup_comoving_units(self):
        z = self["CosmologyCurrentRedshift"]
        h = self["CosmologyHubbleConstantNow"]
        boxcm_cal = self["CosmologyComovingBoxSize"]
        boxcm_uncal = boxcm_cal / h
        box_proper = boxcm_uncal/(1+z)
        self.units['aye']  = (1.0 + self["CosmologyInitialRedshift"])/(z + 1.0)
        if not self.has_key("Time"):
            LengthUnit = 3.086e24 * box_proper
            self.conversion_factors["Time"] = LengthUnit / self["x-velocity"]
        for unit in mpc_conversion:
            self.units[unit] = mpc_conversion[unit] * box_proper
            self.units[unit+'h'] = mpc_conversion[unit] * box_proper * h
            self.units[unit+'hcm'] = mpc_conversion[unit] * boxcm_cal

    def _setup_getunits_units(self):
        # We are given LengthUnits, which is number of cm per box length
        # So we convert that to box-size in Mpc
        box_proper = 3.24077e-25 * self["LengthUnits"]
        self.units['aye']  = 1.0
        for unit in mpc_conversion.keys():
            self.units[unit] = mpc_conversion[unit] * box_proper

    def _setup_nounits_units(self):
        z = 0
        box_proper = ytcfg.getfloat("lagos","nounitslength")
        self.units['aye'] = 1.0
        mylog.warning("No length units.  Setting 1.0 = %0.3e proper Mpc.", box_proper)
        if not self.has_key("TimeUnits"):
            mylog.warning("No time units.  Setting 1.0 = 1 second.")
            self.conversion_factors["Time"] = 1.0
        for unit in mpc_conversion.keys():
            self.units[unit] = mpc_conversion[unit] * box_proper

    def cosmology_get_units(self):
        """
        Return an Enzo-fortran style dictionary of units to feed into custom
        routines.  This is typically only necessary if you are interacting
        with fortran code.
        """
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
