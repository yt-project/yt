"""
This is a library of yt-defined exceptions

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2009 Matthew Turk.  All Rights Reserved.

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

# We don't need to import 'exceptions'
#import exceptions

class YTException(Exception):
    def __init__(self, pf = None):
        Exception.__init__(self)
        self.pf = pf

# Data access exceptions:

class YTOutputNotIdentified(YTException):
    def __init__(self, args, kwargs):
        self.args = args
        self.kwargs = kwargs

    def __str__(self):
        return "Supplied %s %s, but could not load!" % (
            self.args, self.kwargs)

class YTSphereTooSmall(YTException):
    def __init__(self, pf, radius, smallest_cell):
        YTException.__init__(self, pf)
        self.radius = radius
        self.smallest_cell = smallest_cell

    def __str__(self):
        return "%0.5e < %0.5e" % (self.radius, self.smallest_cell)

class YTAxesNotOrthogonalError(YTException):
    def __init__(self, axes):
        self.axes = axes

    def __str__(self):
        return "The supplied axes are not orthogonal.  %s" % (self.axes)

class YTNoDataInObjectError(YTException):
    def __init__(self, obj):
        self.obj_type = getattr(obj, "_type_name", "")

    def __str__(self):
        s = "The object requested has no data included in it."
        if self.obj_type == "slice":
            s += "  It may lie on a grid face.  Try offsetting slightly."
        return s

class YTFieldNotFound(YTException):
    def __init__(self, fname, pf):
        self.fname = fname
        self.pf = pf

    def __str__(self):
        return "Could not find field '%s' in %s." % (self.fname, self.pf)

class YTFieldTypeNotFound(YTException):
    def __init__(self, fname):
        self.fname = fname

    def __str__(self):
        return "Could not find field '%s'." % (self.fname)

class YTSimulationNotIdentified(YTException):
    def __init__(self, sim_type):
        YTException.__init__(self)
        self.sim_type = sim_type

    def __str__(self):
        return "Simulation time-series type %s not defined." % self.sim_type

class AmbiguousOutputs(YTException):
    def __init__(self, pf):
        YTException.__init__(self, pf)

    def __str__(self):
        return "Simulation %s has both dtDataDump and CycleSkipDataDump set.  Unable to calculate datasets." % \
            self.pf

class MissingParameter(YTException):
    def __init__(self, pf, parameter):
        YTException.__init__(self, pf)
        self.parameter = parameter

    def __str__(self):
        return "Parameter file %s is missing %s parameter." % \
            (self.pf, self.parameter)

class NoStoppingCondition(YTException):
    def __init__(self, pf):
        YTException.__init__(self, pf)

    def __str__(self):
        return "Simulation %s has no stopping condition.  StopTime or StopCycle should be set." % \
            self.pf

class YTNotInsideNotebook(YTException):
    def __str__(self):
        return "This function only works from within an IPython Notebook."

class YTNotDeclaredInsideNotebook(YTException):
    def __str__(self):
        return "You have not declared yourself to be inside the IPython" + \
               "Notebook.  Do so with this command:\n\n" + \
               "ytcfg['yt','ipython_notebook'] = 'True'"

class YTGeometryNotSupported(YTException):
    def __init__(self, geom):
        self.geom = geom

    def __str__(self):
        return "We don't currently support %s geometry" % self.geom

class YTCoordinateNotImplemented(YTException):
    def __str__(self):
        return "This coordinate is not implemented for this geometry type."

class YTUnitNotRecognized(YTException):
    def __init__(self, unit):
        self.unit = unit

    def __str__(self):
        return "This parameter file doesn't recognize %s" % self.unit

class YTHubRegisterError(YTException):
    def __str__(self):
        return "You must create an API key before uploading.  See " + \
               "https://data.yt-project.org/getting_started.html"
