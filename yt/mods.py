"""
Very simple convenience function for importing all the modules, setting up
the namespace and getting the last argument on the command line.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008 Matthew Turk.  All Rights Reserved.

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

#
# ALL IMPORTS GO HERE
#

# First module imports
import yt.lagos as lagos
import yt.lagos.hop as hop
import yt.raven as raven
import yt.fido as fido
import numpy as na
import sys, types
from logger import ytLogger as mylog

# Now individual component imports from lagos
from yt.lagos import EnzoStaticOutput, \
    BinnedProfile1D, BinnedProfile2D, BinnedProfile3D, \
    add_field, FieldInfo, EnzoFieldInfo, Enzo2DFieldInfo, OrionFieldInfo, \
    Clump, write_clump_hierarchy, find_clumps, write_clumps

# This is a temporary solution -- in the future, we will allow the user to
# select this via ytcfg.

fieldInfo = lagos.fieldInfo = EnzoFieldInfo

# Now individual component imports from raven
from yt.raven import PlotCollection, PlotCollectionInteractive, \
    QuiverCallback, ParticleCallback, ContourCallback, \
    GridBoundaryCallback, UnitBoundaryCallback, \
    LinePlotCallback, CuttingQuiverCallback, ClumpContourCallback, \
    HopCircleCallback

# Optional component imports from raven
try:
    from yt.raven import VolumeRenderingDataCube, \
        VolumeRendering3DProfile, HaloMassesPositionPlot
except ImportError:
    pass

# Individual imports from Fido
from yt.fido import GrabCollections, OutputCollection

# Some convenience functions to ease our time running scripts
# from the command line

def get_pf():
    return EnzoStaticOutput(sys.argv[-1])

def get_pc():
    return PlotCollection(EnzoStaticOutput(sys.argv[-1]))

# Now the | operator overloading
# (which is totally a stunt)

class _StaticOutputIfier(object):
    def __init__(self):
        pass
    def __ror__(self, other):
        return EnzoStaticOutput(other)
static = _StaticOutputIfier()

class __PlotVM(object):
    def __init__(self, axis = 0, field = "Density", name = None,
                width = None, **kwargs):
        self.axis = axis
        self.field = field
        self.name = name
        self.width = width
        self.kwargs = kwargs

    def __ror__(self, pf):
        if isinstance(pf, types.StringTypes): pf = EnzoStaticOutput(pf)
        pc = PlotCollection(pf)
        self._add_plot(pc)
        if self.name is None: self.name = str(pf)
        if self.width is not None:
            pc.set_width(*self.width)
        return pc.save(self.name)

    def __call__(self, *args, **kwargs):
        return type(self)(*args, **kwargs)

class _PlotSlice(__PlotVM):
    def _add_plot(self, pc):
        pc.add_slice(self.field, self.axis, **self.kwargs)
x_slicer = _PlotSlice(axis=0)
y_slicer = _PlotSlice(axis=1)
z_slicer = _PlotSlice(axis=2)

class _PlotProj(__PlotVM):
    def _add_plot(self, pc):
        pc.add_projection(self.field, self.axis, **self.kwargs)
x_projector = _PlotProj(axis=0)
y_projector = _PlotProj(axis=1)
z_projector = _PlotProj(axis=2)

class __MultiPlotter(object):
    def __init__(self, field = "Density", name = None, width = None, **kwargs):
        self.field = field
        self.name = name
        self.width = width
        self.kwargs = kwargs

    def __ror__(self, pf):
        if isinstance(pf, types.StringTypes): pf = EnzoStaticOutput(pf)
        pc = PlotCollection(pf)
        self._add_plot(pc)
        if self.name is None: self.name = str(pf)
        if self.width is not None:
            pc.set_width(*self.width)
        return pc.save(self.name)

    def __call__(self, *args, **kwargs):
        return type(self)(*args, **kwargs)

class _MultiPlotSlice(__PlotVM):
    def _add_plot(self, pc):
        for ax in range(3): pc.add_slice(self.field, ax, **self.kwargs)
slicer = _MultiPlotSlice()

class _MultiPlotProj(__PlotVM):
    def _add_plot(self, pc):
        for ax in range(3): pc.add_projection(self.field, ax, **self.kwargs)
projector = _MultiPlotProj()

# Now some recipes
#
# NOTE HIGH LEVEL OF MAGIC.
# This is not for advanced users.

def _get_current_pf():
    # We continue until we have 'pf' in the locals space
    import inspect
    for s in inspect.stack()[1:]:
        if 'pf' in s[0].f_locals:
            __pf = s[0].f_locals['pf']
            if isinstance(__pf, types.StringTypes):
                __pf = EnzoStaticOutput(__pf)
            mylog.info("Obtained parameter file %s", __pf)
            return __pf
    
def hop_plot(my_pf = None):
    if my_pf is None: my_pf = _get_current_pf()
    pc = PlotCollection(my_pf, center=[0.5,0.5,0.5])
    center = (my_pf["DomainRightEdge"]-my_pf["DomainLeftEdge"])/2.0
    hop_output = hop.HopList(my_pf.h.sphere(center, 1.0/my_pf["1"]))
    hop_output.write_out("%s.hop" % my_pf)
    for ax in range(3):
        pc.add_projection("Density", ax).add_callback(
                            HopCircleCallback(hop_output, ax))
    pc.save("%s_hop" % my_pf)
