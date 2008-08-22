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


import yt.lagos as lagos
import yt.lagos.hop as hop
import yt.raven as raven
import yt.fido as fido
import numpy as na
import sys, types

from yt.lagos import EnzoStaticOutput, \
    BinnedProfile1D, BinnedProfile2D, BinnedProfile3D, \
    add_field, fieldInfo, \
    Clump, write_clump_hierarchy, find_clumps, write_clumps

from yt.raven import PlotCollection, PlotCollectionInteractive, \
    QuiverCallback, ParticleCallback, ContourCallback, \
    GridBoundaryCallback, UnitBoundaryCallback, \
    LinePlotCallback, CuttingQuiverCallback, ClumpContourCallback

try:
    from yt.raven import VolumeRenderingDataCube, \
        VolumeRendering3DProfile, HaloMassesPositionPlot
except ImportError:
    pass

from yt.fido import GrabCollections, OutputCollection

def get_pf():
    return lagos.EnzoStaticOutput(sys.argv[-1])

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
