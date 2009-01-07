"""
A means of extracting a subset of the hierarchy

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


from yt.mods import *
import tables, os.path

import yt.commands as commands

class DummyHierarchy(object):
    pass

class ConstructedRootGrid(object):
    def __init__(self, pf, level):
        self.pf = pf
        self.hierarchy = DummyHierarchy()
        self.hierarchy.data_style = -1
        self.Level = level
        self.LeftEdge = na.min([grid.LeftEdge for grid in
                            pf.h.select_grids(level)], axis=0).astype('float64')
        self.RightEdge = na.max([grid.RightEdge for grid in
                             pf.h.select_grids(level)], axis=0).astype('float64')
        self.index = na.min([grid.get_global_startindex() for grid in
                             pf.h.select_grids(level)], axis=0).astype('int64')
        self.dx = pf.h.select_grids(level)[0].dx
        dims = (self.RightEdge-self.LeftEdge)/self.dx
        self.ActiveDimensions = dims
        self.cg = pf.h.smoothed_covering_grid(level, self.LeftEdge,
                        self.RightEdge, dims=dims)

    def __getitem__(self, field):
        return self.cg[field]

    def get_global_startindex(self):
        return self.index

class ExtractedHierarchy(object):

    def __init__(self, pf, min_level, max_level = -1, offset = None):
        self.pf = pf
        self.min_level = min_level
        self.int_offset = na.min([grid.get_global_startindex() for grid in
                             pf.h.select_grids(min_level)], axis=0).astype('float64')
        if offset is None:
            self.left_edge_offset = na.min([grid.LeftEdge for grid in
                                       pf.h.select_grids(min_level)], axis=0).astype('float64')
        else:
            self.left_edge_offset = offset
        self.mult_factor = 2**min_level
        if max_level == -1: max_level = pf.h.max_level
        self.max_level = min(max_level, pf.h.max_level)

    def export_output(self, afile, n, field):
        # I prefer dict access, but tables doesn't.
        time_node = afile.createGroup("/", "time-%s" % n)
        time_node._v_attrs.time = self.pf["InitialTime"]
        time_node._v_attrs.numLevels = self.pf.h.max_level+1-self.min_level
        # Can take a while, so let's get a progressbar
        if len(self.pf.h.select_grids(self.min_level)) > 0:
            grids = [ConstructedRootGrid(self.pf, self.min_level)]
            self.export_level(afile, time_node, self.min_level, field, grids)
            self._export_levels(afile, time_node, field, 
                        self.min_level+1, self.max_level)
        else:
            self._export_levels(afile, time_node, field, 
                        self.min_level, self.max_level)

    def _export_levels(self, afile, time_node, field, min_level, max_level):
        pbar = yt.funcs.get_pbar("Exporting levels", max_level+1-min_level)
        for level in range(min_level, max_level+1):
            pbar.update(level-self.min_level)
            self.export_level(afile, time_node, level, field)
        pbar.finish()

    def export_level(self, afile, time_node, level, field, grids = None):
        level_node = afile.createGroup(time_node, "level-%s" % \
                                       (level-self.min_level))
        # Grid objects on this level...
        if grids is None: grids = self.pf.h.select_grids(level)
        level_node._v_attrs.delta = na.array([grids[0].dx]*3,
                                        dtype='float64')*self.mult_factor
        level_node._v_attrs.relativeRefinementFactor = na.array([2]*3, dtype='int32')
        level_node._v_attrs.numGrids = len(grids)
        for i,g in enumerate(grids):
            self.export_grid(afile, level_node, g, i, field)

    def export_grid(self, afile, level_node, grid, i, field):
        grid_node = afile.createGroup(level_node, "grid-%s" % i)
        grid_node._v_attrs.integerOrigin = (grid.get_global_startindex() \
                                         - self.int_offset*2**(grid.Level-self.min_level)).astype('int64')
        grid_node._v_attrs.origin = (grid.LeftEdge - self.left_edge_offset)*self.mult_factor
        grid_node._v_attrs.ghostzoneFlags = na.zeros(6, dtype='int32')
        grid_node._v_attrs.numGhostzones = na.zeros(3, dtype='int32')
        grid_node._v_attrs.dims = grid.ActiveDimensions[::-1].astype('int32')
        if self.pf.h.data_style == 6 and field in self.pf.h.field_list:
            if grid.hierarchy.data_style == -1: # constructed grid
                # if we can get conversion in amira we won't need to do this
                ff = grid[field]
                ff /= self.pf.conversion_factors.get(field, 1.0)
                afile.createArray(grid_node, "grid-data", ff.swapaxes(0,2))
            else:
                tfn = os.path.abspath(afile.filename)
                gfn = os.path.abspath(grid.filename)
                fpn = os.path.commonprefix([tfn, grid.filename])
                fn = grid.filename[len(os.path.commonprefix([tfn, grid.filename])):]
                grid_node._v_attrs.referenceFileName = fn
                grid_node._v_attrs.referenceDataPath = \
                    "/Grid%08i/%s" % (grid.id, field)
        else:
            # Export our array
            afile.createArray(grid_node, "grid-data", grid[field].swapaxes(0,2))

def __get_pf(bn, n):
    bn_try = "%s%04i" % (bn, n)
    return commands._fix_pf(bn_try)


def export_amira():
    parser = commands._get_parser("bn", "field", "skip")

    # We don't like the 'output' I have in the commands module
    parser.add_option("-o", "--output", action="store", type="string",
                      dest="output", default="movie.a5",
                      help="Name of our output file")
    parser.add_option("", "--minlevel", action="store", type="int",
                      dest="min_level", default=0,
                      help="The minimum level to extract (chooses first grid at that level)")
    parser.add_option("", "--maxlevel", action="store", type="int",
                      dest="max_level", default=-1,
                      help="The maximum level to extract (chooses first grid at that level)")
    parser.add_option("-d","--subtract-time", action="store_true",
                      dest="subtract_time", help="Subtract the physical time of " + \
                      "the first timestep (useful for small delta t)", default=False)
    parser.add_option("-r","--recenter", action="store_true",
                      dest="recenter", help="Recenter on maximum density in final output")
    parser.usage = "%prog [options] FIRST_ID LAST_ID"
    opts, args = parser.parse_args()

    first = int(args[0])
    last = int(args[1])

    # Set up our global metadata
    afile = tables.openFile(opts.output, "w")
    md = afile.createGroup("/", "globalMetaData")
    mda = md._v_attrs
    mda.datatype = 0
    mda.staggering = 1
    mda.fieldtype = 1

    mda.minTimeStep = first
    mda.maxTimeStep = last

    times = []
    # Get our staggering correct based on skip
    timesteps = na.arange(first, last+1, opts.skip, dtype='int32')
    time_offset = None
    t2 = []

    offset = None
    if opts.recenter:
        tpf = __get_pf(opts.basename, timesteps[-1])
        offset = tpf.h.find_max("Density")[1]
        del tpf

    for n in timesteps:
        # Try super hard to get the right parameter file
        pf = __get_pf(opts.basename, n)
        hh = pf.h
        times.append(pf["InitialTime"] * pf["years"])
        eh = ExtractedHierarchy(pf, opts.min_level, max_level = opts.max_level, offset=offset)
        eh.export_output(afile, n, opts.field)
        t2.append(pf["InitialTime"])

    # This should be the same
    mda.rootDelta = (pf["unitary"]/pf["TopGridDimensions"]).astype('float64')
    mda.minTime = times[0]
    mda.maxTime = times[-1]
    mda.numTimeSteps = len(timesteps)

    # I think we just want one value here
    rel_times = na.array(times, dtype='float64') - int(opts.subtract_time)*times[0]
    afile.createArray(md, "sorted_times", na.array(rel_times))
    afile.createArray(md, "sorted_timesteps", timesteps)

    afile.close()
