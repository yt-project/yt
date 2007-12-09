"""
Various non-grid data containers.

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

class EnzoData:
    """
    Generic EnzoData container.  By itself, will attempt to
    generate field, read fields (method defined by derived classes)
    """
    _grids = None
    _num_ghost_zones = 0

    def __init__(self, pf, fields):
        """
        Sets up EnzoData instance

        @param hierarchy: hierarchy we're associated with
        @type hierarchy: L{EnzoHierarchy<EnzoHierarchy>}
        @param fields: Fields represented in the data
        @type fields: list of strings
        """
        if pf != None:
            self.pf = pf
            self.hierarchy = pf.hierarchy
        if fields == None: fields = []
        self.fields = ensure_list(fields)
        self.data = {}
        self.field_parameters = {}
        self.__set_default_field_parameters()

    def __set_default_field_parameters(self):
        self.set_field_parameter("center",na.zeros(3,dtype='float64'))
        self.set_field_parameter("bulk_velocity",na.zeros(3,dtype='float64'))

    def get_field_parameter(self, name, default=None):
        if self.field_parameters.has_key(name):
            return self.field_parameters[name]
        else:
            return default

    def set_field_parameter(self, name, val):
        self.field_parameters[name] = val

    def convert(self, datatype):
        return self.hierarchy[datatype]

    def clear_data(self):
        """
        Clears out all data from the EnzoData instance, freeing memory.
        """
        del self.data
        self.data = {}

    def has_key(self, key):
        return self.data.has_key(key)

    def _refresh_data(self):
        """
        Wipes data and rereads/regenerates it from the self.fields.
        """
        self.clear_data()
        self.get_data()

    def __getitem__(self, key):
        """
        Returns a single field.  Will add if necessary.
        """
        if not self.data.has_key(key):
            if key not in self.fields:
                self.fields.append(key)
            self.get_data(key)
        return self.data[key]

    def __setitem__(self, key, val):
        """
        Sets a field to be some other value.
        """
        if key not in self.fields: self.fields.append(key)
        self.data[key] = val

    def __delitem__(self, key):
        """
        Sets a field to be some other value.
        """
        try:
            del self.fields[self.fields.index(key)]
        except ValueError:
            pass
        del self.data[key]

    def _generate_field_in_grids(self, fieldName):
        pass

class Enzo2DData(EnzoData):
    """
    Class to represent a set of EnzoData that's 2-D in nature.  Slices and
    projections, primarily.
    """
    _spatial = False
    def __init__(self, axis, fields, pf=None):
        """
        Prepares the Enzo2DData.

        @param hierarchy: hierarchy associated with this data
        @type hierarchy: L{EnzoHierarchy<EnzoHierarchy>}
        @param axis: axis to which this data is parallel
        @type axis: integer
        @param fields: fields to be processed or generated
        @type fields: list of strings
        """
        self.axis = axis
        EnzoData.__init__(self, pf, fields)

    def _generate_field(self, field):
        """
        Generates, or attempts to generate, a field not found in the data file

        @param field: field to generate
        @type field: string
        """
        if fieldInfo.has_key(field):
            # First we check the validator
            try:
                fieldInfo[field].check_available(self)
            except NeedsGridType, ngt_exception:
                # We leave this to be implementation-specific
                self._generate_field_in_grids(field, ngt_exception.ghost_zones)
                return False
            else:
                self[field] = fieldInfo[field](self)
                return True
        else: # Can't find the field, try as it might
            raise exceptions.KeyError(field)


    def interpolate_discretize(self, LE, RE, field, side, logSpacing=True):
        """
        This returns a uniform grid of points, interpolated using the nearest
        neighbor method.

        @note: Requires Delaunay triangulation, which is not included in
        most/all scipy binaries.
        @param LE: Left Edge
        @type LE: array of Floats
        @param RE: Right Edge
        @type RE: array of Floats
        @param field: The field to discretize
        @type field: string
        @param side: The number of points on a side
        @type side: int
        """
        import scipy.sandbox.delaunay as de
        if logSpacing:
            zz = na.log10(self[field])
        else:
            zz = self[field]
        xi, yi = na.array( \
                 na.mgrid[LE[0]:RE[0]:side*1j, \
                          LE[1]:RE[1]:side*1j], 'float64')
        zi = de.Triangulation(self['px'],self['py']).nn_interpolator(zz)\
                 [LE[0]:RE[0]:side*1j, \
                  LE[1]:RE[1]:side*1j]
        if logSpacing:
            zi = 10**(zi)
        return [xi,yi,zi]
        #return [xx,yy,zz]

class EnzoSliceBase(Enzo2DData):
    """
    A slice at a given coordinate along a given axis through the entire grid
    hierarchy.
    """

    @time_execution
    def __init__(self, axis, coord, fields = None, center=None, pf=None):
        """
        We are not mandating a field be passed in
        The field and coordinate we want to be able to change -- however, the
        axis we do NOT want to change.  So think of EnzoSlice as defining a slice
        operator, rather than a set piece of data.

        Arguments:
        @param hierarchy: hierarchy associated with this data
        @type hierarchy: L{EnzoHierarchy<EnzoHierarchy>}
        @param axis: axis to which this data is parallel
        @type axis: integer
        @param coord: three points defining the center
        @type coord: na.array
        @keyword fields: fields to be processed or generated
        @type fields: list of strings
        """
        Enzo2DData.__init__(self, axis, fields, pf)
        self.center = center
        self.coord = coord
        self._refresh_data()

    def reslice(self, coord):
        mylog.debug("Setting coordinate to %0.5e" % coord)
        self.coord = coord
        self._refresh_data()

    def shift(self, val):
        """
        Moves the slice coordinate up by either a floating point value, or an
        integer number of indices of the finest grid.

        @param val: shift amount
        @type val: integer (number of cells) or float (distance)
        """
        if isinstance(val, types.FloatType):
            # We add the dx
            self.coord += val
        elif isinstance(val, types.IntType):
            # Here we assume that the grid is the max level
            level = self.hierarchy.max_level
            self.coord
            dx = self.hierarchy.gridDxs[self.hierarchy.levelIndices[level][0]]
            self.coord += dx * val
        else:
            raise ValueError(val)
        self.refreshData()

    def _generate_coords(self):
        points = []
        for grid in self._grids:
            points.append(self._generate_grid_coords(grid))
        t = na.concatenate(points)
        self['px'] = t[:,0]
        self['py'] = t[:,1]
        self['pz'] = t[:,2]
        self['pdx'] = t[:,3]
        self['pdy'] = t[:,3]
        self['pdz'] = t[:,3]

        # Now we set the *actual* coordinates
        self[axis_names[x_dict[self.axis]]] = t[:,0]
        self[axis_names[y_dict[self.axis]]] = t[:,1]
        self[axis_names[self.axis]] = t[:,2]
        self['dx'] = t[:,3]
        self['dy'] = t[:,3]
        self['dz'] = t[:,3]

        self.ActiveDimensions = (t.shape[0], 1, 1)

    @time_execution
    def get_data(self, fields = None):
        """
        Iterates over the list of fields and generates/reads them all.

        @keyword field: field to add (list or single)
        @type field: string or list of strings
        """
        # We get it for the values in fields and coords
        # We take a 3-tuple of the coordinate we want to slice through, as well
        # as the axis we're slicing along
        self._grids, ind = self.hierarchy.find_slice_grids(self.coord, self.axis)
        if not self.has_key('dx'):
            self._generate_coords()
        if fields == None:
            fields_to_get = self.fields
        else:
            fields_to_get = ensure_list(fields)
        for field in fields_to_get:
            if self.data.has_key(field):
                continue
            rvs=[]
            if field not in self.hierarchy.fieldList:
                if self._generate_field(field):
                    continue # A "True" return means we did it
            self[field] = na.concatenate(
                [self._get_data_from_grid(grid, field)
                 for grid in self._grids])

    def _generate_grid_coords(self, grid):
        xaxis = x_dict[self.axis]
        yaxis = y_dict[self.axis]
        wantedIndex = int(((self.coord-grid.LeftEdge[self.axis])/grid.dx))
        sl = [slice(None), slice(None), slice(None)]
        sl[self.axis] = slice(wantedIndex, wantedIndex + 1)
        #sl.reverse()
        sl = tuple(sl)
        nx = grid.child_mask.shape[xaxis]
        ny = grid.child_mask.shape[yaxis]
        cm = na.where(grid.child_mask[sl].ravel() == 1)
        cmI = na.indices((nx,ny))
        xind = cmI[0,:].ravel()
        xpoints = na.ones(cm[0].shape, 'float64')
        xpoints *= xind[cm]*grid.dx+(grid.LeftEdge[xaxis] + 0.5*grid.dx)
        yind = cmI[1,:].ravel()
        ypoints = na.ones(cm[0].shape, 'float64')
        ypoints *= yind[cm]*grid.dx+(grid.LeftEdge[yaxis] + 0.5*grid.dx)
        zpoints = na.ones(xpoints.shape, 'float64') * self.coord
        dx = na.ones(xpoints.shape, 'float64') * grid.dx/2.0
        t = na.array([xpoints, ypoints, zpoints, dx]).swapaxes(0,1)
        return t

    def _get_data_from_grid(self, grid, field):
        # So what's our index of slicing?  This is what we need to figure out
        # first, so we can deal with our data in the fastest way.
        wantedIndex = int(((self.coord-grid.LeftEdge[self.axis])/grid.dx))
        sl = [slice(None), slice(None), slice(None)]
        sl[self.axis] = slice(wantedIndex, wantedIndex + 1)
        slHERE = tuple(sl)
        sl.reverse()
        slHDF = tuple(sl)
        if not grid.has_key(field):
            conv_factor = 1.0
            if fieldInfo.has_key(field):
                conv_factor = fieldInfo[field]._convert_function(self)
            dv = self._read_data_slice(grid, field, slHDF) * conv_factor
        else:
            dv = grid[field][slHERE]
        cm = na.where(grid.child_mask[slHERE].ravel() == 1)
        dataVals = dv.ravel()[cm]
        return dataVals

    def _generate_field_in_grids(self, field, num_ghost_zones=0):
        for grid in self._grids:
            temp = grid[field]

class EnzoCuttingPlaneBase(Enzo2DData):
    _plane = None
    def __init__(self, normal, point, fields = None):
        Enzo2DData.__init__(self, 4, fields, pf)
        # Let's set up our plane equation
        # ax + by + cz + d = 0
        self._norm_vec = normal/na.sqrt(na.dot(normal,normal))
        self._d = na.dot(self._norm_vec, point)

    def _identify_grids(self):
        # We want to examine each grid and see if any opposing corner vertices
        # are on opposite sides of the plane
        # Check these pairs: [ (LLL,RRR), (LLR,RRL), (LRR,RLL), (LRL,RLR) ]
        LE = self.pf.h.gridLeftEdge
        RE = self.pf.h.gridRightEdge
        vertices = na.array([[[LE[:,0],LE[:,1],LE[:,2]],
                              [RE[:,0],RE[:,1],RE[:,2]]],
                             [[LE[:,0],LE[:,1],RE[:,2]],
                              [RE[:,0],RE[:,1],LE[:,2]]],
                             [[LE[:,0],RE[:,1],RE[:,2]],
                              [RE[:,0],LE[:,1],LE[:,2]]],
                             [[LE[:,0],RE[:,1],LE[:,2]],
                              [RE[:,0],LE[:,1],RE[:,2]]]])
        # This gives us shape: 4, 2, 3, n_grid
        intersect
        #for i in [0,1,2,3]:
            # Iterate over each pair


class EnzoProjBase(Enzo2DData):
    def __init__(self, axis, field, weight_field = None,
                 max_level = None, center = None, pf = None,
                 type=0):
        """
        Returns an instance of EnzoProj.

        @todo: Set up for multiple fields
        @param hierarchy: the hierarchy we are projecting
        @type hierarchy: L{EnzoHierarchy<EnzoHierarchy>}
        @param axis: axis to project along
        @type axis: integer
        @param field: the field to project (NOT multiple)
        @type field: string
        @keyword weightField: the field to weight by
        @keyword max_level: the maximum level to project through
        @keyword type: The type of projection: 0 for sum, 1 for MIP
        """
        Enzo2DData.__init__(self, axis, field, pf)
        if max_level == None:
            max_level = self.hierarchy.maxLevel
        self._max_level = max_level
        self._weight = weight_field
        self.center = center
        if type == 1:
            self.type="MIP"
            self.func = na.max
        else:
            self.type="SUM"
            self.func = na.sum
        self.__calculate_memory()
        self._refresh_data()

    @time_execution
    def __calculate_memory(self):
        """
        Here we simply calculate how much memory is needed, which speeds up
        allocation later.  Additionally, overlap masks get pre-generated.
        """
        level_mem = {}
        h = self.hierarchy
        i = 0
        mylog.info("Calculating memory usage")
        for level in range(self._max_level+1):
            level_mem[level] = 0
            mylog.debug("Examining level %s", level)
            grids = h.levelIndices[level]
            num_grids = len(grids)
            RE = h.gridRightEdge[grids].copy()
            LE = h.gridLeftEdge[grids].copy()
            for grid in h.grids[grids]:
                if (i%1e3) == 0:
                    mylog.debug("Reading and masking %s / %s", i, h.num_grids)
                for ax in [0,1,2]:
                    grid._generate_overlap_masks(ax, LE, RE)
                    grid._overlap_grids[ax] = \
                      h.grids[grids[na.where(grid.overlap_masks[ax] == 1)]]
                level_mem[level] += \
                          grid.ActiveDimensions.prod() / \
                          grid.ActiveDimensions[ax]
                i += 1
        for level in range(self._max_level+1):
            gI = h.levelIndices[level]
            mylog.debug("%s cells and %s grids for level %s", \
                level_mem[level], len(gI), level)
        mylog.debug("We need %s cells total",
                    na.add.reduce(level_mem.values()))
        self.__memory_per_level = level_mem

    def _serialize(self):
        mylog.info("Serializing data...")
        nodeName = "%s_%s_%s" % (self.fields[0], self.weightField, self.axis)
        mylog.info("nodeName: %s", nodeName)
        projArray = na.array([self.x, self.y, self.dx, self.dy, self[self.fields[0]]])
        self.hierarchy.saveData(projArray, "/Projections", nodeName)
        mylog.info("Done serializing...")

    def _deserialize(self):
        nodeName = "%s_%s_%s" % (self.fields[0], self.weightField, self.axis)
        array=self.hierarchy.getData("/Projections", nodeName)
        if array == None:
            return
        self['px'] = array[0,:]
        self['py'] = array[1,:]
        self['pdx'] = array[2,:]
        self['pdy']= array[3,:]
        self[self.fields[0]] = array[4,:]
        return True

    def __project_level(self, level, field):
        grids_to_project = self.hierarchy.select_grids(level)
        zeroOut = (level != self._max_level)
        pbar = get_pbar('Projecting level % 2i / % 2i ' \
                          % (level, self._max_level), len(grids_to_project))
        for pi, grid in enumerate(grids_to_project):
            grid.retVal = grid._get_projection(self.axis, field, zeroOut,
                                               weight=self._weight,
                                               func = self.func)
            pbar.update(pi)
        pbar.finish()
        self.__combine_grids(level) # In-place
        all_data = [ [grid.retVal[j] for grid in grids_to_project] for j in range(5)]
        for grid in grids_to_project:
            cI = na.where(grid.retVal[3]==0) # Where childmask = 0
            grid.coarseData = [grid.retVal[j][cI] for j in range(5)]
        levelData = [na.concatenate(all_data[i]) for i in range(5)]
        mylog.debug("All done combining and refining with a final %s points",
                    levelData[0].shape[0])
        dblI = na.where((levelData[0]>-1) & (levelData[3]==1))
        if self._weight != None:
            weightedData = levelData[2][dblI] / levelData[4][dblI]
        else:
            weightedData = levelData[2][dblI]
        mylog.debug("Level %s done: %s final of %s", \
                   level, len(dblI[0]), \
                   levelData[0].shape[0])
        dx = grids_to_project[0].dx * na.ones(len(dblI[0]), dtype='float64')
        return [levelData[0][dblI], levelData[1][dblI], weightedData, dx]

    def __combine_grids(self, level):
        grids = self.hierarchy.select_grids(level)
        pbar = get_pbar('Combining level % 2i / % 2i ' \
                          % (level, self._max_level), len(grids))
        # We have an N^2 check, so we try to be as quick as possible
        # and to skip as many as possible
        for pi, grid1 in enumerate(grids):
            pbar.update(pi)
            if grid1.retVal[0].shape[0] == 0: continue
            for grid2 in grid1._overlap_grids[self.axis]:
                if grid2.retVal[0].shape[0] == 0 \
                  or grid1.id == grid2.id:
                    continue
                args = grid1.retVal + grid2.retVal + [0]
                PointCombine.CombineData(*args)
                goodI = na.where(grid2.retVal[0] > -1)
                grid2.retVal = [grid2.retVal[i][goodI] for i in range(5)]
            numRefined = 0
            if level <= self._max_level and level > 0:
                for grid2 in grid1.Parent._overlap_grids[self.axis]:
                    if grid2.coarseData[0].shape[0] == 0: continue # Already refined
                    args = grid1.retVal[:3] + [grid1.retVal[4]] + \
                           grid2.coarseData[:3] + [grid2.coarseData[4]] + [2]
                    numRefined += PointCombine.RefineCoarseData(*args)
        pbar.finish()

    @time_execution
    def get_data(self, field = None):
        if not field: field = self.fields[0]
        all_data = []
        h = self.hierarchy
        for level in range(0, self._max_level+1):
            all_data.append(self.__project_level(level, field))
        all_data = na.concatenate(all_data, axis=1)
        # We now convert to half-widths and center-points
        self['pdx'] = all_data[3,:]
        self['px'] = (all_data[0,:]+0.5) * self['pdx']
        self['py'] = (all_data[1,:]+0.5) * self['pdx']
        self['pdx'] *= 0.5
        self['pdy'] = self['pdx'].copy()
        self.data[field] = all_data[2,:]
        # Now, we should clean up after ourselves...
        [grid.clear_all for grid in h.grids]

class Enzo3DData(EnzoData):
    """
    Class describing a cluster of data points, not necessarily sharing a
    coordinate.
    """
    _spatial = False
    _num_ghost_zones = 0
    def __init__(self, center, fields, pf = None):
        """
        Returns an instance of Enzo3DData, or prepares one.

        @param hierarchy: hierarchy we are associated with
        @type hierarchy: L{EnzoHierarchy<EnzoHierarchy>}
        @param center: center of the region
        @type center: array of floats
        @param fields: fields to read/generate
        @type fields: list of strings
        """
        EnzoData.__init__(self, pf, fields)
        self.center = center
        self.set_field_parameter("center",center)
        self.coords = None
        self.dx = None

    def write_out(self, fields, filename):
        """
        Doesn't work yet -- allPoints doesn't quite exist
        """
        if not isinstance(fields, types.ListType):
            fields = [fields]
        f = open(filename,"w")
        f.write("x\ty\tz\tdx\tdy\tdz")
        for field in fields:
            f.write("\t%s" % (field))
        f.write("\n")
        for i in range(self.x.shape[0]):
            f.write("%0.20f %0.20f %0.20f %0.20f %0.20f %0.20f" % \
                    (self.x[i], \
                     self.y[i], \
                     self.z[i], \
                     self.dx[i], \
                     self.dy[i], \
                     self.dz[i]))
            for field in fields:
                f.write("\t%0.7e" % (self[field][i]))
            f.write("\n")
        f.close()

    def make_profile(self, fields, num_bins, bounds, bin_field="RadiusCode",
                     take_log = True):
        """
        Here we make a profile.  Note that, for now, we do log-spacing in the
        bins, and we also assume we are given the radii in code units.
        field_weights defines the weighting for a given field.  If not given,
        assumed to be mass-weighted.

        Units might be screwy, but I don't think so.  Unless otherwise stated,
        returned in code units.  CellMass, by the way, is Msun.  And
        accumulated.

        Returns the bin outer radii and a dict of the profiles.

        @param fields: fields to get profiles of
        @type fields: list of strings
        @param nBins: number of bins
        @type nBins: int
        @param rInner: inner radius, code units
        @param rOuter: outer radius, code units
        """
        fields = ensure_list(fields)
        r_outer = min(bounds[1], self["Radius"].max())
        r_inner = min(bounds[0], self["Radius"].min())
        # Let's make the bins
        if logIt:
            bins = na.logspace(log10(r_inner), log10(r_outer), num_bins)
        else:
            bins = na.logspace(r_inner, r_outer, num_bins)
        radiiOrder = na.argsort(self[bin_field])
        fieldCopies = {} # We double up our memory usage here for sorting
        radii = self[binBy][radiiOrder]
        binIndices = na.searchsorted(bins, radii)
        nE = self[binBy].shape[0]
        defaultWeight = self["CellMass"][radiiOrder]
        fieldProfiles = {}
        if "CellMass" not in fields:
            fields.append("CellMass")
        for field in fields:
            code = WeaveStrings.ProfileBinningWeighted
            fc = self[field][radiiOrder]
            fp = na.zeros(nBins,'float64')
            if field_weights.has_key(field):
                if field_weights[field] == -999:
                    ss = "Accumulation weighting"
                    code = WeaveStrings.ProfileBinningAccumulation
                    weight = na.ones(nE, 'float64')
                elif field_weights[field] != None:
                    ww = field_weights[field]
                    ss="Weighting with %s" % (ww)
                    weight = self[ww][radiiOrder]
                elif field_weights[field] == None:
                    ss="Not weighted"
                    weight = na.ones(nE, 'float64')
                else:
                    mylog.warning("UNDEFINED weighting for %s; defaulting to unweighted", field)
                    ss="Undefined weighting"
                    weight = na.ones(nE, 'float64')
            else:
                ss="Weighting with default"
                weight = defaultWeight
            mylog.info("Getting profile of %s (%s)", field, ss)
            #print fc.dtype, binIndices.dtype, fp.dtype, weight.dtype
            #PointCombine.BinProfile(fc, binIndices, \
                                   #fp, weight)
            ld = { 'num_bins' : nBins,
                   'weightvalues' : weight,
                   'profilevalues' : fp,
                   'binindices' : binIndices,
                   'num_elements' : fc.shape[0],
                   'fieldvalues' : fc }
            weave.inline(code, ld.keys(), local_dict=ld, compiler='gcc',
                         type_converters=converters.blitz, auto_downcast = 0, verbose=2)
            fieldProfiles[field] = fp
        fieldProfiles[binBy] = bins[:nBins]
        co = AnalyzeClusterOutput(fieldProfiles)
        return co

    def _generate_coords(self):
        points = []
        for grid in self._grids:
            grid._generate_coords()
            grid.set_field_parameter("center", self.center)
            points.append(self._generate_grid_coords(grid))
        t = na.concatenate(points)
        self['x'] = t[:,0]
        self['y'] = t[:,1]
        self['z'] = t[:,2]
        self['RadiusCode'] = t[:,3]
        self['dx'] = t[:,4]
        self['dy'] = t[:,4]
        self['dz'] = t[:,4]

    def _generate_grid_coords(self, grid):
        pointI = self._get_point_indices(grid)
        dx = na.ones(pointI[0].shape[0], 'float64') * grid.dx
        tr = na.array([grid['x'][pointI].ravel(), \
                grid['y'][pointI].ravel(), \
                grid['z'][pointI].ravel(), \
                grid["RadiusCode"][pointI].ravel(),
                dx], 'float64').swapaxes(0,1)
        return tr

    def get_data(self, fields=None):
        self._get_list_of_grids()
        points = []
        if not self.has_key('dx'):
            self._generate_coords()
        if not fields:
            fields_to_get = self.fields
        else:
            fields_to_get = ensure_list(fields)
        for field in fields_to_get:
            if self.data.has_key(field):
                continue
            mylog.info("Getting field %s from %s", field, len(self._grids))
            if field not in self.hierarchy.fieldList:
                if self._generate_field(field):
                    continue # True means we already assigned it
            self[field] = na.concatenate(
                [self._get_data_from_grid(grid, field)
                 for grid in self._grids])

    def _get_data_from_grid(self, grid, field):
        pointI = self._get_point_indices(grid)
        return grid[field][pointI].ravel()

    def _generate_field(self, field):
        if fieldInfo.has_key(field):
            # First we check the validator
            try:
                fieldInfo[field].check_available(self)
            except NeedsGridType, ngt_exception:
                # We leave this to be implementation-specific
                self._generate_field_in_grids(field, ngt_exception.ghost_zones)
                return False
            else:
                self[field] = fieldInfo[field](self)
                return True
        else: # Can't find the field, try as it might
            raise exceptions.KeyError(field)

    def _generate_field_in_grids(self, field, num_ghost_zones=0):
        for grid in self._grids:
            temp = grid[field]

    def _get_point_indices(self, grid):
        pointI = na.where( self._get_cut_mask(grid) & grid.child_mask )
        return pointI

class EnzoRegionBase(Enzo3DData):
    def __init__(self, center, left_edge, right_edge, fields = None, pf = None):
        """
        Initializer for a rectangular prism of data.

        @note: Center does not have to be (rightEdge - leftEdge) / 2.0
        """
        Enzo3DData.__init__(self, center, fields, pf)
        self.fields = ["Radius"] + self.fields
        self.left_edge = left_edge
        self.right_edge = right_edge
        self._refresh_data()

    def _get_list_of_grids(self):
        self._grids, ind = self.pf.hierarchy.get_box_grids(self.left_edge, self.right_edge)

    def _get_cut_mask(self, grid):
        pointI = \
               ( (grid['x'] <= self.right_edge[0])
               & (grid['x'] >= self.left_edge[0])
               & (grid['y'] <= self.right_edge[1])
               & (grid['y'] >= self.left_edge[1])
               & (grid['z'] <= self.right_edge[2])
               & (grid['z'] >= self.left_edge[2]) )
        return pointI

class EnzoSphereBase(Enzo3DData):
    """
    A sphere of points
    """
    def __init__(self, center, radius, fields = None, pf = None):
        """
        @param hierarchy: hierarchy we are associated with
        @type hierarchy: L{EnzoHierarchy<EnzoHierarchy>}
        @param center: center of the region
        @type center: array of floats
        @param radius: radius of the sphere
        @type radius: float
        @keyword fields: fields to read/generate
        @type fields: list of strings
        """
        Enzo3DData.__init__(self, center, fields, pf)
        self.fields = ["Radius"] + self.fields
        self.radius = radius
        self._refresh_data()

    def _get_list_of_grids(self, field = None):
        self._grids,ind = self.hierarchy.find_sphere_grids(self.center, self.radius)

    def _get_cut_mask(self, grid):
        pointI = ((grid["RadiusCode"]<=self.radius) &
                  (grid.child_mask==1))
        return pointI

class EnzoCoveringGrid(Enzo3DData):
    def __init__(self, level, left_edge, right_edge, dims, fields = None,
                 pf = None):
        """
        """

        Enzo3DData.__init__(self, center=None, fields=fields, pf=pf)
        self.left_edge = na.array(left_edge)
        self.right_edge = na.array(right_edge)
        self.level = level
        self.ActiveDimensions = na.array(dims)
        self.dx, self.dy, self.dz = (self.right_edge-self.left_edge) \
                                  / self.ActiveDimensions
        self._refresh_data()

    def _get_list_of_grids(self):
        grids, ind = self.pf.hierarchy.get_box_grids(self.left_edge, self.right_edge)
        ind = na.where(self.pf.hierarchy.grids[ind] <= self.level)
        self._grids = self.pf.hierarchy.grids[ind]

    def __setup_weave_dict(self, grid):
        return {
            'nxc': int(grid.ActiveDimensions[0]),
            'nyc': int(grid.ActiveDimensions[1]),
            'nzc': int(grid.ActiveDimensions[2]),
            'leftEdgeCoarse': na.array(grid.LeftEdge),
            'rf': int(grid.dx / self.dx),
            'dx_c': float(grid.dx),
            'dy_c': float(grid.dy),
            'dz_c': float(grid.dz),
            'dx_f': float(self.dx),
            'dy_f': float(self.dy),
            'dz_f': float(self.dz),
            'cubeLeftEdge': na.array(self.left_edge),
            'cubeRightEdge': na.array(self.right_edge),
            'nxf': int(self.ActiveDimensions[0]),
            'nyf': int(self.ActiveDimensions[1]),
            'nzf': int(self.ActiveDimensions[2]),
            'childMask' : grid.child_mask
        }

    def get_data(self, field=None):
        self._get_list_of_grids()
        # We don't generate coordinates here.
        if field == None:
            fields_to_get = self.fields
        else:
            fields_to_get = ensure_list(field)
        for field in fields_to_get:
            if self.data.has_key(field):
                continue
            mylog.info("Getting field %s from %s", field, len(self._grids))
            self[field] = na.zeros(self.ActiveDimensions, dtype='float32')
            for grid in self._grids:
                self._get_data_from_grid(grid, field)

    def _get_data_from_grid(self, grid, fields):
        last_level = (grid.dx / 2.0) < self.dx
        for field in ensure_list(fields):
            locals_dict = self.__setup_weave_dict(grid)
            locals_dict['fieldData'] = grid[field]
            locals_dict['cubeData'] = self[field]
            locals_dict['lastLevel'] = int(last_level)
            weave.inline(DataCubeRefineCoarseData,
                         locals_dict.keys(), local_dict=locals_dict,
                         compiler='gcc',
                         type_converters=converters.blitz,
                         auto_downcast=0, verbose=2)
