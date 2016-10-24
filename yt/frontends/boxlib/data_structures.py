"""
Data structures for BoxLib Codes



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import inspect
import os
import re

from stat import ST_CTIME

import numpy as np

from yt.funcs import \
    ensure_tuple, \
    mylog, \
    setdefaultattr
from yt.data_objects.grid_patch import AMRGridPatch
from yt.extern.six.moves import zip as izip
from yt.geometry.grid_geometry_handler import GridIndex
from yt.data_objects.static_output import Dataset

from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_root_only
from yt.utilities.lib.misc_utilities import \
    get_box_grids_level
from yt.utilities.io_handler import \
    io_registry

from .fields import \
    BoxlibFieldInfo, \
    MaestroFieldInfo, \
    CastroFieldInfo


# This is what we use to find scientific notation that might include d's
# instead of e's.
_scinot_finder = re.compile(r"[-+]?[0-9]*\.?[0-9]+([eEdD][-+]?[0-9]+)?")
# This is the dimensions in the Cell_H file for each level
# It is different for different dimensionalities, so we make a list
_dim_finder = [
    re.compile(r"\(\((\d+)\) \((\d+)\) \(\d+\)\)$"),
    re.compile(r"\(\((\d+,\d+)\) \((\d+,\d+)\) \(\d+,\d+\)\)$"),
    re.compile(r"\(\((\d+,\d+,\d+)\) \((\d+,\d+,\d+)\) \(\d+,\d+,\d+\)\)$")]
# This is the line that prefixes each set of data for a FAB in the FAB file
# It is different for different dimensionalities, so we make a list
_endian_regex = r"^FAB \(\(\d+, \([0-9 ]+\)\),\((\d+), \(([0-9 ]+)\)\)\)"
_header_pattern = [
    re.compile(_endian_regex +
               r"\(\((\d+)\) \((\d+)\) \((\d+)\)\) (\d+)\n"),
    re.compile(_endian_regex +
               r"\(\((\d+,\d+)\) \((\d+,\d+)\) \((\d+,\d+)\)\) (\d+)\n"),
    re.compile(_endian_regex +
               r"\(\((\d+,\d+,\d+)\) \((\d+,\d+,\d+)\) \((\d+,\d+,\d+)\)\) (\d+)\n")]


class BoxlibGrid(AMRGridPatch):
    _id_offset = 0
    _offset = -1

    def __init__(self, grid_id, offset, filename=None,
                 index=None):
        super(BoxlibGrid, self).__init__(grid_id, filename, index)
        self._base_offset = offset
        self._parent_id = []
        self._children_ids = []

    def _prepare_grid(self):
        super(BoxlibGrid, self)._prepare_grid()
        my_ind = self.id - self._id_offset
        self.start_index = self.index.grid_start_index[my_ind]

    def get_global_startindex(self):
        return self.start_index

    def _setup_dx(self):
        # has already been read in and stored in index
        self.dds = self.index.ds.arr(self.index.level_dds[self.Level, :], 'code_length')
        self.field_data['dx'], self.field_data['dy'], self.field_data['dz'] = self.dds

    def __repr__(self):
        return "BoxlibGrid_%04i" % (self.id)

    @property
    def Parent(self):
        if len(self._parent_id) == 0:
            return None
        return [self.index.grids[pid - self._id_offset]
                for pid in self._parent_id]

    @property
    def Children(self):
        return [self.index.grids[cid - self._id_offset]
                for cid in self._children_ids]

    def _get_offset(self, f):
        # This will either seek to the _offset or figure out the correct
        # _offset.
        if self._offset == -1:
            f.seek(self._base_offset, os.SEEK_SET)
            f.readline()
            self._offset = f.tell()
        return self._offset

    # We override here because we can have varying refinement levels
    def select_ires(self, dobj):
        mask = self._get_selector_mask(dobj.selector)
        if mask is None: return np.empty(0, dtype='int64')
        coords = np.empty(self._last_count, dtype='int64')
        coords[:] = self.Level + self.ds.level_offsets[self.Level]
        return coords

    # Override this as well, since refine_by can vary
    def _fill_child_mask(self, child, mask, tofill, dlevel=1):
        rf = self.ds.ref_factors[self.Level]
        if dlevel != 1:
            raise NotImplementedError
        gi, cgi = self.get_global_startindex(), child.get_global_startindex()
        startIndex = np.maximum(0, cgi // rf - gi)
        endIndex = np.minimum((cgi + child.ActiveDimensions) // rf - gi,
                              self.ActiveDimensions)
        endIndex += (startIndex == endIndex)
        mask[startIndex[0]:endIndex[0],
             startIndex[1]:endIndex[1],
             startIndex[2]:endIndex[2]] = tofill


class BoxlibHierarchy(GridIndex):
    grid = BoxlibGrid

    def __init__(self, ds, dataset_type='boxlib_native'):
        self.dataset_type = dataset_type
        self.header_filename = os.path.join(ds.output_dir, 'Header')
        self.directory = ds.output_dir

        GridIndex.__init__(self, ds, dataset_type)
        self._cache_endianness(self.grids[-1])

    def _parse_index(self):
        """
        read the global header file for an Boxlib plotfile output.
        """
        self.max_level = self.dataset._max_level
        header_file = open(self.header_filename, 'r')

        self.dimensionality = self.dataset.dimensionality
        _our_dim_finder = _dim_finder[self.dimensionality-1]
        DRE = self.dataset.domain_right_edge  # shortcut
        DLE = self.dataset.domain_left_edge   # shortcut

        # We can now skip to the point in the file we want to start parsing.
        header_file.seek(self.dataset._header_mesh_start)

        dx = []
        for i in range(self.max_level + 1):
            dx.append([float(v) for v in next(header_file).split()])
            # account for non-3d data sets
            if self.dimensionality < 2:
                dx[i].append(DRE[1] - DLE[1])
            if self.dimensionality < 3:
                dx[i].append(DRE[2] - DLE[1])
        self.level_dds = np.array(dx, dtype="float64")
        next(header_file)
        if self.ds.geometry == "cartesian":
            default_ybounds = (0.0, 1.0)
            default_zbounds = (0.0, 1.0)
        elif self.ds.geometry == "cylindrical":
            # Now we check for dimensionality issues
            if self.dimensionality != 2:
                raise RuntimeError("yt needs cylindrical to be 2D")
            self.level_dds[:,2] = 2*np.pi
            default_zbounds = (0.0, 2*np.pi)
        elif self.ds.geometry == "spherical":
            # BoxLib only supports 1D spherical, so ensure
            # the other dimensions have the right extent.
            self.level_dds[:,1] = np.pi
            self.level_dds[:,2] = 2*np.pi
            default_ybounds = (0.0, np.pi)
            default_zbounds = (0.0, 2*np.pi)
        else:
            raise RuntimeError("Unknown BoxLib coordinate system.")
        if int(next(header_file)) != 0:
            raise RuntimeError("INTERNAL ERROR! This should be a zero.")

        # each level is one group with ngrids on it.
        # each grid has self.dimensionality number of lines of 2 reals
        self.grids = []
        grid_counter = 0
        for level in range(self.max_level + 1):
            vals = next(header_file).split()
            lev, ngrids = int(vals[0]), int(vals[1])
            assert(lev == level)
            nsteps = int(next(header_file))  # NOQA
            for gi in range(ngrids):
                xlo, xhi = [float(v) for v in next(header_file).split()]
                if self.dimensionality > 1:
                    ylo, yhi = [float(v) for v in next(header_file).split()]
                else:
                    ylo, yhi = default_ybounds
                if self.dimensionality > 2:
                    zlo, zhi = [float(v) for v in next(header_file).split()]
                else:
                    zlo, zhi = default_zbounds
                self.grid_left_edge[grid_counter + gi, :] = [xlo, ylo, zlo]
                self.grid_right_edge[grid_counter + gi, :] = [xhi, yhi, zhi]
            # Now we get to the level header filename, which we open and parse.
            fn = os.path.join(self.dataset.output_dir,
                              next(header_file).strip())
            level_header_file = open(fn + "_H")
            level_dir = os.path.dirname(fn)
            # We skip the first two lines, which contain BoxLib header file
            # version and 'how' the data was written
            next(level_header_file)
            next(level_header_file)
            # Now we get the number of components
            ncomp_this_file = int(next(level_header_file))  # NOQA
            # Skip the next line, which contains the number of ghost zones
            next(level_header_file)
            # To decipher this next line, we expect something like:
            # (8 0
            # where the first is the number of FABs in this level.
            ngrids = int(next(level_header_file).split()[0][1:])
            # Now we can iterate over each and get the indices.
            for gi in range(ngrids):
                # components within it
                start, stop = _our_dim_finder.match(next(level_header_file)).groups()
                # fix for non-3d data 
                # note we append '0' to both ends b/c of the '+1' in dims below
                start += ',0'*(3-self.dimensionality)
                stop += ',0'*(3-self.dimensionality)
                start = np.array(start.split(","), dtype="int64")
                stop = np.array(stop.split(","), dtype="int64")
                dims = stop - start + 1
                self.grid_dimensions[grid_counter + gi,:] = dims
                self.grid_start_index[grid_counter + gi,:] = start
            # Now we read two more lines.  The first of these is a close
            # parenthesis.
            next(level_header_file)
            # The next is again the number of grids
            next(level_header_file)
            # Now we iterate over grids to find their offsets in each file.
            for gi in range(ngrids):
                # Now we get the data file, at which point we're ready to
                # create the grid.
                dummy, filename, offset = next(level_header_file).split()
                filename = os.path.join(level_dir, filename)
                go = self.grid(grid_counter + gi, int(offset), filename, self)
                go.Level = self.grid_levels[grid_counter + gi,:] = level
                self.grids.append(go)
            grid_counter += ngrids
            # already read the filenames above...
        self.float_type = 'float64'

    def _cache_endianness(self, test_grid):
        """
        Cache the endianness and bytes perreal of the grids by using a
        test grid and assuming that all grids have the same
        endianness. This is a pretty safe assumption since Boxlib uses
        one file per processor, and if you're running on a cluster
        with different endian processors, then you're on your own!
        """
        # open the test file & grab the header
        with open(os.path.expanduser(test_grid.filename), 'rb') as f:
            header = f.readline().decode("ascii", "ignore")

        bpr, endian, start, stop, centering, nc = \
            _header_pattern[self.dimensionality-1].search(header).groups()
        # Note that previously we were using a different value for BPR than we
        # use now.  Here is an example set of information directly from BoxLib:
        #  * DOUBLE data
        #  * FAB ((8, (64 11 52 0 1 12 0 1023)),(8, (1 2 3 4 5 6 7 8)))((0,0) (63,63) (0,0)) 27
        #  * FLOAT data
        #  * FAB ((8, (32 8 23 0 1 9 0 127)),(4, (1 2 3 4)))((0,0) (63,63) (0,0)) 27
        if bpr == endian[0]:
            dtype = '<f%s' % bpr
        elif bpr == endian[-1]:
            dtype = '>f%s' % bpr
        else:
            raise ValueError("FAB header is neither big nor little endian. Perhaps the file is corrupt?")

        mylog.debug("FAB header suggests dtype of %s", dtype)
        self._dtype = np.dtype(dtype)

    def _populate_grid_objects(self):
        mylog.debug("Creating grid objects")
        self.grids = np.array(self.grids, dtype='object')
        self._reconstruct_parent_child()
        for i, grid in enumerate(self.grids):
            if (i % 1e4) == 0: mylog.debug("Prepared % 7i / % 7i grids", i,
                                           self.num_grids)
            grid._prepare_grid()
            grid._setup_dx()
        mylog.debug("Done creating grid objects")

    def _reconstruct_parent_child(self):
        mask = np.empty(len(self.grids), dtype='int32')
        mylog.debug("First pass; identifying child grids")
        for i, grid in enumerate(self.grids):
            get_box_grids_level(self.grid_left_edge[i,:],
                                self.grid_right_edge[i,:],
                                self.grid_levels[i] + 1,
                                self.grid_left_edge, self.grid_right_edge,
                                self.grid_levels, mask)
            ids = np.where(mask.astype("bool"))  # where is a tuple
            grid._children_ids = ids[0] + grid._id_offset
        mylog.debug("Second pass; identifying parents")
        for i, grid in enumerate(self.grids):  # Second pass
            for child in grid.Children:
                child._parent_id.append(i + grid._id_offset)

    def _count_grids(self):
        # We can get everything from the Header file, but note that we're
        # duplicating some work done elsewhere.  In a future where we don't
        # pre-allocate grid arrays, this becomes unnecessary.
        header_file = open(self.header_filename, 'r')
        header_file.seek(self.dataset._header_mesh_start)
        # Skip over the level dxs, geometry and the zero:
        [next(header_file) for i in range(self.dataset._max_level + 3)]
        # Now we need to be very careful, as we've seeked, and now we iterate.
        # Does this work?  We are going to count the number of places that we
        # have a three-item line.  The three items would be level, number of
        # grids, and then grid time.
        self.num_grids = 0
        for line in header_file:
            if len(line.split()) != 3: continue
            self.num_grids += int(line.split()[1])

    def _initialize_grid_arrays(self):
        super(BoxlibHierarchy, self)._initialize_grid_arrays()
        self.grid_start_index = np.zeros((self.num_grids, 3), 'int64')

    def _initialize_state_variables(self):
        """override to not re-initialize num_grids in AMRHierarchy.__init__

        """
        self._parallel_locking = False
        self._data_file = None
        self._data_mode = None

    def _detect_output_fields(self):
        # This is all done in _parse_header_file
        self.field_list = [("boxlib", f) for f in
                           self.dataset._field_list]
        self.field_indexes = dict((f[1], i)
                                  for i, f in enumerate(self.field_list))
        # There are times when field_list may change.  We copy it here to
        # avoid that possibility.
        self.field_order = [f for f in self.field_list]

    def _setup_data_io(self):
        self.io = io_registry[self.dataset_type](self.dataset)


class BoxlibDataset(Dataset):
    """
    This class is a stripped down class that simply reads and parses
    *filename*, without looking at the Boxlib index.
    """
    _index_class = BoxlibHierarchy
    _field_info_class = BoxlibFieldInfo
    _output_prefix = None

    # THIS SHOULD BE FIXED:
    periodicity = (True, True, True)

    def __init__(self, output_dir,
                 cparam_filename="inputs",
                 fparam_filename="probin",
                 dataset_type='boxlib_native',
                 storage_filename=None,
                 units_override=None,
                 unit_system="cgs"):
        """
        The paramfile is usually called "inputs"
        and there may be a fortran inputs file usually called "probin"
        plotname here will be a directory name
        as per BoxLib, dataset_type will be Native (implemented here), IEEE (not
        yet implemented) or ASCII (not yet implemented.)
        """
        self.fluid_types += ("boxlib",)
        self.output_dir = os.path.abspath(os.path.expanduser(output_dir))
        self.cparam_filename = self._localize_check(cparam_filename)
        self.fparam_filename = self._localize_check(fparam_filename)
        self.storage_filename = storage_filename

        Dataset.__init__(self, output_dir, dataset_type,
                         units_override=units_override,
                         unit_system=unit_system)

        # These are still used in a few places.
        if "HydroMethod" not in self.parameters.keys():
            self.parameters["HydroMethod"] = 'boxlib'
        self.parameters["Time"] = 1.     # default unit is 1...
        self.parameters["EOSType"] = -1  # default
        self.parameters["gamma"] = self.parameters.get(
            "materials.gamma", 1.6667)

    def _localize_check(self, fn):
        # If the file exists, use it.  If not, set it to None.
        root_dir = os.path.dirname(self.output_dir)
        full_fn = os.path.join(root_dir, fn)
        if os.path.exists(full_fn):
            return full_fn
        return None

    @classmethod
    def _is_valid(cls, *args, **kwargs):
        # fill our args
        output_dir = args[0]
        # boxlib datasets are always directories
        if not os.path.isdir(output_dir): return False
        header_filename = os.path.join(output_dir, "Header")
        jobinfo_filename = os.path.join(output_dir, "job_info")
        if not os.path.exists(header_filename):
            # We *know* it's not boxlib if Header doesn't exist.
            return False
        args = inspect.getcallargs(cls.__init__, args, kwargs)
        # This might need to be localized somehow
        inputs_filename = os.path.join(
            os.path.dirname(os.path.abspath(output_dir)),
            args['cparam_filename'])
        if not os.path.exists(inputs_filename) and \
           not os.path.exists(jobinfo_filename):
            return True  # We have no parameters to go off of
        # If we do have either inputs or jobinfo, we should be deferring to a
        # different frontend.
        return False

    def _parse_parameter_file(self):
        """
        Parses the parameter file and establishes the various
        dictionaries.
        """
        self._parse_header_file()
        # Let's read the file
        hfn = os.path.join(self.output_dir, 'Header')
        self.unique_identifier = int(os.stat(hfn)[ST_CTIME])
        # the 'inputs' file is now optional
        self._parse_cparams()
        self._parse_fparams()

    def _parse_cparams(self):
        if self.cparam_filename is None:
            return
        for line in (line.split("#")[0].strip() for line in
                     open(self.cparam_filename)):
            if "=" not in line: continue
            if len(line) == 0: continue
            param, vals = [s.strip() for s in line.split("=")]
            if param == "amr.n_cell":
                vals = self.domain_dimensions = np.array(vals.split(), dtype='int32')

                # For 1D and 2D simulations in BoxLib usually only the relevant dimensions
                # have a specified number of zones, but yt requires domain_dimensions to 
                # have three elements, with 1 in the additional slots if we're not in 3D, 
                # so append them as necessary.

                if (len(vals) == 1):
                    vals = self.domain_dimensions = np.array([vals[0], 1, 1])
                elif (len(vals) == 2):
                    vals = self.domain_dimensions = np.array([vals[0], vals[1], 1])
            elif param == "amr.ref_ratio":
                vals = self.refine_by = int(vals[0])
            elif param == "Prob.lo_bc":
                vals = self.periodicity = ensure_tuple([i == 0 for i in vals.split()])
            elif param == "castro.use_comoving":
                vals = self.cosmological_simulation = int(vals)
            else:
                vals = _guess_pcast(vals)
            self.parameters[param] = vals

        if getattr(self, "cosmological_simulation", 0) == 1:
            self.omega_lambda = self.parameters["comoving_OmL"]
            self.omega_matter = self.parameters["comoving_OmM"]
            self.hubble_constant = self.parameters["comoving_h"]
            a_file = open(os.path.join(self.output_dir, 'comoving_a'))
            line = a_file.readline().strip()
            a_file.close()
            self.current_redshift = 1/float(line) - 1
        else:
            self.current_redshift = self.omega_lambda = self.omega_matter = \
                self.hubble_constant = self.cosmological_simulation = 0.0

    def _parse_fparams(self):
        """
        Parses the fortran parameter file for Orion. Most of this will
        be useless, but this is where it keeps mu = mass per
        particle/m_hydrogen.
        """
        if self.fparam_filename is None:
            return
        for line in (l for l in open(self.fparam_filename) if "=" in l):
            param, vals = [v.strip() for v in line.split("=")]
            # Now, there are a couple different types of parameters.
            # Some will be where you only have floating point values, others
            # will be where things are specified as string literals.
            # Unfortunately, we're also using Fortran values, which will have
            # things like 1.d-2 which is pathologically difficult to parse if
            # your C library doesn't include 'd' in its locale for strtod.
            # So we'll try to determine this.
            vals = vals.split()
            if any(_scinot_finder.match(v) for v in vals):
                vals = [float(v.replace("D", "e").replace("d", "e"))
                        for v in vals]
            if len(vals) == 1:
                vals = vals[0]
            self.parameters[param] = vals

    def _parse_header_file(self):
        """
        We parse the Boxlib header, which we use as our basis.  Anything in the
        inputs file will override this, but the inputs file is not strictly
        necessary for orientation of the data in space.
        """

        # Note: Python uses a read-ahead buffer, so using next(), which would
        # be my preferred solution, won't work here.  We have to explicitly
        # call readline() if we want to end up with an offset at the very end.
        # Fortunately, elsewhere we don't care about the offset, so we're fine
        # everywhere else using iteration exclusively.
        header_file = open(os.path.join(self.output_dir, 'Header'))
        self.orion_version = header_file.readline().rstrip()
        n_fields = int(header_file.readline())

        self._field_list = [header_file.readline().strip()
                            for i in range(n_fields)]

        self.dimensionality = int(header_file.readline())
        self.current_time = float(header_file.readline())
        # This is traditionally a index attribute, so we will set it, but
        # in a slightly hidden variable.
        self._max_level = int(header_file.readline())
        self.domain_left_edge = np.array(header_file.readline().split(),
                                         dtype="float64")
        self.domain_right_edge = np.array(header_file.readline().split(),
                                          dtype="float64")
        ref_factors = np.array([int(i) for i in
                                header_file.readline().split()])
        if ref_factors.size == 0:
            # We use a default of two, as Nyx doesn't always output this value
            ref_factors = [2] * (self._max_level + 1)
        # We can't vary refinement factors based on dimension, or whatever else
        # they are vaied on.  In one curious thing, I found that some Castro 3D
        # data has only two refinement factors, which I don't know how to
        # understand.
        self.ref_factors = ref_factors
        if np.unique(ref_factors).size > 1:
            # We want everything to be a multiple of this.
            self.refine_by = min(ref_factors)
            # Check that they're all multiples of the minimum.
            if not all(float(rf)/self.refine_by ==
                       int(float(rf)/self.refine_by) for rf in ref_factors):
                raise RuntimeError
            base_log = np.log2(self.refine_by)
            self.level_offsets = [0]  # level 0 has to have 0 offset
            lo = 0
            for lm1, rf in enumerate(self.ref_factors):
                lo += int(np.log2(rf) / base_log) - 1
                self.level_offsets.append(lo)
        # assert(np.unique(ref_factors).size == 1)
        else:
            self.refine_by = ref_factors[0]
            self.level_offsets = [0 for l in range(self._max_level + 1)]
        # Now we read the global index space, to get
        index_space = header_file.readline()
        # This will be of the form:
        #  ((0,0,0) (255,255,255) (0,0,0)) ((0,0,0) (511,511,511) (0,0,0))
        # So note that if we split it all up based on spaces, we should be
        # fine, as long as we take the first two entries, which correspond to
        # the root level.  I'm not 100% pleased with this solution.
        root_space = index_space.replace("(", "").replace(")", "").split()[:2]
        start = np.array(root_space[0].split(","), dtype="int64")
        stop = np.array(root_space[1].split(","), dtype="int64")
        self.domain_dimensions = stop - start + 1
        # Skip timesteps per level
        header_file.readline()
        self._header_mesh_start = header_file.tell()
        # Skip the cell size information per level - we'll get this later
        for i in range(self._max_level+1): header_file.readline()
        # Get the geometry
        next_line = header_file.readline()
        if len(next_line.split()) == 1:
            coordinate_type = int(next_line)
        else:
            coordinate_type = 0
        if coordinate_type == 0:
            self.geometry = "cartesian"
        elif coordinate_type == 1:
            self.geometry = "cylindrical"
        elif coordinate_type == 2:
            self.geometry = "spherical"
        else:
            raise RuntimeError("Unknown BoxLib coord_type")

        # overrides for 1/2-dimensional data
        if self.dimensionality == 1:
            self._setup1d()
        elif self.dimensionality == 2:
            self._setup2d()

    def _set_code_unit_attributes(self):
        setdefaultattr(self, 'length_unit', self.quan(1.0, "cm"))
        setdefaultattr(self, 'mass_unit', self.quan(1.0, "g"))
        setdefaultattr(self, 'time_unit', self.quan(1.0, "s"))
        setdefaultattr(self, 'velocity_unit', self.quan(1.0, "cm/s"))

    def _setup1d(self):
        # self._index_class = BoxlibHierarchy1D
        # self._fieldinfo_fallback = Orion1DFieldInfo
        self.domain_left_edge = \
            np.concatenate([self.domain_left_edge, [0.0, 0.0]])
        self.domain_right_edge = \
            np.concatenate([self.domain_right_edge, [1.0, 1.0]])
        tmp = self.domain_dimensions.tolist()
        tmp.extend((1, 1))
        self.domain_dimensions = np.array(tmp)
        tmp = list(self.periodicity)
        tmp[1] = False
        tmp[2] = False
        self.periodicity = ensure_tuple(tmp)

    def _setup2d(self):
        self.domain_left_edge = \
            np.concatenate([self.domain_left_edge, [0.0]])
        self.domain_right_edge = \
            np.concatenate([self.domain_right_edge, [1.0]])
        if self.geometry == "cylindrical":
            self.domain_right_edge[2] = 2.0 * np.pi
        tmp = self.domain_dimensions.tolist()
        tmp.append(1)
        self.domain_dimensions = np.array(tmp)
        tmp = list(self.periodicity)
        tmp[2] = False
        self.periodicity = ensure_tuple(tmp)

    @parallel_root_only
    def print_key_parameters(self):
        for a in ["current_time", "domain_dimensions", "domain_left_edge",
                  "domain_right_edge"]:
            if not hasattr(self, a):
                mylog.error("Missing %s in parameter file definition!", a)
                continue
            v = getattr(self, a)
            mylog.info("Parameters: %-25s = %s", a, v)

    def relative_refinement(self, l0, l1):
        offset = self.level_offsets[l1] - self.level_offsets[l0]
        return self.refine_by**(l1-l0 + offset)


class OrionHierarchy(BoxlibHierarchy):

    def __init__(self, ds, dataset_type='orion_native'):
        BoxlibHierarchy.__init__(self, ds, dataset_type)
        self._read_particles()
        # self.io = IOHandlerOrion

    def _detect_output_fields(self):
        # This is all done in _parse_header_file
        self.field_list = [("boxlib", f) for f in
                           self.dataset._field_list]
        self.field_indexes = dict((f[1], i)
                                  for i, f in enumerate(self.field_list))
        # There are times when field_list may change.  We copy it here to
        # avoid that possibility.
        self.field_order = [f for f in self.field_list]

        # look for particle fields
        self.particle_filename = None
        for particle_filename in ["StarParticles", "SinkParticles"]:
            fn = os.path.join(self.ds.output_dir, particle_filename)
            if os.path.exists(fn):
                self.particle_filename = fn

        if self.particle_filename is None:
            return

        pfield_list = [("io", c) for c in self.io.particle_field_index.keys()]
        self.field_list.extend(pfield_list)

    def _read_particles(self):
        """
        reads in particles and assigns them to grids. Will search for
        Star particles, then sink particles if no star particle file
        is found, and finally will simply note that no particles are
        found if neither works. To add a new Orion particle type,
        simply add it to the if/elif/else block.

        """
        self.grid_particle_count = np.zeros(len(self.grids))

        if self.particle_filename is not None:
            self._read_particle_file(self.particle_filename)

    def _read_particle_file(self, fn):
        """actually reads the orion particle data file itself.

        """
        if not os.path.exists(fn): return
        with open(fn, 'r') as f:
            lines = f.readlines()
            self.num_stars = int(lines[0].strip()[0])
            for num, line in enumerate(lines[1:]):
                particle_position_x = float(line.split(' ')[1])
                particle_position_y = float(line.split(' ')[2])
                particle_position_z = float(line.split(' ')[3])
                coord = [particle_position_x, particle_position_y, particle_position_z]
                # for each particle, determine which grids contain it
                # copied from object_finding_mixin.py
                mask = np.ones(self.num_grids)
                for i in range(len(coord)):
                    np.choose(np.greater(self.grid_left_edge.d[:,i],coord[i]), (mask,0), mask)
                    np.choose(np.greater(self.grid_right_edge.d[:,i],coord[i]), (0,mask), mask)
                ind = np.where(mask == 1)
                selected_grids = self.grids[ind]
                # in orion, particles always live on the finest level.
                # so, we want to assign the particle to the finest of
                # the grids we just found
                if len(selected_grids) != 0:
                    grid = sorted(selected_grids, key=lambda grid: grid.Level)[-1]
                    ind = np.where(self.grids == grid)[0][0]
                    self.grid_particle_count[ind] += 1
                    self.grids[ind].NumberOfParticles += 1

                    # store the position in the particle file for fast access.
                    try:
                        self.grids[ind]._particle_line_numbers.append(num + 1)
                    except AttributeError:
                        self.grids[ind]._particle_line_numbers = [num + 1]
        return True


class OrionDataset(BoxlibDataset):

    _index_class = OrionHierarchy

    def __init__(self, output_dir,
                 cparam_filename="inputs",
                 fparam_filename="probin",
                 dataset_type='orion_native',
                 storage_filename=None,
                 units_override=None,
                 unit_system="cgs"):

        BoxlibDataset.__init__(self, output_dir,
                               cparam_filename, fparam_filename,
                               dataset_type, units_override=units_override,
                               unit_system=unit_system)

    @classmethod
    def _is_valid(cls, *args, **kwargs):
        # fill our args
        output_dir = args[0]
        # boxlib datasets are always directories
        if not os.path.isdir(output_dir): return False
        header_filename = os.path.join(output_dir, "Header")
        jobinfo_filename = os.path.join(output_dir, "job_info")
        if not os.path.exists(header_filename):
            # We *know* it's not boxlib if Header doesn't exist.
            return False
        args = inspect.getcallargs(cls.__init__, args, kwargs)
        # This might need to be localized somehow
        inputs_filename = os.path.join(
            os.path.dirname(os.path.abspath(output_dir)),
            args['cparam_filename'])
        if not os.path.exists(inputs_filename):
            return False
        if os.path.exists(jobinfo_filename):
            return False
        # Now we check for all the others
        lines = open(inputs_filename).readlines()
        if any(("castro." in line for line in lines)): return False
        if any(("nyx." in line for line in lines)): return False
        if any(("maestro" in line.lower() for line in lines)): return False
        if any(("geometry.prob_lo" in line for line in lines)): return True
        return False


class CastroDataset(BoxlibDataset):

    _field_info_class = CastroFieldInfo

    @classmethod
    def _is_valid(cls, *args, **kwargs):
        # fill our args
        output_dir = args[0]
        # boxlib datasets are always directories
        if not os.path.isdir(output_dir): return False
        header_filename = os.path.join(output_dir, "Header")
        jobinfo_filename = os.path.join(output_dir, "job_info")
        if not os.path.exists(header_filename):
            # We *know* it's not boxlib if Header doesn't exist.
            return False
        if not os.path.exists(jobinfo_filename):
            return False
        # Now we check for all the others
        lines = open(jobinfo_filename).readlines()
        if any(line.startswith("Castro   ") for line in lines): return True
        return False

    def _parse_parameter_file(self):
        super(CastroDataset, self)._parse_parameter_file()
        jobinfo_filename = os.path.join(self.output_dir, "job_info")
        line = ""
        with open(jobinfo_filename, "r") as f:
            while not line.startswith(" Inputs File Parameters"):
                # boundary condition info starts with -x:, etc.
                bcs = ["-x:", "+x:", "-y:", "+y:", "-z:", "+z:"]
                if any(b in line for b in bcs):
                    p, v = line.strip().split(":")
                    self.parameters[p] = v.strip()
                if "git hash" in line:
                    # line format: codename git hash:  the-hash
                    fields = line.split(":")
                    self.parameters[fields[0]] = fields[1].strip()
                line = next(f)
            
            # runtime parameters that we overrode follow "Inputs File
            # Parameters"
            # skip the "====..." line
            line = next(f)
            for line in f:
                p, v = line.strip().split("=")
                self.parameters[p] = v.strip()

            
        # hydro method is set by the base class -- override it here
        self.parameters["HydroMethod"] = "Castro"

        # set the periodicity based on the runtime parameters
        periodicity = [True, True, True]
        if not self.parameters['-x'] == "interior": periodicity[0] = False
        if self.dimensionality >= 2:
            if not self.parameters['-y'] == "interior": periodicity[1] = False
        if self.dimensionality == 3:
            if not self.parameters['-z'] == "interior": periodicity[2] = False

        self.periodicity = ensure_tuple(periodicity)
    

class MaestroDataset(BoxlibDataset):

    _field_info_class = MaestroFieldInfo

    @classmethod
    def _is_valid(cls, *args, **kwargs):
        # fill our args
        output_dir = args[0]
        # boxlib datasets are always directories
        if not os.path.isdir(output_dir): return False
        header_filename = os.path.join(output_dir, "Header")
        jobinfo_filename = os.path.join(output_dir, "job_info")
        if not os.path.exists(header_filename):
            # We *know* it's not boxlib if Header doesn't exist.
            return False
        if not os.path.exists(jobinfo_filename):
            return False
        # Now we check the job_info for the mention of maestro
        lines = open(jobinfo_filename).readlines()
        if any(line.startswith("MAESTRO   ") for line in lines): return True
        return False

    def _parse_parameter_file(self):
        super(MaestroDataset, self)._parse_parameter_file()
        jobinfo_filename = os.path.join(self.output_dir, "job_info")
        line = ""
        with open(jobinfo_filename, "r") as f:
            while not line.startswith(" [*] indicates overridden default"):
                # get the code git hashes
                if "git hash" in line:
                    # line format: codename git hash:  the-hash
                    fields = line.split(":")
                    self.parameters[fields[0]] = fields[1].strip()
                line = next(f)
            # get the runtime parameters
            for line in f:
                p, v = (_.strip() for _ in line[4:].split("=", 1))
                if len(v) == 0:
                    self.parameters[p] = ""
                else:
                    self.parameters[p] = _guess_pcast(v)
            # hydro method is set by the base class -- override it here
            self.parameters["HydroMethod"] = "Maestro"

        # set the periodicity based on the integer BC runtime parameters
        periodicity = [True, True, True]
        if not self.parameters['bcx_lo'] == -1:
            periodicity[0] = False

        if not self.parameters['bcy_lo'] == -1:
            periodicity[1] = False

        if not self.parameters['bcz_lo'] == -1:
            periodicity[2] = False

        self.periodicity = ensure_tuple(periodicity)


class NyxHierarchy(BoxlibHierarchy):

    def __init__(self, ds, dataset_type='nyx_native'):
        super(NyxHierarchy, self).__init__(ds, dataset_type)
        self._read_particle_header()

    def _read_particle_header(self):
        if not self.ds.parameters["particles"]:
            self.pgrid_info = np.zeros((self.num_grids, 3), dtype='int64')
            return
        for fn in ['particle_position_%s' % ax for ax in 'xyz'] + \
                  ['particle_mass'] +  \
                  ['particle_velocity_%s' % ax for ax in 'xyz']:
            self.field_list.append(("io", fn))
        header = open(os.path.join(self.ds.output_dir, "DM", "Header"))
        version = header.readline()  # NOQA
        ndim = header.readline()  # NOQA
        nfields = header.readline()  # NOQA
        ntotalpart = int(header.readline())  # NOQA
        nextid = header.readline()  # NOQA
        maxlevel = int(header.readline())  # NOQA

        # Skip over how many grids on each level; this is degenerate
        for i in range(maxlevel + 1):
            header.readline()

        grid_info = np.fromiter((int(i) for line in header.readlines()
                                 for i in line.split()),
                                dtype='int64',
                                count=3*self.num_grids).reshape((self.num_grids, 3))
        # we need grid_info in `populate_grid_objects`, so save it to self

        for g, pg in izip(self.grids, grid_info):
            g.particle_filename = os.path.join(self.ds.output_dir, "DM",
                                               "Level_%s" % (g.Level),
                                               "DATA_%04i" % pg[0])
            g.NumberOfParticles = pg[1]
            g._particle_offset = pg[2]

        self.grid_particle_count[:, 0] = grid_info[:, 1]


class NyxDataset(BoxlibDataset):

    _index_class = NyxHierarchy

    @classmethod
    def _is_valid(cls, *args, **kwargs):
        # fill our args
        output_dir = args[0]
        # boxlib datasets are always directories
        if not os.path.isdir(output_dir): return False
        header_filename = os.path.join(output_dir, "Header")
        jobinfo_filename = os.path.join(output_dir, "job_info")
        if not os.path.exists(header_filename):
            # We *know* it's not boxlib if Header doesn't exist.
            return False
        if not os.path.exists(jobinfo_filename):
            return False
        # Now we check the job_info for the mention of maestro
        lines = open(jobinfo_filename).readlines()
        if any(line.startswith("Nyx  ") for line in lines): return True
        if any(line.startswith("nyx.") for line in lines): return True
        return False

    def _parse_parameter_file(self):
        super(NyxDataset, self)._parse_parameter_file()

        # Nyx is always cosmological.
        self.cosmological_simulation = 1

        jobinfo_filename = os.path.join(self.output_dir, "job_info")
        line = ""
        with open(jobinfo_filename, "r") as f:
            while not line.startswith(" Cosmology Information"):
                # get the code git hashes
                if "git hash" in line:
                    # line format: codename git hash:  the-hash
                    fields = line.split(":")
                    self.parameters[fields[0]] = fields[1].strip()
                line = next(f)

            # get the cosmology
            for line in f:
                if "Omega_m (comoving)" in line:
                    self.omega_matter = float(line.split(":")[1])
                elif "Omega_lambda (comoving)" in line:
                    self.omega_lambda = float(line.split(":")[1])
                elif "h (comoving)" in line:
                    self.hubble_constant = float(line.split(":")[1])


        # Read in the `comoving_a` file and parse the value. We should fix this
        # in the new Nyx output format...
        a_file = open(os.path.join(self.output_dir, "comoving_a"))
        a_string = a_file.readline().strip()
        a_file.close()

        # Set the scale factor and redshift
        self.cosmological_scale_factor = float(a_string)
        self.parameters["CosmologyCurrentRedshift"] = 1 / float(a_string) - 1

        # alias
        self.current_redshift = self.parameters["CosmologyCurrentRedshift"]
        if os.path.isdir(os.path.join(self.output_dir, "DM")):
            # we have particles
            self.parameters["particles"] = 1 
            self.particle_types = ("io",)
            self.particle_types_raw = self.particle_types

    def _set_code_unit_attributes(self):
        setdefaultattr(self, 'mass_unit', self.quan(1.0, "Msun"))
        setdefaultattr(self, 'time_unit', self.quan(1.0 / 3.08568025e19, "s"))
        setdefaultattr(self, 'length_unit',
                       self.quan(1.0 / (1 + self.current_redshift), "Mpc"))
        setdefaultattr(self, 'velocity_unit', self.length_unit / self.time_unit)

def _guess_pcast(vals):
    # Now we guess some things about the parameter and its type
    # Just in case there are multiple; we'll go
    # back afterward to using vals.
    v = vals.split()[0]
    try:
        float(v.upper().replace("D", "E"))
    except:
        pcast = str
        if v in ("F", "T"):
            pcast = bool
    else:
        syms = (".", "D+", "D-", "E+", "E-", "E", "D")
        if any(sym in v.upper() for sym in syms for v in vals.split()):
            pcast = float
        else:
            pcast = int
    vals = [pcast(value) for value in vals.split()]
    if len(vals) == 1:
        vals = vals[0]
    return vals
