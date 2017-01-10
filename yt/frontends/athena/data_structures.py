"""
Data structures for Athena.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
import os
import weakref
import glob

from yt.funcs import \
    mylog, \
    ensure_tuple
from yt.data_objects.grid_patch import \
    AMRGridPatch
from yt.geometry.grid_geometry_handler import \
    GridIndex
from yt.data_objects.static_output import \
    Dataset
from yt.utilities.lib.misc_utilities import \
    get_box_grids_level
from yt.geometry.geometry_handler import \
    YTDataChunk
from yt.extern.six import PY2

from .fields import AthenaFieldInfo
from yt.utilities.decompose import \
    decompose_array, get_psize

def chk23(strin):
    if PY2:
        return strin
    else:
        return strin.encode('utf-8')

def str23(strin):
    if PY2:
        return strin
    else:
        if isinstance(strin, list):
            return [s.decode('utf-8') for s in strin]
        else:
            return strin.decode('utf-8')

def check_readline(fl):
    line = fl.readline()
    chk = chk23("SCALARS")
    if chk in line and not line.startswith(chk):
        line = line[line.find(chk):]
    chk = chk23("VECTORS")
    if chk in line and not line.startswith(chk):
        line = line[line.find(chk):]
    return line

def check_break(line):
    splitup = line.strip().split()
    do_break = chk23('SCALAR') in splitup
    do_break = (chk23('VECTOR') in splitup) & do_break
    do_break = (chk23('TABLE') in splitup) & do_break
    do_break = (len(line) == 0) & do_break
    return do_break

def _get_convert(fname):
    def _conv(data):
        return data.convert(fname)
    return _conv

class AthenaGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, id, index, level, start, dimensions,
                 file_offset, read_dims):
        gname = index.grid_filenames[id]
        AMRGridPatch.__init__(self, id, filename = gname,
                              index = index)
        self.filename = gname
        self.Parent = []
        self.Children = []
        self.Level = level
        self.start_index = start.copy()
        self.stop_index = self.start_index + dimensions
        self.ActiveDimensions = dimensions.copy()
        self.file_offset = file_offset
        self.read_dims = read_dims

    def _setup_dx(self):
        # So first we figure out what the index is.  We don't assume
        # that dx=dy=dz , at least here.  We probably do elsewhere.
        id = self.id - self._id_offset
        if len(self.Parent) > 0:
            self.dds = self.Parent[0].dds / self.ds.refine_by
        else:
            LE, RE = self.index.grid_left_edge[id,:], \
                     self.index.grid_right_edge[id,:]
            self.dds = self.ds.arr((RE-LE)/self.ActiveDimensions, "code_length")
        if self.ds.dimensionality < 2: self.dds[1] = 1.0
        if self.ds.dimensionality < 3: self.dds[2] = 1.0
        self.field_data['dx'], self.field_data['dy'], self.field_data['dz'] = self.dds

    def __repr__(self):
        return "AthenaGrid_%04i (%s)" % (self.id, self.ActiveDimensions)

def parse_line(line, grid):
    # grid is a dictionary
    splitup = line.strip().split()
    if chk23("vtk") in splitup:
        grid['vtk_version'] = str23(splitup[-1])
    elif chk23("time=") in splitup:
        time_index = splitup.index(chk23("time="))
        grid['time'] = float(str23(splitup[time_index+1]).rstrip(','))
        grid['level'] = int(str23(splitup[time_index+3]).rstrip(','))
        grid['domain'] = int(str23(splitup[time_index+5]).rstrip(','))
    elif chk23("DIMENSIONS") in splitup:
        grid['dimensions'] = np.array(str23(splitup[-3:])).astype('int')
    elif chk23("ORIGIN") in splitup:
        grid['left_edge'] = np.array(str23(splitup[-3:])).astype('float64')
    elif chk23("SPACING") in splitup:
        grid['dds'] = np.array(str23(splitup[-3:])).astype('float64')
    elif chk23("CELL_DATA") in splitup or chk23("POINT_DATA") in splitup:
        grid["ncells"] = int(str23(splitup[-1]))
    elif chk23("SCALARS") in splitup:
        field = str23(splitup[1])
        grid['read_field'] = field
        grid['read_type'] = 'scalar'
    elif chk23("VECTORS") in splitup:
        field = str23(splitup[1])
        grid['read_field'] = field
        grid['read_type'] = 'vector'
    elif chk23("time") in splitup:
        time_index = splitup.index(chk23("time"))
        grid['time'] = float(str23(splitup[time_index+1]))

class AthenaHierarchy(GridIndex):

    grid = AthenaGrid
    _dataset_type='athena'
    _data_file = None

    def __init__(self, ds, dataset_type='athena'):
        self.dataset = weakref.proxy(ds)
        self.directory = os.path.dirname(self.dataset.filename)
        self.dataset_type = dataset_type
        # for now, the index file is the dataset!
        self.index_filename = os.path.join(os.getcwd(), self.dataset.filename)
        self._fhandle = open(self.index_filename,'rb')
        GridIndex.__init__(self, ds, dataset_type)

        self._fhandle.close()

    def _detect_output_fields(self):
        field_map = {}
        f = open(self.index_filename,'rb')
        line = check_readline(f)
        chkwhile = chk23('')
        while line != chkwhile:
            splitup = line.strip().split()
            chkd = chk23("DIMENSIONS")
            chkc = chk23("CELL_DATA")
            chkp = chk23("POINT_DATA")
            if chkd in splitup:
                field = str23(splitup[-3:])
                grid_dims = np.array(field).astype('int')
                line = check_readline(f)
            elif chkc in splitup or chkp in splitup:
                grid_ncells = int(str23(splitup[-1]))
                line = check_readline(f)
                if np.prod(grid_dims) != grid_ncells:
                    grid_dims -= 1
                    grid_dims[grid_dims==0]=1
                if np.prod(grid_dims) != grid_ncells:
                    mylog.error('product of dimensions %i not equal to number of cells %i' %
                          (np.prod(grid_dims), grid_ncells))
                    raise TypeError
                break
            else:
                line = check_readline(f)
        read_table = False
        read_table_offset = f.tell()
        while line != chkwhile:
            splitup = line.strip().split()
            chks = chk23('SCALARS')
            chkv = chk23('VECTORS')
            if chks in line and chks not in splitup:
                splitup = str23(line[line.find(chks):].strip().split())
            if chkv in line and chkv not in splitup:
                splitup = str23(line[line.find(chkv):].strip().split())
            if chks in splitup:
                field = ("athena", str23(splitup[1]))
                dtype = str23(splitup[-1]).lower()
                if not read_table:
                    line = check_readline(f) # Read the lookup table line
                    read_table = True
                field_map[field] = ('scalar', f.tell() - read_table_offset, dtype)
                read_table=False
            elif chkv in splitup:
                field = str23(splitup[1])
                dtype = str23(splitup[-1]).lower()
                for ax in 'xyz':
                    field_map[("athena","%s_%s" % (field, ax))] =\
                            ('vector', f.tell() - read_table_offset, dtype)
            line = check_readline(f)

        f.close()

        self.field_list = list(field_map.keys())
        self._field_map = field_map

    def _count_grids(self):
        self.num_grids = self.dataset.nvtk*self.dataset.nprocs

    def _parse_index(self):
        f = open(self.index_filename,'rb')
        grid = {}
        grid['read_field'] = None
        grid['read_type'] = None
        line = f.readline()
        while grid['read_field'] is None:
            parse_line(line, grid)
            if check_break(line): break
            line = f.readline()
        f.close()

        # It seems some datasets have a mismatch between ncells and
        # the actual grid dimensions.
        if np.prod(grid['dimensions']) != grid['ncells']:
            grid['dimensions'] -= 1
            grid['dimensions'][grid['dimensions']==0]=1
        if np.prod(grid['dimensions']) != grid['ncells']:
            mylog.error('product of dimensions %i not equal to number of cells %i' %
                  (np.prod(grid['dimensions']), grid['ncells']))
            raise TypeError

        # Need to determine how many grids: self.num_grids
        dataset_dir = os.path.dirname(self.index_filename)
        dname = os.path.split(self.index_filename)[-1]
        if dataset_dir.endswith("id0"):
            dname = "id0/"+dname
            dataset_dir = dataset_dir[:-3]

        gridlistread = glob.glob(os.path.join(dataset_dir, 'id*/%s-id*%s' % (dname[4:-9],dname[-9:])))
        gridlistread.insert(0,self.index_filename)
        if 'id0' in dname:
            gridlistread += glob.glob(os.path.join(dataset_dir, 'id*/lev*/%s*-lev*%s' % (dname[4:-9],dname[-9:])))
        else :
            gridlistread += glob.glob(os.path.join(dataset_dir, 'lev*/%s*-lev*%s' % (dname[:-9],dname[-9:])))
        ndots = dname.count(".")
        gridlistread = [fn for fn in gridlistread if os.path.basename(fn).count(".") == ndots]
        self.num_grids = len(gridlistread)
        dxs=[]
        levels = np.zeros(self.num_grids, dtype='int32')
        glis = np.empty((self.num_grids,3), dtype='float64')
        gdds = np.empty((self.num_grids,3), dtype='float64')
        gdims = np.ones_like(glis)
        j = 0
        self.grid_filenames = gridlistread
        while j < (self.num_grids):
            f = open(gridlistread[j],'rb')
            gridread = {}
            gridread['read_field'] = None
            gridread['read_type'] = None
            line = f.readline()
            while gridread['read_field'] is None:
                parse_line(line, gridread)
                splitup = line.strip().split()
                if chk23('X_COORDINATES') in splitup:
                    gridread['left_edge'] = np.zeros(3)
                    gridread['dds'] = np.zeros(3)
                    v = np.fromfile(f, dtype='>f8', count=2)
                    gridread['left_edge'][0] = v[0]-0.5*(v[1]-v[0])
                    gridread['dds'][0] = v[1]-v[0]
                if chk23('Y_COORDINATES') in splitup:
                    v = np.fromfile(f, dtype='>f8', count=2)
                    gridread['left_edge'][1] = v[0]-0.5*(v[1]-v[0])
                    gridread['dds'][1] = v[1]-v[0]
                if chk23('Z_COORDINATES') in splitup:
                    v = np.fromfile(f, dtype='>f8', count=2)
                    gridread['left_edge'][2] = v[0]-0.5*(v[1]-v[0])
                    gridread['dds'][2] = v[1]-v[0]
                if check_break(line): break
                line = f.readline()
            f.close()
            levels[j] = gridread.get('level', 0)
            glis[j,0] = gridread['left_edge'][0]
            glis[j,1] = gridread['left_edge'][1]
            glis[j,2] = gridread['left_edge'][2]
            # It seems some datasets have a mismatch between ncells and
            # the actual grid dimensions.
            if np.prod(gridread['dimensions']) != gridread['ncells']:
                gridread['dimensions'] -= 1
                gridread['dimensions'][gridread['dimensions']==0]=1
            if np.prod(gridread['dimensions']) != gridread['ncells']:
                mylog.error('product of dimensions %i not equal to number of cells %i' %
                      (np.prod(gridread['dimensions']), gridread['ncells']))
                raise TypeError
            gdims[j,0] = gridread['dimensions'][0]
            gdims[j,1] = gridread['dimensions'][1]
            gdims[j,2] = gridread['dimensions'][2]
            # Setting dds=1 for non-active dimensions in 1D/2D datasets
            gridread['dds'][gridread['dimensions']==1] = 1.
            gdds[j,:] = gridread['dds']

            j=j+1

        gres = glis + gdims*gdds
        # Now we convert the glis, which were left edges (floats), to indices
        # from the domain left edge.  Then we do a bunch of fixing now that we
        # know the extent of all the grids.
        glis = np.round((glis - self.dataset.domain_left_edge.ndarray_view())/gdds).astype('int')
        new_dre = np.max(gres,axis=0)
        dre_units = self.dataset.domain_right_edge.uq
        self.dataset.domain_right_edge = np.round(new_dre, decimals=12)*dre_units
        self.dataset.domain_width = \
                (self.dataset.domain_right_edge -
                 self.dataset.domain_left_edge)
        self.dataset.domain_center = \
                0.5*(self.dataset.domain_left_edge +
                     self.dataset.domain_right_edge)
        self.dataset.domain_dimensions = \
                np.round(self.dataset.domain_width/gdds[0]).astype('int')

        if self.dataset.dimensionality <= 2 :
            self.dataset.domain_dimensions[2] = np.int(1)
        if self.dataset.dimensionality == 1 :
            self.dataset.domain_dimensions[1] = np.int(1)

        dle = self.dataset.domain_left_edge
        dre = self.dataset.domain_right_edge
        dx_root = (self.dataset.domain_right_edge-
                   self.dataset.domain_left_edge)/self.dataset.domain_dimensions

        if self.dataset.nprocs > 1:
            gle_all = []
            gre_all = []
            shapes_all = []
            levels_all = []
            new_gridfilenames = []
            file_offsets = []
            read_dims = []
            for i in range(levels.shape[0]):
                dx = dx_root/self.dataset.refine_by**(levels[i])
                gle_orig = self.ds.arr(np.round(dle + dx*glis[i], decimals=12),
                                       "code_length")
                gre_orig = self.ds.arr(np.round(gle_orig + dx*gdims[i], decimals=12),
                                       "code_length")
                bbox = np.array([[le,re] for le, re in zip(gle_orig, gre_orig)])
                psize = get_psize(self.ds.domain_dimensions, self.ds.nprocs)
                gle, gre, shapes, slices = decompose_array(gdims[i], psize, bbox)
                gle_all += gle
                gre_all += gre
                shapes_all += shapes
                levels_all += [levels[i]]*self.dataset.nprocs
                new_gridfilenames += [self.grid_filenames[i]]*self.dataset.nprocs
                file_offsets += [[slc[0].start, slc[1].start, slc[2].start] for slc in slices]
                read_dims += [np.array([gdims[i][0], gdims[i][1], shape[2]], dtype="int") for shape in shapes]
            self.num_grids *= self.dataset.nprocs
            self.grids = np.empty(self.num_grids, dtype='object')
            self.grid_filenames = new_gridfilenames
            self.grid_left_edge = self.ds.arr(gle_all, "code_length")
            self.grid_right_edge = self.ds.arr(gre_all, "code_length")
            self.grid_dimensions = np.array([shape for shape in shapes_all],
                                            dtype="int32")
            gdds = (self.grid_right_edge-self.grid_left_edge)/self.grid_dimensions
            glis = np.round((self.grid_left_edge - self.ds.domain_left_edge)/gdds).astype('int')
            for i in range(self.num_grids):
                self.grids[i] = self.grid(i,self,levels_all[i],
                                          glis[i], shapes_all[i],
                                          file_offsets[i], read_dims[i])
        else:
            self.grids = np.empty(self.num_grids, dtype='object')
            for i in range(levels.shape[0]):
                self.grids[i] = self.grid(i,self,levels[i],
                                          glis[i], gdims[i], [0]*3,
                                          gdims[i])
                dx = dx_root/self.dataset.refine_by**(levels[i])
                dxs.append(dx)

            dx = self.ds.arr(dxs, "code_length")
            self.grid_left_edge = self.ds.arr(np.round(dle + dx*glis, decimals=12),
                                              "code_length")
            self.grid_dimensions = gdims.astype("int32")
            self.grid_right_edge = self.ds.arr(np.round(self.grid_left_edge +
                                                        dx*self.grid_dimensions,
                                                        decimals=12),
                                               "code_length")
        if self.dataset.dimensionality <= 2:
            self.grid_right_edge[:,2] = dre[2]
        if self.dataset.dimensionality == 1:
            self.grid_right_edge[:,1:] = dre[1:]
        self.grid_particle_count = np.zeros([self.num_grids, 1], dtype='int64')

    def _populate_grid_objects(self):
        for g in self.grids:
            g._prepare_grid()
            g._setup_dx()
        self._reconstruct_parent_child()
        self.max_level = self.grid_levels.max()

    def _reconstruct_parent_child(self):
        mask = np.empty(len(self.grids), dtype='int32')
        mylog.debug("First pass; identifying child grids")
        for i, grid in enumerate(self.grids):
            get_box_grids_level(self.grid_left_edge[i,:],
                                self.grid_right_edge[i,:],
                                self.grid_levels[i] + 1,
                                self.grid_left_edge, self.grid_right_edge,
                                self.grid_levels, mask)
            grid.Children = [g for g in self.grids[mask.astype("bool")] if g.Level == grid.Level + 1]
        mylog.debug("Second pass; identifying parents")
        for i, grid in enumerate(self.grids): # Second pass
            for child in grid.Children:
                child.Parent.append(grid)

    def _get_grid_children(self, grid):
        mask = np.zeros(self.num_grids, dtype='bool')
        grids, grid_ind = self.get_box_grids(grid.LeftEdge, grid.RightEdge)
        mask[grid_ind] = True
        return [g for g in self.grids[mask] if g.Level == grid.Level + 1]

    def _chunk_io(self, dobj, cache = True, local_only = False):
        gobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for subset in gobjs:
            yield YTDataChunk(dobj, "io", [subset],
                              self._count_selection(dobj, [subset]),
                              cache = cache)

class AthenaDataset(Dataset):
    _index_class = AthenaHierarchy
    _field_info_class = AthenaFieldInfo
    _dataset_type = "athena"

    def __init__(self, filename, dataset_type='athena',
                 storage_filename=None, parameters=None,
                 units_override=None, nprocs=1, unit_system="cgs"):
        self.fluid_types += ("athena",)
        self.nprocs = nprocs
        if parameters is None:
            parameters = {}
        self.specified_parameters = parameters.copy()
        if units_override is None:
            units_override = {}
        # This is for backwards-compatibility
        already_warned = False
        for k, v in list(self.specified_parameters.items()):
            if k.endswith("_unit") and k not in units_override:
                if not already_warned:
                    mylog.warning("Supplying unit conversions from the parameters dict is deprecated, "+
                                  "and will be removed in a future release. Use units_override instead.")
                    already_warned = True
                units_override[k] = self.specified_parameters.pop(k)
        Dataset.__init__(self, filename, dataset_type, units_override=units_override,
                         unit_system=unit_system)
        self.filename = filename
        if storage_filename is None:
            storage_filename = '%s.yt' % filename.split('/')[-1]
        self.storage_filename = storage_filename
        self.backup_filename = self.filename[:-4] + "_backup.gdf"
        # Unfortunately we now have to mandate that the index gets
        # instantiated so that we can make sure we have the correct left
        # and right domain edges.
        self.index

    def _set_code_unit_attributes(self):
        """
        Generates the conversion to various physical _units based on the
        parameter file
        """
        if "length_unit" not in self.units_override:
            self.no_cgs_equiv_length = True
        for unit, cgs in [("length", "cm"), ("time", "s"), ("mass", "g")]:
            # We set these to cgs for now, but they may have been overriden
            if getattr(self, unit+'_unit', None) is not None:
                continue
            mylog.warning("Assuming 1.0 = 1.0 %s", cgs)
            setattr(self, "%s_unit" % unit, self.quan(1.0, cgs))
        self.magnetic_unit = np.sqrt(4*np.pi * self.mass_unit /
                                     (self.time_unit**2 * self.length_unit))
        self.magnetic_unit.convert_to_units("gauss")
        self.velocity_unit = self.length_unit / self.time_unit

    def _parse_parameter_file(self):
        self._handle = open(self.parameter_filename, "rb")
        # Read the start of a grid to get simulation parameters.
        grid = {}
        grid['read_field'] = None
        line = self._handle.readline()
        while grid['read_field'] is None:
            parse_line(line, grid)
            splitup = line.strip().split()
            if chk23('X_COORDINATES') in splitup:
                grid['left_edge'] = np.zeros(3)
                grid['dds'] = np.zeros(3)
                v = np.fromfile(self._handle, dtype='>f8', count=2)
                grid['left_edge'][0] = v[0]-0.5*(v[1]-v[0])
                grid['dds'][0] = v[1]-v[0]
            if chk23('Y_COORDINATES') in splitup:
                v = np.fromfile(self._handle, dtype='>f8', count=2)
                grid['left_edge'][1] = v[0]-0.5*(v[1]-v[0])
                grid['dds'][1] = v[1]-v[0]
            if chk23('Z_COORDINATES') in splitup:
                v = np.fromfile(self._handle, dtype='>f8', count=2)
                grid['left_edge'][2] = v[0]-0.5*(v[1]-v[0])
                grid['dds'][2] = v[1]-v[0]
            if check_break(line): break
            line = self._handle.readline()

        self.domain_left_edge = grid['left_edge']
        mylog.info("Temporarily setting domain_right_edge = -domain_left_edge."+
                  " This will be corrected automatically if it is not the case.")
        self.domain_right_edge = -self.domain_left_edge
        self.domain_width = self.domain_right_edge-self.domain_left_edge
        self.domain_dimensions = np.round(self.domain_width/grid['dds']).astype('int32')
        refine_by = None
        if refine_by is None: refine_by = 2
        self.refine_by = refine_by
        dimensionality = 3
        if grid['dimensions'][2] == 1 :
            dimensionality = 2
        if grid['dimensions'][1] == 1 :
            dimensionality = 1
        if dimensionality <= 2 : self.domain_dimensions[2] = np.int32(1)
        if dimensionality == 1 : self.domain_dimensions[1] = np.int32(1)
        if dimensionality != 3 and self.nprocs > 1:
            raise RuntimeError("Virtual grids are only supported for 3D outputs!")
        self.dimensionality = dimensionality
        self.current_time = grid["time"]
        self.unique_identifier = self.parameter_filename.__hash__()
        self.cosmological_simulation = False
        self.num_ghost_zones = 0
        self.field_ordering = 'fortran'
        self.boundary_conditions = [1]*6
        if 'periodicity' in self.specified_parameters:
            self.periodicity = ensure_tuple(self.specified_parameters['periodicity'])
        else:
            self.periodicity = (True,True,True,)
        if 'gamma' in self.specified_parameters:
            self.gamma = float(self.specified_parameters['gamma'])
        else:
            self.gamma = 5./3.
        dataset_dir = os.path.dirname(self.parameter_filename)
        dname = os.path.split(self.parameter_filename)[-1]
        if dataset_dir.endswith("id0"):
            dname = "id0/"+dname
            dataset_dir = dataset_dir[:-3]

        gridlistread = glob.glob(os.path.join(dataset_dir, 'id*/%s-id*%s' % (dname[4:-9],dname[-9:])))
        if 'id0' in dname :
            gridlistread += glob.glob(os.path.join(dataset_dir, 'id*/lev*/%s*-lev*%s' % (dname[4:-9],dname[-9:])))
        else :
            gridlistread += glob.glob(os.path.join(dataset_dir, 'lev*/%s*-lev*%s' % (dname[:-9],dname[-9:])))
        ndots = dname.count(".")
        gridlistread = [fn for fn in gridlistread if os.path.basename(fn).count(".") == ndots]
        self.nvtk = len(gridlistread)+1

        self.current_redshift = self.omega_lambda = self.omega_matter = \
            self.hubble_constant = self.cosmological_simulation = 0.0
        self.parameters['Time'] = self.current_time # Hardcode time conversion for now.
        self.parameters["HydroMethod"] = 0 # Hardcode for now until field staggering is supported.
        if "gamma" in self.specified_parameters:
            self.parameters["Gamma"] = self.specified_parameters["gamma"]
        else:
            self.parameters["Gamma"] = 5./3.
        self.geometry = self.specified_parameters.get("geometry", "cartesian")
        self._handle.close()

    @classmethod
    def _is_valid(self, *args, **kwargs):
        try:
            if 'vtk' in args[0]:
                return True
        except:
            pass
        return False

    @property
    def _skip_cache(self):
        return True

    def __repr__(self):
        return self.basename.rsplit(".", 1)[0]
