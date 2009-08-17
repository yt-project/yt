"""
AMR hierarchy container class

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2009 Matthew Turk.  All Rights Reserved.

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
import string, re, gc, time
import cPickle
from itertools import chain, izip
#import yt.enki

_data_style_funcs = \
   { 'enzo_hdf4': (readDataHDF4,readAllDataHDF4, getFieldsHDF4, readDataSliceHDF4,
         getExceptionHDF4, DataQueueHDF4),
     'enzo_hdf4_2d': (readDataHDF4, readAllDataHDF4, getFieldsHDF4, readDataSliceHDF4_2D,
         getExceptionHDF4, DataQueueHDF4_2D),
     'enzo_hdf5': (readDataHDF5, readAllDataHDF5, getFieldsHDF5, readDataSliceHDF5,
         getExceptionHDF5, DataQueueHDF5),
     'enzo_packed_3d': (readDataPacked, readAllDataPacked, getFieldsPacked, readDataSlicePacked,
         getExceptionHDF5, DataQueuePackedHDF5),
     'orion_native': (readDataNative, readAllDataNative, None, readDataSliceNative,
         getExceptionHDF5, DataQueueNative), \
     'enzo_inline': (readDataInMemory, readAllDataInMemory, getFieldsInMemory, readDataSliceInMemory,
         getExceptionInMemory, DataQueueInMemory),
     'enzo_packed_2d': (readDataPacked, readAllDataPacked, getFieldsPacked, readDataSlicePacked2D,
         getExceptionHDF5, DataQueuePacked2D),
     'enzo_packed_1d': (readDataPacked, readAllDataPacked, getFieldsPacked, readDataSlicePacked1D,
         getExceptionHDF5, DataQueuePacked1D),
   }

class OldAMRHierarchy:
    _data_mode = None # Default
    def __init__(self, pf):
        self.parameter_file = weakref.proxy(pf)
        self._max_locations = {}
        self._data_file = None
        self._setup_classes()
        self._initialize_grids()

        # For use with derived quantities depending on centers
        # Although really, I think perhaps we should take a closer look
        # at how "center" is used.
        self.center = None

        self._initialize_level_stats()

        mylog.debug("Initializing data file")
        self._initialize_data_file()
        mylog.debug("Populating hierarchy")
        self._populate_hierarchy()
        mylog.debug("Done populating hierarchy")
        self._add_detected_fields()
        mylog.debug("Done adding auto-detected fields")

    def __getitem__(self, item):
        return self.parameter_file[item]

    @time_execution
    def export_particles_pb(self, filename, filter = 1, indexboundary = 0, fields = None, scale=1.0):
        """
        Exports all the star particles, or a subset, to pb-format *filename*
        for viewing in partiview.  Filters based on particle_type=*filter*,
        particle_index>=*indexboundary*, and exports *fields*, if supplied.
        Otherwise, index, position(x,y,z).  Optionally *scale* by a given
        factor before outputting.
        """
        import struct
        pbf_magic = 0xffffff98
        header_fmt = 'Iii'
        fmt = 'ifff'
        f = open(filename,"w")
        if fields:
            fmt += len(fields)*'f'
            padded_fields = string.join(fields,"\0") + "\0"
            header_fmt += "%ss" % len(padded_fields)
            args = [pbf_magic, struct.calcsize(header_fmt), len(fields), padded_fields]
            fields = ["particle_index","particle_position_x","particle_position_y","particle_position_z"] \
                   + fields
            format = 'Int32,Float32,Float32,Float32' + ',Float32'*(len(fields)-4)
        else:
            args = [pbf_magic, struct.calcsize(header_fmt), 0]
            fields = ["particle_index","particle_position_x","particle_position_y","particle_position_z"]
            format = 'Int32,Float32,Float32,Float32'
        f.write(struct.pack(header_fmt, *args))
        tot = 0
        sc = na.array([1.0] + [scale] * 3 + [1.0]*(len(fields)-4))
        gI = na.where(self.gridNumberOfParticles.ravel() > 0)
        for g in self.grids[gI]:
            pI = na.where(na.logical_and((g["particle_type"] == filter),(g["particle_index"] >= indexboundary)) == 1)
            tot += pI[0].shape[0]
            toRec = []
            for field, scale in zip(fields, sc):
                toRec.append(scale*g[field][pI])
            particle_info = rec.array(toRec,formats=format)
            particle_info.tofile(f)
        f.close()
        mylog.info("Wrote %s particles to %s", tot, filename)

    @time_execution
    def export_boxes_pv(self, filename):
        """
        Exports the grid structure in partiview text format.
        """
        f=open(filename,"w")
        for l in xrange(self.maxLevel):
            f.write("add object g%s = l%s\n" % (l,l))
            ind = self._select_level(l)
            for i in ind:
                f.write("add box -n %s -l %s %s,%s %s,%s %s,%s\n" % \
                    (i+1, self.gridLevels.ravel()[i],
                     self.gridLeftEdge[i,0], self.gridRightEdge[i,0],
                     self.gridLeftEdge[i,1], self.gridRightEdge[i,1],
                     self.gridLeftEdge[i,2], self.gridRightEdge[i,2]))

    def __initialize_octree_list(self, fields):
        import DepthFirstOctree as dfo
        o_length = r_length = 0
        grids = []
        levels_finest, levels_all = defaultdict(lambda: 0), defaultdict(lambda: 0)
        for g in self.grids:
            ff = na.array([g[f] for f in fields])
            grids.append(dfo.OctreeGrid(
                            g.child_index_mask.astype('int32'),
                            ff.astype("float64"),
                            g.LeftEdge.astype('float64'),
                            g.ActiveDimensions.astype('int32'),
                            na.ones(1,dtype='float64') * g.dds[0], g.Level))
            levels_all[g.Level] += g.ActiveDimensions.prod()
            levels_finest[g.Level] += g.child_mask.ravel().sum()
            g.clear_data()
        ogl = dfo.OctreeGridList(grids)
        return ogl, levels_finest, levels_all

    def _generate_flat_octree(self, fields):
        """
        Generates two arrays, one of the actual values in a depth-first flat
        octree array, and the other of the values describing the refinement.
        This allows for export to a code that understands this.  *field* is the
        field used in the data array.
        """
        import DepthFirstOctree as dfo
        fields = ensure_list(fields)
        ogl, levels_finest, levels_all = self.__initialize_octree_list(fields)
        o_length = na.sum(levels_finest.values())
        r_length = na.sum(levels_all.values())
        output = na.zeros((o_length,len(fields)), dtype='float64')
        refined = na.zeros(r_length, dtype='int32')
        position = dfo.position()
        dfo.RecurseOctreeDepthFirst(0, 0, 0,
                ogl[0].dimensions[0],
                ogl[0].dimensions[1],
                ogl[0].dimensions[2],
                position, 1,
                output, refined, ogl)
        dd = {}
        for i,field in enumerate(fields): dd[field] = output[:,i]
        return dd, refined

    def _generate_levels_octree(self, fields):
        import DepthFirstOctree as dfo
        fields = ensure_list(fields) + ["Ones", "Ones"]
        ogl, levels_finest, levels_all = self.__initialize_octree_list(fields)
        o_length = na.sum(levels_finest.values())
        r_length = na.sum(levels_all.values())
        output = na.zeros((r_length,len(fields)), dtype='float64')
        genealogy = na.zeros((r_length, 3), dtype='int32') - 1 # init to -1
        corners = na.zeros((r_length, 3), dtype='float64')
        position = na.add.accumulate(
                    na.array([0] + [levels_all[v] for v in
                        sorted(levels_all)[:-1]], dtype='int32'))
        pp = position.copy()
        dfo.RecurseOctreeByLevels(0, 0, 0,
                ogl[0].dimensions[0],
                ogl[0].dimensions[1],
                ogl[0].dimensions[2],
                position.astype('int32'), 1,
                output, genealogy, corners, ogl)
        return output, genealogy, levels_all, levels_finest, pp, corners

class AMRHierarchy(ObjectFindingMixin):
    def __init__(self, pf, data_style):
        self.parameter_file = weakref.proxy(pf)
        self.pf = self.parameter_file

        self._initialize_state_variables()

        mylog.debug("Initializing data storage.")
        self._initialize_data_storage()

        mylog.debug("Counting grids.")
        self._count_grids()

        # Must be defined in subclass
        mylog.debug("Setting up classes.")
        self._setup_classes()

        mylog.debug("Counting grids.")
        self._initialize_grid_arrays()

        mylog.debug("Parsing hierarchy.")
        self._parse_hierarchy()

        mylog.debug("Constructing grid objects.")
        self._populate_grid_objects()

        mylog.debug("Detecting fields.")
        self.field_list = []
        self._detect_fields()

        mylog.debug("Adding unknown detected fields")
        self._setup_unknown_fields()

        mylog.debug("Setting up derived fields")
        self._setup_derived_fields()

        mylog.debug("Initializing data grid data IO")
        self._setup_data_queue()

        mylog.debug("Re-examining hierarchy")
        self._initialize_level_stats()

    def _get_parameters(self):
        return self.parameter_file.parameters
    parameters=property(_get_parameters)

    def select_grids(self, level):
        """
        Returns an array of grids at *level*.
        """
        return self.grids[self.grid_levels.flat == level]

    def _initialize_state_variables(self):
        self._parallel_locking = False
        self._data_file = None
        self._data_mode = None
        self._max_locations = {}

    def _setup_classes(self, dd):
        # Called by subclass
        self.object_types = []
        self.objects = []
        for name, cls in sorted(data_object_registry.items()):
            cname = cls.__name__
            if cname.endswith("Base"): cname = cname[:-4]
            self._add_object_class(name, cname, cls, dd)
        self.object_types.sort()

    # Now all the object related stuff

    def all_data(self, find_max=False):
        pf = self.parameter_file
        if find_max: c = self.find_max("Density")[1]
        else: c = (pf["DomainRightEdge"] + pf["DomainLeftEdge"])/2.0
        return self.region(c, 
            pf["DomainLeftEdge"], pf["DomainRightEdge"])

    def _join_field_lists(self, field_list):
        return field_list

    def _get_data_reader_dict(self):
        dd = { 'readDataFast' : _data_style_funcs[self.data_style][0],
               'readAllData' : _data_style_funcs[self.data_style][1],
               'getFields' : _data_style_funcs[self.data_style][2],
               'readDataSlice' : _data_style_funcs[self.data_style][3],
               '_read_data' : _data_style_funcs[self.data_style][0],
               '_read_all_data' : _data_style_funcs[self.data_style][1],
               '_read_field_names' : _data_style_funcs[self.data_style][2],
               '_read_data_slice' : _data_style_funcs[self.data_style][3],
               '_read_exception' : _data_style_funcs[self.data_style][4](),
               'pf' : self.parameter_file, # Already weak
               'hierarchy': weakref.proxy(self) }
        return dd

    def _initialize_data_storage(self):
        if not ytcfg.getboolean('lagos','serialize'): return
        if os.path.isfile(os.path.join(self.directory,
                            "%s.yt" % self["CurrentTimeIdentifier"])):
            fn = os.path.join(self.directory,"%s.yt" % self["CurrentTimeIdentifier"])
        else:
            fn = os.path.join(self.directory,
                    "%s.yt" % self.parameter_file.basename)
        if os.path.isfile(fn):
            writable = os.access(fn, os.W_OK)
        else:
            writable = os.access(self.directory, os.W_OK)
        if ytcfg.getboolean('lagos','onlydeserialize') or not writable:
            self._data_mode = mode = 'r'
        else:
            self._data_mode = mode = 'a'
        self.__create_data_file(fn)
        self.__data_filename = fn
        self._data_file = h5py.File(fn, self._data_mode)

    @parallel_root_only
    def __create_data_file(self, fn):
        f = h5py.File(fn, 'a')
        f.close()

    def _setup_data_queue(self):
        self.queue = _data_style_funcs[self.data_style][5]()

    def _save_data(self, array, node, name, set_attr=None, force=False, passthrough = False):
        """
        Arbitrary numpy data will be saved to the region in the datafile
        described by *node* and *name*.  If data file does not exist, it throws
        no error and simply does not save.
        """

        if self._data_mode != 'a': return
        if "ArgsError" in dir(h5py.h5):
            exception = h5py.h5.ArgsError
        else:
            exception = h5py.h5.H5Error
        try:
            node_loc = self._data_file[node]
            if name in node_loc.listnames() and force:
                mylog.info("Overwriting node %s/%s", node, name)
                del self._data_file[node][name]
            elif name in node_loc.listnames() and passthrough:
                return
        except exception:
            pass
        myGroup = self._data_file['/']
        for q in node.split('/'):
            if q: myGroup = myGroup.require_group(q)
        arr = myGroup.create_dataset(name,data=array)
        if set_attr is not None:
            for i, j in set_attr.items(): arr.attrs[i] = j
        self._data_file.flush()

    def _reload_data_file(self, *args, **kwargs):
        if self._data_file is None: return
        self._data_file.close()
        del self._data_file
        self._data_file = h5py.File(self.__data_filename, self._data_mode)

    def _reset_save_data(self,round_robin=False):
        if round_robin:
            self.save_data = self._save_data
        else:
            self.save_data = parallel_splitter(self._save_data, self._reload_data_file)
    
    save_data = parallel_splitter(_save_data, _reload_data_file)

    def save_object(self, obj, name):
        s = cPickle.dumps(obj, protocol=-1)
        self.save_data(s, "/Objects", name, force = True)

    def load_object(self, name):
        obj = self.get_data("/Objects", name)
        if obj is None:
            return
        obj = cPickle.loads(obj.value)
        if iterable(obj) and len(obj) == 2:
            obj = obj[1] # Just the object, not the pf
        if hasattr(obj, '_fix_pickle'): obj._fix_pickle()
        return obj

    def get_data(self, node, name):
        """
        Return the dataset with a given *name* located at *node* in the
        datafile.
        """
        if self._data_file == None:
            return None
        if node[0] != "/": node = "/%s" % node

        myGroup = self._data_file['/']
        for group in node.split('/'):
            if group:
                if group not in myGroup.listnames():
                    return None
                myGroup = myGroup[group]
        if name not in myGroup.listnames():
            return None

        full_name = "%s/%s" % (node, name)
        return self._data_file[full_name]

    def _close_data_file(self):
        if self._data_file:
            self._data_file.close()
            del self._data_file
            self._data_file = None

    def get_smallest_dx(self):
        """
        Returns (in code units) the smallest cell size in the simulation.
        """
        return self.select_grids(self.grid_levels.max())[0].dds[0]

    def _add_object_class(self, name, class_name, base, dd):
        self.object_types.append(name)
        obj = classobj(class_name, (base,), dd)
        setattr(self, name, obj)

    def _initialize_level_stats(self):
        # Now some statistics:
        #   0 = number of grids
        #   1 = number of cells
        #   2 = blank
        desc = {'names': ['numgrids','numcells','level'],
                'formats':['Int32']*3}
        self.level_stats = blankRecordArray(desc, MAXLEVEL)
        self.level_stats['level'] = [i for i in range(MAXLEVEL)]
        self.level_stats['numgrids'] = [0 for i in range(MAXLEVEL)]
        self.level_stats['numcells'] = [0 for i in range(MAXLEVEL)]
        for level in xrange(self.max_level+1):
            self.level_stats[level]['numgrids'] = na.sum(self.grid_levels == level)
            li = (self.grid_levels[:,0] == level)
            self.level_stats[level]['numcells'] = self.grid_dimensions[li,:].prod(axis=1).sum()

    @property
    def grid_corners(self):
        return na.array([
          [self.grid_left_edge[:,0], self.grid_left_edge[:,1], self.grid_left_edge[:,2]],
          [self.grid_right_edge[:,0], self.grid_left_edge[:,1], self.grid_left_edge[:,2]],
          [self.grid_right_edge[:,0], self.grid_right_edge[:,1], self.grid_left_edge[:,2]],
          [self.grid_right_edge[:,0], self.grid_right_edge[:,1], self.grid_right_edge[:,2]],
          [self.grid_left_edge[:,0], self.grid_right_edge[:,1], self.grid_right_edge[:,2]],
          [self.grid_left_edge[:,0], self.grid_left_edge[:,1], self.grid_right_edge[:,2]],
          [self.grid_right_edge[:,0], self.grid_left_edge[:,1], self.grid_right_edge[:,2]],
          [self.grid_left_edge[:,0], self.grid_right_edge[:,1], self.grid_left_edge[:,2]],
        ], dtype='float64')

    def print_stats(self):
        """
        Prints out (stdout) relevant information about the simulation
        """
        for level in xrange(MAXLEVEL):
            if (self.level_stats['numgrids'][level]) == 0:
                break
            print "% 3i\t% 6i\t% 11i" % \
                  (level, self.level_stats['numgrids'][level],
                   self.level_stats['numcells'][level])
            dx = self.select_grids(level)[0].dds[0]
        print "-" * 28
        print "   \t% 6i\t% 11i" % (self.level_stats['numgrids'].sum(), self.level_stats['numcells'].sum())
        print "\n"
        try:
            print "z = %0.8f" % (self["CosmologyCurrentRedshift"])
        except:
            pass
        t_s = self.pf["InitialTime"] * self.pf["Time"]
        print "t = %0.8e = %0.8e s = %0.8e years" % \
            (self.pf["InitialTime"], \
             t_s, t_s / (365*24*3600.0) )
        print "\nSmallest Cell:"
        u=[]
        for item in self.parameter_file.units.items():
            u.append((item[1],item[0]))
        u.sort()
        for unit in u:
            print "\tWidth: %0.3e %s" % (dx*unit[0], unit[1])


class EnzoHierarchy(AMRHierarchy):

    _strip_path = False
    def __init__(self, pf, data_style):
        
        self.data_style = data_style
        self.hierarchy_filename = os.path.abspath(
            "%s.hierarchy" % (pf.parameter_filename))
        if os.path.getsize(self.hierarchy_filename) == 0:
            raise IOError(-1,"File empty", self.hierarchy_filename)
        self.directory = os.path.dirname(self.hierarchy_filename)

        # For some reason, r8 seems to want Float64
        if pf.has_key("CompilerPrecision") \
            and pf["CompilerPrecision"] == "r4":
            self.float_type = 'float32'
        else:
            self.float_type = 'float64'

        AMRHierarchy.__init__(self, pf, data_style)
        # sync it back
        self.parameter_file.data_style = self.data_style

    def _initialize_data_storage(self):
        pass

    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        AMRHierarchy._setup_classes(self, dd)
        self._add_object_class('grid', "EnzoGrid", EnzoGridBase, dd)
        self.object_types.sort()

    def _count_grids(self):
        test_grid = test_grid_id = None
        for line in rlines(open(self.hierarchy_filename, "rb")):
            if line.startswith("BaryonFileName") or \
               line.startswith("FileName "):
                test_grid = line.split("=")[-1].strip().rstrip()
            if line.startswith("Grid "):
                self.num_grids = test_grid_id = int(line.split("=")[-1])
                break
        self._guess_data_style(self.pf["TopGridRank"], test_grid, test_grid_id)

    def _initialize_grid_arrays(self):
        mylog.debug("Allocating arrays for %s grids", self.num_grids)
        self.grid_dimensions = na.ones((self.num_grids,3), 'int32')
        self.grid_left_edge = na.zeros((self.num_grids,3), self.float_type)
        self.grid_right_edge = na.ones((self.num_grids,3), self.float_type)
        self.grid_levels = na.zeros((self.num_grids,1), 'int32')
        self.grid_particle_count = na.zeros((self.num_grids,1), 'int32')

    def _guess_data_style(self, rank, test_grid, test_grid_id):
        if self.data_style is not None: return
        if test_grid[0] != os.path.sep:
            test_grid = os.path.join(self.directory, test_grid)
        if not os.path.exists(test_grid):
            test_grid = os.path.join(self.directory,
                                    os.path.basename(test_grid))
            mylog.debug("Your data uses the annoying hardcoded path.")
            self._strip_path = True
        try:
            a = SD.SD(test_grid)
            self.data_style = 'enzo_hdf4'
            mylog.debug("Detected HDF4")
        except:
            list_of_sets = HDF5LightReader.ReadListOfDatasets(test_grid, "/")
            if len(list_of_sets) == 0 and rank == 3:
                mylog.debug("Detected packed HDF5")
                self.data_style = 'enzo_packed_3d'
            elif len(list_of_sets) > 0 and rank == 3:
                mylog.debug("Detected unpacked HDF5")
                self.data_style = 'enzo_hdf5'
            elif len(list_of_sets) == 0 and rank == 2:
                mylog.debug("Detect packed 2D")
                self.data_style = 'enzo_packed_2d'
            elif len(list_of_sets) == 0 and rank == 1:
                mylog.debug("Detect packed 1D")
                self.data_style = 'enzo_packed_1d'
            else:
                raise TypeError

    # Sets are sorted, so that won't work!
    __parse_tokens = ("GridStartIndex", "GridEndIndex",
            "GridLeftEdge", "GridRightEdge", "FileName")
    def _parse_hierarchy(self):
        def _line_yielder(f):
            for token in self.__parse_tokens:
                line = f.readline()
                while token not in line: line = f.readline()
                yield line.split()[2:]
        t1 = time.time()
        pattern = r"Pointer: Grid\[(\d*)\]->NextGrid(Next|This)Level = (\d*)$"
        patt = re.compile(pattern)
        f = open(self.hierarchy_filename, "rb")
        self.grids = [self.grid(1)]
        self.grids[0].Level = 0
        self.wrefs = [weakref.proxy(self.grids[0])]
        si, ei, LE, RE, fn, np = [], [], [], [], [], []
        all = [si, ei, LE, RE, fn]
        f.readline() # Blank at top
        for grid_id in xrange(self.num_grids):
            for a, v in izip(all, _line_yielder(f)):
                a.append(v)
            line = f.readline()
            np.append("0") # this gets replaced if it finds something...
            while len(line) > 2:
                if line.startswith("Pointer:"):
                    vv = patt.findall(line)[0]
                    self.__pointer_handler(vv)
                    line = f.readline()
                    continue
                params = line.split()
                if "NumberOfParticles" == params[0]:
                    np[-1] = params[2]
                line = f.readline()
        self.grid_dimensions.flat[:] = ei
        self.grid_dimensions -= na.array(si, self.float_type)
        self.grid_dimensions += 1
        self.grid_left_edge.flat[:] = LE
        self.grid_right_edge.flat[:] = RE
        self.grid_particle_count.flat[:] = np
        self.grids = na.array(self.grids, dtype='object')
        self.filenames = fn
        t2 = time.time()

    def __pointer_handler(self, m):
        sgi = int(m[2])-1
        if sgi == -1: return # if it's 0, then we're done with that lineage
        # Okay, so, we have a pointer.  We make a new grid, with an id of the length+1
        # (recall, Enzo grids are 1-indexed)
        self.grids.append(self.grid(len(self.grids)+1))
        # We'll just go ahead and make a weakref to cache
        self.wrefs.append(weakref.proxy(self.grids[-1]))
        second_grid = self.wrefs[sgi] # zero-indexex already
        first_grid = self.wrefs[int(m[0])-1]
        if m[1] == "Next":
            first_grid.Children.append(second_grid)
            second_grid.Parent = first_grid
            second_grid.Level = first_grid.Level + 1
        elif m[1] == "This":
            if first_grid.Parent is not None:
                first_grid.Parent.Children.append(second_grid)
                second_grid.Parent = first_grid.Parent
            second_grid.Level = first_grid.Level

    def _populate_grid_objects(self):
        for g,f in izip(self.grids, self.filenames):
            g._prepare_grid()
            g._setup_dx()
            g.set_filename(f[0])
        del self.filenames # No longer needed.
        self.max_level = self.grid_levels.max()

    def _detect_fields(self):
        field_list = self.get_data("/", "DataFields")
        if field_list is None:
            mylog.info("Gathering a field list (this may take a moment.)")
            field_list = set()
            random_sample = self._generate_random_grids()
            for grid in random_sample:
                if not hasattr(grid, 'filename'): continue
                try:
                    gf = grid.getFields()
                except grid._read_exception:
                    mylog.debug("Grid %s is a bit funky?", grid.id)
                    continue
                mylog.debug("Grid %s has: %s", grid.id, gf)
                field_list = field_list.union(gf)
            field_list = self._join_field_lists(field_list)
            self.save_data(list(field_list),"/","DataFields")
        self.field_list = list(field_list)

    def _setup_unknown_fields(self):
        for field in self.field_list:
            if field in self.parameter_file.field_info: continue
            mylog.info("Adding %s to list of fields", field)
            cf = None
            if self.parameter_file.has_key(field):
                cf = lambda d: d.convert(field)
            add_field(field, lambda a, b: None,
                      convert_function=cf, take_log=False)

    def _setup_derived_fields(self):
        self.derived_field_list = []
        for field in self.parameter_file.field_info:
            try:
                fd = self.parameter_file.field_info[field].get_dependencies(
                            pf = self.parameter_file)
            except:
                continue
            available = na.all([f in self.field_list for f in fd.requested])
            if available: self.derived_field_list.append(field)
        for field in self.field_list:
            if field not in self.derived_field_list:
                self.derived_field_list.append(field)

    def _generate_random_grids(self):
        if self.num_grids > 40:
            starter = na.random.randint(0, 20)
            random_sample = na.mgrid[starter:len(self.grids)-1:20j].astype("int32")
            mylog.debug("Checking grids: %s", random_sample.tolist())
        else:
            random_sample = na.mgrid[0:max(len(self.grids)-1,1)].astype("int32")
        return self.grids[(random_sample,)]

class EnzoHierarchyInMemory(EnzoHierarchy):
    _data_style = 'enzo_inline'
    def _obtain_enzo(self):
        import enzo; return enzo

    def __init__(self, pf, data_style = None):
        self.parameter_file = weakref.proxy(pf) # for _obtain_enzo
        if data_style is None: data_style = self._data_style
        enzo = self._obtain_enzo()
        self.float_type = 'float64'
        self.data_style = data_style # Mandated
        self.directory = os.getcwd()
        self.num_grids = enzo.hierarchy_information["GridDimensions"].shape[0]
        self._setup_data_queue()
        AMRHierarchy.__init__(self, pf)

    def _initialize_data_file(self):
        pass

    def _join_field_lists(self, field_list):
        from mpi4py import MPI
        MPI.COMM_WORLD.Barrier()
        data = list(field_list)
        if MPI.COMM_WORLD.rank == 0:
            for i in range(1, MPI.COMM_WORLD.size):
                data += MPI.COMM_WORLD.recv(source=i, tag=0)
            data = list(set(data))
        else:
            MPI.COMM_WORLD.send(data, dest=0, tag=0)
        MPI.COMM_WORLD.Barrier()
        return MPI.COMM_WORLD.bcast(data, root=0)

    def _populate_hierarchy(self):
        self._copy_hierarchy_structure()
        enzo = self._obtain_enzo()
        mylog.debug("Copying reverse tree")
        self.gridReverseTree = enzo.hierarchy_information["GridParentIDs"].ravel().tolist()
        # Initial setup:
        mylog.debug("Reconstructing parent-child relationships")
        #self.gridTree = [ [] for i in range(self.num_grids) ]
        for id,pid in enumerate(self.gridReverseTree):
            if pid > 0:
                self.gridTree[pid-1].append(
                    weakref.proxy(self.grids[id]))
            else:
                self.gridReverseTree[id] = None
        self.max_level = self.gridLevels.max()
        self.maxLevel = self.max_level
        mylog.debug("Preparing grids")
        for i, grid in enumerate(self.grids):
            if (i%1e4) == 0: mylog.debug("Prepared % 7i / % 7i grids", i, self.num_grids)
            grid.filename = None
            grid._prepare_grid()
            grid.proc_num = self.gridProcs[i,0]
        self._setup_grid_dxs()
        mylog.debug("Prepared")
        self._setup_field_lists()
        self.levelIndices = {}
        self.levelNum = {}
        ad = self.gridEndIndices - self.gridStartIndices + 1
        mylog.debug("Hierarchy fully populated.")

    def _copy_hierarchy_structure(self):
        enzo = self._obtain_enzo()
        self.gridDimensions[:] = enzo.hierarchy_information["GridDimensions"][:]
        self.gridStartIndices[:] = enzo.hierarchy_information["GridStartIndices"][:]
        self.gridEndIndices[:] = enzo.hierarchy_information["GridEndIndices"][:]
        self.gridLeftEdge[:] = enzo.hierarchy_information["GridLeftEdge"][:]
        self.gridRightEdge[:] = enzo.hierarchy_information["GridRightEdge"][:]
        self.gridLevels[:] = enzo.hierarchy_information["GridLevels"][:]
        self.gridTimes[:] = enzo.hierarchy_information["GridTimes"][:]
        self.gridProcs = enzo.hierarchy_information["GridProcs"].copy()
        self.gridNumberOfParticles[:] = enzo.hierarchy_information["GridNumberOfParticles"][:]

    def _generate_random_grids(self):
        my_proc = ytcfg.getint("yt","__parallel_rank")
        gg = self.grids[self.gridProcs[:,0] == my_proc]
        if len(gg) > 40:
            starter = na.random.randint(0, 20)
            random_sample = na.mgrid[starter:len(gg)-1:20j].astype("int32")
            mylog.debug("Checking grids: %s", random_sample.tolist())
        else:
            random_sample = na.mgrid[0:max(len(gg)-1,1)].astype("int32")
        return gg[(random_sample,)]

    def save_data(self, *args, **kwargs):
        pass

class EnzoHierarchy1D(EnzoHierarchy):
    def __init__(self, *args, **kwargs):
        EnzoHierarchy.__init__(self, *args, **kwargs)
        self.gridLeftEdge[:,1:3] = 0.0
        self.gridRightEdge[:,1:3] = 1.0
        self.gridDimensions[:,1:3] = 1.0
        self.gridDys[:,0] = 1.0
        self.gridDzs[:,0] = 1.0
        for g in self.grids:
            g._prepare_grid()
            g._setup_dx()

class EnzoHierarchy2D(EnzoHierarchy):
    def __init__(self, *args, **kwargs):
        EnzoHierarchy.__init__(self, *args, **kwargs)
        self.gridLeftEdge[:,2] = 0.0
        self.gridRightEdge[:,2] = 1.0
        self.gridDimensions[:,2] = 1.0
        self.gridDzs[:,0] = 1.0
        for g in self.grids:
            g._prepare_grid()
            g._setup_dx()

scanf_regex = {}
scanf_regex['e'] = r"[-+]?\d+\.?\d*?|\.\d+[eE][-+]?\d+?"
scanf_regex['g'] = scanf_regex['e']
scanf_regex['f'] = scanf_regex['e']
scanf_regex['F'] = scanf_regex['e']
#scanf_regex['g'] = r"[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?"
#scanf_regex['f'] = r"[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?"
#scanf_regex['F'] = r"[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?"
scanf_regex['i'] = r"[-+]?(0[xX][\dA-Fa-f]+|0[0-7]*|\d+)"
scanf_regex['d'] = r"[-+]?\d+"
scanf_regex['s'] = r"\S+"

def constructRegularExpressions(param, toReadTypes):
    re_e=r"[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?"
    re_i=r"[-+]?(0[xX][\dA-Fa-f]+|0[0-7]*|\d+)"
    rs = "^%s\s*=\s*" % (param)
    for t in toReadTypes:
        rs += "(%s)\s*" % (scanf_regex[t])
    rs +="$"
    return re.compile(rs,re.M)

# These next two functions are taken from
# http://www.reddit.com/r/Python/comments/6hj75/reverse_file_iterator/c03vms4
# Credit goes to "Brian" on Reddit

def rblocks(f, blocksize=4096):
    """Read file as series of blocks from end of file to start.

    The data itself is in normal order, only the order of the blocks is reversed.
    ie. "hello world" -> ["ld","wor", "lo ", "hel"]
    Note that the file must be opened in binary mode.
    """
    if 'b' not in f.mode.lower():
        raise Exception("File must be opened using binary mode.")
    size = os.stat(f.name).st_size
    fullblocks, lastblock = divmod(size, blocksize)

    # The first(end of file) block will be short, since this leaves 
    # the rest aligned on a blocksize boundary.  This may be more 
    # efficient than having the last (first in file) block be short
    f.seek(-lastblock,2)
    yield f.read(lastblock)

    for i in range(fullblocks-1,-1, -1):
        f.seek(i * blocksize)
        yield f.read(blocksize)

def rlines(f, keepends=False):
    """Iterate through the lines of a file in reverse order.

    If keepends is true, line endings are kept as part of the line.
    """
    buf = ''
    for block in rblocks(f):
        buf = block + buf
        lines = buf.splitlines(keepends)
        # Return all lines except the first (since may be partial)
        if lines:
            lines.reverse()
            buf = lines.pop() # Last line becomes end of new first line.
            for line in lines:
                yield line
    yield buf  # First line.

class OrionHierarchy(AMRHierarchy):
    def __init__(self, pf, data_style='orion_native'):
        self.field_info = OrionFieldContainer()
        self.field_indexes = {}
        self.parameter_file = weakref.proxy(pf)
        header_filename = os.path.join(pf.fullplotdir,'Header')
        self.directory = pf.fullpath
        self.data_style = data_style
        self._setup_classes()
        self.readGlobalHeader(header_filename,self.parameter_file.paranoid_read) # also sets up the grid objects
        self.__cache_endianness(self.levels[-1].grids[-1])
        AMRHierarchy.__init__(self,pf)
        self._setup_data_queue()
        self._setup_field_list()

    def readGlobalHeader(self,filename,paranoid_read):
        """
        read the global header file for an Orion plotfile output.
        """
        counter = 0
        header_file = open(filename,'r')
        self.__global_header_lines = header_file.readlines()

        # parse the file
        self.orion_version = self.__global_header_lines[0].rstrip()
        self.n_fields      = int(self.__global_header_lines[1])

        counter = self.n_fields+2
        self.field_list = []
        for i,line in enumerate(self.__global_header_lines[2:counter]):
            self.field_list.append(line.rstrip())

        # this is unused...eliminate it?
        #for f in self.field_indexes:
        #    self.field_list.append(orion2ytFieldsDict.get(f,f))

        self.dimension = int(self.__global_header_lines[counter])
        if self.dimension != 3:
            raise RunTimeError("Orion must be in 3D to use yt.")
        counter += 1
        self.Time = float(self.__global_header_lines[counter])
        counter += 1
        self.finest_grid_level = int(self.__global_header_lines[counter])
        self.n_levels = self.finest_grid_level + 1
        counter += 1
        # quantities with _unnecessary are also stored in the inputs
        # file and are not needed.  they are read in and stored in
        # case in the future we want to enable a "backwards" way of
        # taking the data out of the Header file and using it to fill
        # in in the case of a missing inputs file
        self.domainLeftEdge_unnecessary = na.array(map(float,self.__global_header_lines[counter].split()))
        counter += 1
        self.domainRightEdge_unnecessary = na.array(map(float,self.__global_header_lines[counter].split()))
        counter += 1
        self.refinementFactor_unnecessary = self.__global_header_lines[counter].split() #na.array(map(int,self.__global_header_lines[counter].split()))
        counter += 1
        self.globalIndexSpace_unnecessary = self.__global_header_lines[counter]
        #domain_re.search(self.__global_header_lines[counter]).groups()
        counter += 1
        self.timestepsPerLevel_unnecessary = self.__global_header_lines[counter]
        counter += 1
        self.dx = na.zeros((self.n_levels,3))
        for i,line in enumerate(self.__global_header_lines[counter:counter+self.n_levels]):
            self.dx[i] = na.array(map(float,line.split()))
        counter += self.n_levels
        self.geometry = int(self.__global_header_lines[counter])
        if self.geometry != 0:
            raise RunTimeError("yt only supports cartesian coordinates.")
        counter += 1

        # this is just to debug. eventually it should go away.
        linebreak = int(self.__global_header_lines[counter])
        if linebreak != 0:
            raise RunTimeError("INTERNAL ERROR! This should be a zero.")
        counter += 1

        # each level is one group with ngrids on it. each grid has 3 lines of 2 reals
        self.levels = []
        grid_counter = 0
        file_finder_pattern = r"FabOnDisk: (\w+_D_[0-9]{4}) (\d+)\n"
        re_file_finder = re.compile(file_finder_pattern)
        dim_finder_pattern = r"\(\((\d+,\d+,\d+)\) \((\d+,\d+,\d+)\) \(\d+,\d+,\d+\)\)\n"
        re_dim_finder = re.compile(dim_finder_pattern)
        data_files_pattern = r"Level_[\d]/"
        data_files_finder = re.compile(data_files_pattern)

        for level in range(0,self.n_levels):
            tmp = self.__global_header_lines[counter].split()
            # should this be grid_time or level_time??
            lev,ngrids,grid_time = int(tmp[0]),int(tmp[1]),float(tmp[2])
            counter += 1
            nsteps = int(self.__global_header_lines[counter])
            counter += 1
            self.levels.append(OrionLevel(lev,ngrids))
            # open level header, extract file names and offsets for
            # each grid
            # read slightly out of order here: at the end of the lo,hi
            # pairs for x,y,z is a *list* of files types in the Level
            # directory. each type has Header and a number of data
            # files (one per processor)
            tmp_offset = counter + 3*ngrids
            nfiles = 0
            key_off = 0
            files =   {} # dict(map(lambda a: (a,[]),self.field_list))
            offsets = {} # dict(map(lambda a: (a,[]),self.field_list))
            while nfiles+tmp_offset < len(self.__global_header_lines) and data_files_finder.match(self.__global_header_lines[nfiles+tmp_offset]):
                filen = os.path.join(self.parameter_file.fullplotdir, \
                                     self.__global_header_lines[nfiles+tmp_offset].strip())
                # open each "_H" header file, and get the number of
                # components within it
                level_header_file = open(filen+'_H','r').read()
                start_stop_index = re_dim_finder.findall(level_header_file) # just take the last one
                grid_file_offset = re_file_finder.findall(level_header_file)
                ncomp_this_file = int(level_header_file.split('\n')[2])
                for i in range(ncomp_this_file):
                    key = self.field_list[i+key_off]
                    f,o = zip(*grid_file_offset)
                    files[key] = f
                    offsets[key] = o
                    self.field_indexes[key] = i
                key_off += ncomp_this_file
                nfiles += 1
            # convert dict of lists to list of dicts
            fn = []
            off = []
            lead_path = os.path.join(self.parameter_file.fullplotdir,'Level_%i'%level)
            for i in range(ngrids):
                fi = [os.path.join(lead_path,files[key][i]) for key in self.field_list]
                of = [int(offsets[key][i]) for key in self.field_list]
                fn.append(dict(zip(self.field_list,fi)))
                off.append(dict(zip(self.field_list,of)))

            for grid in range(0,ngrids):
                gfn = fn[grid]  # filename of file containing this grid
                gfo = off[grid] # offset within that file
                xlo,xhi = map(float,self.__global_header_lines[counter].split())
                counter+=1
                ylo,yhi = map(float,self.__global_header_lines[counter].split())
                counter+=1
                zlo,zhi = map(float,self.__global_header_lines[counter].split())
                counter+=1
                lo = na.array([xlo,ylo,zlo])
                hi = na.array([xhi,yhi,zhi])
                dims,start,stop = self.__calculate_grid_dimensions(start_stop_index[grid])
                self.levels[-1].grids.append(self.grid(lo,hi,grid_counter,level,gfn, gfo, dims,start,stop,paranoia=paranoid_read))
                grid_counter += 1 # this is global, and shouldn't be reset
                                  # for each level

            # already read the filenames above...
            counter+=nfiles
            self.num_grids = grid_counter
            self.float_type = 'float64'

        self.maxLevel = self.n_levels - 1 
        self.max_level = self.n_levels - 1
        header_file.close()

    def __cache_endianness(self,test_grid):
        """
        Cache the endianness and bytes perreal of the grids by using a
        test grid and assuming that all grids have the same
        endianness. This is a pretty safe assumption since Orion uses
        one file per processor, and if you're running on a cluster
        with different endian processors, then you're on your own!
        """
        # open the test file & grab the header
        inFile = open(os.path.expanduser(test_grid.filename[self.field_list[0]]),'rb')
        header = inFile.readline()
        inFile.close()
        header.strip()
        
        # parse it. the patter is in OrionDefs.py
        headerRe = re.compile(orion_FAB_header_pattern)
        bytesPerReal,endian,start,stop,centerType,nComponents = headerRe.search(header).groups()
        self._bytesPerReal = int(bytesPerReal)
        if self._bytesPerReal == int(endian[0]):
            dtype = '<'
        elif self._bytesPerReal == int(endian[-1]):
            dtype = '>'
        else:
            raise ValueError("FAB header is neither big nor little endian. Perhaps the file is corrupt?")

        dtype += ('f%i' % self._bytesPerReal) # always a floating point
        self._dtype = dtype

    def __calculate_grid_dimensions(self,start_stop):
        start = na.array(map(int,start_stop[0].split(',')))
        stop = na.array(map(int,start_stop[1].split(',')))
        dimension = stop - start + 1
        return dimension,start,stop
        

    def _initialize_grids(self):
        mylog.debug("Allocating memory for %s grids", self.num_grids)
        self.gridDimensions = na.zeros((self.num_grids,3), 'int32')
        self.gridStartIndices = na.zeros((self.num_grids,3), 'int32')
        self.gridEndIndices = na.zeros((self.num_grids,3), 'int32')
        self.gridTimes = na.zeros((self.num_grids,1), 'float64')
        self.gridNumberOfParticles = na.zeros((self.num_grids,1))
        mylog.debug("Done allocating")
        mylog.debug("Creating grid objects")
        self.grids = na.concatenate([level.grids for level in self.levels])
        self.gridLevels = na.concatenate([level.ngrids*[level.level] for level in self.levels])
        self.gridLevels = self.gridLevels.reshape((self.num_grids,1))
        gridDcs = na.concatenate([level.ngrids*[self.dx[level.level]] for level in self.levels],axis=0)
        self.gridDxs = gridDcs[:,0].reshape((self.num_grids,1))
        self.gridDys = gridDcs[:,1].reshape((self.num_grids,1))
        self.gridDzs = gridDcs[:,2].reshape((self.num_grids,1))
        left_edges = []
        right_edges = []
        for level in self.levels:
            left_edges += [g.LeftEdge for g in level.grids]
            right_edges += [g.RightEdge for g in level.grids]
        self.gridLeftEdge = na.array(left_edges)
        self.gridRightEdge = na.array(right_edges)
        self.gridReverseTree = [] * self.num_grids
        self.gridReverseTree = [ [] for i in range(self.num_grids)]
        self.gridTree = [ [] for i in range(self.num_grids)]
        mylog.debug("Done creating grid objects")

    def _populate_hierarchy(self):
        self.__setup_grid_tree()
        self._setup_grid_corners()
        for i, grid in enumerate(self.grids):
            if (i%1e4) == 0: mylog.debug("Prepared % 7i / % 7i grids", i, self.num_grids)
            grid._prepare_grid()
            grid._setup_dx()

    def __setup_grid_tree(self):
        for i, grid in enumerate(self.grids):
            children = self._get_grid_children(grid)
            for child in children:
                self.gridReverseTree[child.id].append(i)
                self.gridTree[i].append(weakref.proxy(child))

    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        dd["field_indexes"] = self.field_indexes
        AMRHierarchy._setup_classes(self, dd)
        self._add_object_class('grid', "OrionGrid", OrionGridBase, dd)
        self.object_types.sort()

    def _get_grid_children(self, grid):
        mask = na.zeros(self.num_grids, dtype='bool')
        grids, grid_ind = self.get_box_grids(grid.LeftEdge, grid.RightEdge)
        mask[grid_ind] = True
        mask = na.logical_and(mask, (self.gridLevels == (grid.Level+1)).flat)
        return self.grids[mask]

    def _setup_field_list(self):
        self.derived_field_list = []
        for field in self.field_info:
            try:
                fd = self.field_info[field].get_dependencies(pf = self.parameter_file)
            except:
                continue
            available = na.all([f in self.field_list for f in fd.requested])
            if available: self.derived_field_list.append(field)
        for field in self.field_list:
            if field not in self.derived_field_list:
                self.derived_field_list.append(field)


class OrionLevel:
    def __init__(self,level,ngrids):
        self.level = level
        self.ngrids = ngrids
        self.grids = []
    
