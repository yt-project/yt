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

class AMRHierarchy(ObjectFindingMixin, ParallelAnalysisInterface):
    float_type = 'float64'

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

        mylog.debug("Initializing data grid data IO")
        self._setup_data_io()

        mylog.debug("Detecting fields.")
        self._detect_fields()

        mylog.debug("Adding unknown detected fields")
        self._setup_unknown_fields()

        mylog.debug("Setting up derived fields")
        self._setup_derived_fields()

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
        self.num_grids = None

    def _initialize_grid_arrays(self):
        mylog.debug("Allocating arrays for %s grids", self.num_grids)
        self.grid_dimensions = na.ones((self.num_grids,3), 'int32')
        self.grid_left_edge = na.zeros((self.num_grids,3), self.float_type)
        self.grid_right_edge = na.ones((self.num_grids,3), self.float_type)
        self.grid_levels = na.zeros((self.num_grids,1), 'int32')
        self.grid_particle_count = na.zeros((self.num_grids,1), 'int32')

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

    def clear_all_data(self):
        """
        This routine clears all the data currently being held onto by the grids
        and the data io handler.
        """
        for g in self.grids: g.clear_data()
        self.io.queue.clear()

    def _get_data_reader_dict(self):
        dd = { 'pf' : self.parameter_file, # Already weak
               'hierarchy': weakref.proxy(self) }
        return dd

    def _initialize_data_storage(self):
        if not ytcfg.getboolean('lagos','serialize'): return
        if os.path.isfile(os.path.join(self.directory,
                            "%s.yt" % self.pf["CurrentTimeIdentifier"])):
            fn = os.path.join(self.directory,"%s.yt" % self.pf["CurrentTimeIdentifier"])
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

    def _setup_data_io(self):
        self.io = io_registry[self.data_style]()

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
        """
        Save an object (*obj*) to the data_file using the Pickle protocol,
        under the name *name* on the node /Objects.
        """
        s = cPickle.dumps(obj, protocol=-1)
        self.save_data(s, "/Objects", name, force = True)

    def load_object(self, name):
        """
        Load and return and object from the data_file using the Pickle protocol,
        under the name *name* on the node /Objects.
        """
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
        return self._data_file[full_name][:]

    def _close_data_file(self):
        if self._data_file:
            self._data_file.close()
            del self._data_file
            self._data_file = None

    def _deserialize_hierarchy(self, harray):
        # THIS IS BROKEN AND NEEDS TO BE FIXED
        mylog.debug("Cached entry found.")
        self.gridDimensions[:] = harray[:,0:3]
        self.gridStartIndices[:] = harray[:,3:6]
        self.gridEndIndices[:] = harray[:,6:9]
        self.gridLeftEdge[:] = harray[:,9:12]
        self.gridRightEdge[:] = harray[:,12:15]
        self.gridLevels[:] = harray[:,15:16]
        self.gridTimes[:] = harray[:,16:17]
        self.gridNumberOfParticles[:] = harray[:,17:18]

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
    grid = EnzoGrid

    def __init__(self, pf, data_style):
        
        self.data_style = data_style
        self.hierarchy_filename = os.path.abspath(
            "%s.hierarchy" % (pf.parameter_filename))
        harray_fn = self.hierarchy_filename[:-9] + "harrays"
        if os.path.exists(harray_fn):
            harray_fp = h5py.File(harray_fn)
            self.num_grids = harray_fp["/Level"].len()
        elif os.path.getsize(self.hierarchy_filename) == 0:
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

    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        AMRHierarchy._setup_classes(self, dd)
        self.object_types.sort()

    def _count_grids(self):
        test_grid = test_grid_id = None
        self.num_stars = 0
        for line in rlines(open(self.hierarchy_filename, "rb")):
            if line.startswith("BaryonFileName") or \
               line.startswith("FileName "):
                test_grid = line.split("=")[-1].strip().rstrip()
            if line.startswith("NumberOfStarParticles"):
                self.num_stars = int(line.split("=")[-1])
            if line.startswith("Grid "):
                self.num_grids = test_grid_id = int(line.split("=")[-1])
                break
        self._guess_data_style(self.pf["TopGridRank"], test_grid, test_grid_id)

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
    def _parse_hierarchy(self):
        def _next_token_line(token, f):
            line = f.readline()
            while token not in line:
                line = f.readline()
            return line.split()[2:]
        if os.path.exists(self.hierarchy_filename[:-9] + "harrays"):
            if self._parse_binary_hierarchy(): return
        t1 = time.time()
        pattern = r"Pointer: Grid\[(\d*)\]->NextGrid(Next|This)Level = (\d*)$"
        patt = re.compile(pattern)
        f = open(self.hierarchy_filename, "rb")
        self.grids = [self.grid(1, self)]
        self.grids[0].Level = 0
        si, ei, LE, RE, fn, np = [], [], [], [], [], []
        all = [si, ei, LE, RE, fn]
        f.readline() # Blank at top
        for grid_id in xrange(self.num_grids):
            # We will unroll this list
            si.append(_next_token_line("GridStartIndex", f))
            ei.append(_next_token_line("GridEndIndex", f))
            LE.append(_next_token_line("GridLeftEdge", f))
            RE.append(_next_token_line("GridRightEdge", f))
            nb = int(_next_token_line("NumberOfBaryonFields", f)[0])
            fn.append(["-1"])
            if nb > 0: fn[-1] = _next_token_line("BaryonFileName", f)
            np.append(int(_next_token_line("NumberOfParticles", f)[0]))
            if nb == 0 and np[-1] > 0: fn[-1] = _next_token_line("FileName", f)
            line = f.readline()
            while len(line) > 2:
                if line.startswith("Pointer:"):
                    vv = patt.findall(line)[0]
                    self.__pointer_handler(vv)
                    line = f.readline()
                    continue
                params = line.split()
                line = f.readline()
        self._fill_arrays(ei, si, LE, RE, np)
        self.grids = na.array(self.grids, dtype='object')
        self.filenames = fn
        self._store_binary_hierarchy()
        t2 = time.time()

    def _fill_arrays(self, ei, si, LE, RE, np):
        self.grid_dimensions.flat[:] = ei
        self.grid_dimensions -= na.array(si, self.float_type)
        self.grid_dimensions += 1
        self.grid_left_edge.flat[:] = LE
        self.grid_right_edge.flat[:] = RE
        self.grid_particle_count.flat[:] = np

    def __pointer_handler(self, m):
        sgi = int(m[2])-1
        if sgi == -1: return # if it's 0, then we're done with that lineage
        # Okay, so, we have a pointer.  We make a new grid, with an id of the length+1
        # (recall, Enzo grids are 1-indexed)
        self.grids.append(self.grid(len(self.grids)+1, self))
        # We'll just go ahead and make a weakref to cache
        second_grid = self.grids[sgi] # zero-indexed already
        first_grid = self.grids[int(m[0])-1]
        if m[1] == "Next":
            first_grid._children_ids.append(second_grid.id)
            second_grid._parent_id = first_grid.id
            second_grid.Level = first_grid.Level + 1
        elif m[1] == "This":
            if first_grid.Parent is not None:
                first_grid.Parent._children_ids.append(second_grid.id)
                second_grid._parent_id = first_grid._parent_id
            second_grid.Level = first_grid.Level
        self.grid_levels[sgi] = second_grid.Level

    def _parse_binary_hierarchy(self):
        mylog.info("Getting the binary hierarchy")
        f = h5py.File(self.hierarchy_filename[:-9] + "harrays")
        self.grid_dimensions[:] = f["/ActiveDimensions"][:]
        self.grid_left_edge[:] = f["/LeftEdges"][:]
        self.grid_right_edge[:] = f["/RightEdges"][:]
        self.grid_particle_count[:,0] = f["/NumberOfParticles"][:]
        levels = f["/Level"][:]
        parents = f["/ParentIDs"][:]
        procs = f["/Processor"][:]
        grids = []
        self.filenames = []
        grids = [self.grid(gi+1, self) for gi in xrange(self.num_grids)]
        giter = izip(grids, levels, procs, parents)
        bn = "%s.cpu%%04i" % (self.pf)
        pmap = [(bn % P,) for P in xrange(procs.max()+1)]
        for grid,L,P,Pid in giter:
            grid.Level = L
            grid._parent_id = Pid
            if Pid > -1:
                grids[Pid-1]._children_ids.append(grid.id)
            self.filenames.append(pmap[P])
        self.grids = na.array(grids, dtype='object')
        f.close()
        mylog.info("Finished with binary hierarchy reading")
        return True

    @parallel_blocking_call
    def _store_binary_hierarchy(self):
        # We don't do any of the logic here, we just check if the data file
        # is open...
        if self._data_file is None: return
        if self.data_style != "enzo_packed_3d": return
        mylog.info("Storing the binary hierarchy")
        f = h5py.File(self.hierarchy_filename[:-9] + "harrays", "w")
        f.create_dataset("/LeftEdges", data=self.grid_left_edge)
        f.create_dataset("/RightEdges", data=self.grid_right_edge)
        parents, procs, levels = [], [], []
        for i,g in enumerate(self.grids):
            if g.Parent is not None:
                parents.append(g.Parent.id)
            else:
                parents.append(-1)
            procs.append(int(self.filenames[i][0][-4:]))
            levels.append(g.Level)

        parents = na.array(parents, dtype='int64')
        procs = na.array(procs, dtype='int64')
        levels = na.array(levels, dtype='int64')
        f.create_dataset("/ParentIDs", data=parents)
        f.create_dataset("/Processor", data=procs)
        f.create_dataset("/Level", data=levels)

        f.create_dataset("/ActiveDimensions", data=self.grid_dimensions)
        f.create_dataset("/NumberOfParticles", data=self.grid_particle_count[:,0])

        f.close()

    def _rebuild_top_grids(self, level = 0):
        #for level in xrange(self.max_level+1):
        mylog.info("Rebuilding grids on level %s", level)
        cmask = (self.grid_levels.flat == (level + 1))
        cmsum = cmask.sum()
        mask = na.zeros(self.num_grids, dtype='bool')
        for grid in self.select_grids(level):
            mask[:] = 0
            LE = self.grid_left_edge[grid.id - grid._id_offset]
            RE = self.grid_right_edge[grid.id - grid._id_offset]
            grids, grid_i = self.get_box_grids(LE, RE)
            mask[grid_i] = 1
            grid._children_ids = []
            cgrids = self.grids[ ( mask * cmask).astype('bool') ]
            mylog.info("%s: %s / %s", grid, len(cgrids), cmsum)
            for cgrid in cgrids:
                grid._children_ids.append(cgrid.id)
                cgrid._parent_id = grid.id
        mylog.info("Finished rebuilding")

    def _populate_grid_objects(self):
        for g,f in izip(self.grids, self.filenames):
            g._prepare_grid()
            g._setup_dx()
            g.set_filename(f[0])
            #if g.Parent is not None: g._guess_properties_from_parent()
        del self.filenames # No longer needed.
        self.max_level = self.grid_levels.max()

    def _detect_fields(self):
        self.field_list = []
        # Do this only on the root processor to save disk work.
        if self._mpi_get_rank() == 0 or self._mpi_get_rank() == None:
            field_list = self.get_data("/", "DataFields")
            if field_list is None:
                mylog.info("Gathering a field list (this may take a moment.)")
                field_list = set()
                random_sample = self._generate_random_grids()
                for grid in random_sample:
                    if not hasattr(grid, 'filename'): continue
                    try:
                        gf = self.io._read_field_names(grid)
                    except self.io._read_exception:
                        mylog.debug("Grid %s is a bit funky?", grid.id)
                        continue
                    mylog.debug("Grid %s has: %s", grid.id, gf)
                    field_list = field_list.union(gf)
        else:
            field_list = None
        field_list = self._mpi_bcast_pickled(field_list)
        self.save_data(list(field_list),"/","DataFields",passthrough=True)
        self.field_list = list(field_list)

    def _setup_unknown_fields(self):
        for field in self.field_list:
            if field in self.parameter_file.field_info: continue
            mylog.info("Adding %s to list of fields", field)
            cf = None
            if self.parameter_file.has_key(field):
                def external_wrapper(f):
                    def _convert_function(data):
                        return data.convert(f)
                    return _convert_function
                cf = external_wrapper(field)
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

    def find_particles_by_type(self, ptype, max_num=None, additional_fields=None):
        """
        Returns a structure of arrays with all of the particles'
        positions, velocities, masses, types, IDs, and attributes for
        a particle type **ptype** for a maximum of **max_num**
        particles.  If non-default particle fields are used, provide
        them in **additional_fields**.
        """
        # Not sure whether this routine should be in the general HierarchyType.
        if self.gridNumberOfParticles.sum() == 0:
            mylog.info("Data contains no particles.");
            return None
        if additional_fields is None:
            additional_fields = ['metallicity_fraction', 'creation_time',
                                 'dynamical_time']
        pfields = [f for f in self.field_list if f.startswith('particle_')]
        nattr = self.parameter_file['NumberOfParticleAttributes']
        if nattr > 0:
            pfields += additional_fields[:nattr]
        # Find where the particles reside and count them
        if max_num is None: max_num = 1e100
        total = 0
        pstore = []
        for level in range(self.max_level, -1, -1):
            for grid in self.select_grids(level):
                index = na.where(grid['particle_type'] == ptype)[0]
                total += len(index)
                pstore.append(index)
                if total >= max_num: break
            if total >= max_num: break
        result = None
        if total > 0:
            result = {}
            for p in pfields:
                result[p] = na.zeros(total, 'float64')
            # Now we retrieve data for each field
            ig = count = 0
            for level in range(self.max_level, -1, -1):
                for grid in self.select_grids(level):
                    nidx = len(pstore[ig])
                    if nidx > 0:
                        for p in pfields:
                            result[p][count:count+nidx] = grid[p][pstore[ig]]
                        count += nidx
                    ig += 1
                    if count >= total: break
                if count >= total: break
            # Crop data if retrieved more than max_num
            if count > max_num:
                for p in pfields:
                    result[p] = result[p][0:max_num]
        return result


class EnzoHierarchyInMemory(EnzoHierarchy):

    grid = EnzoGridInMemory
    _enzo = None

    @property
    def enzo(self):
        if self._enzo is None:
            import enzo
            self._enzo = enzo
        return self._enzo

    def __init__(self, pf, data_style = None):
        self.data_style = data_style
        self.float_type = 'float64'
        self.parameter_file = weakref.proxy(pf) # for _obtain_enzo
        self.float_type = self.enzo.hierarchy_information["GridLeftEdge"].dtype
        self.directory = os.getcwd()
        AMRHierarchy.__init__(self, pf, data_style)

    def _initialize_data_storage(self):
        pass

    def _count_grids(self):
        self.num_grids = self.enzo.hierarchy_information["GridDimensions"].shape[0]

    def _parse_hierarchy(self):
        self._copy_hierarchy_structure()
        mylog.debug("Copying reverse tree")
        reverse_tree = self.enzo.hierarchy_information["GridParentIDs"].ravel().tolist()
        # Initial setup:
        mylog.debug("Reconstructing parent-child relationships")
        self.grids = []
        # We enumerate, so it's 0-indexed id and 1-indexed pid
        self.filenames = ["-1"] * self.num_grids
        for id,pid in enumerate(reverse_tree):
            self.grids.append(self.grid(id+1, self))
            self.grids[-1].Level = self.grid_levels[id]
            if pid > 0:
                self.grids[-1]._parent_id = pid
                self.grids[pid-1]._children_ids.append(self.grids[-1].id)
        self.max_level = self.grid_levels.max()
        mylog.debug("Preparing grids")
        for i, grid in enumerate(self.grids):
            if (i%1e4) == 0: mylog.debug("Prepared % 7i / % 7i grids", i, self.num_grids)
            grid.filename = None
            grid._prepare_grid()
            grid.proc_num = self.grid_procs[i,0]
        self.grids = na.array(self.grids, dtype='object')
        mylog.debug("Prepared")

    def _initialize_grid_arrays(self):
        EnzoHierarchy._initialize_grid_arrays(self)
        self.grid_procs = na.zeros((self.num_grids,1),'int32')

    def _copy_hierarchy_structure(self):
        # Dimensions are important!
        self.grid_dimensions[:] = self.enzo.hierarchy_information["GridEndIndices"][:]
        self.grid_dimensions -= self.enzo.hierarchy_information["GridStartIndices"][:]
        self.grid_dimensions += 1
        self.grid_left_edge[:] = self.enzo.hierarchy_information["GridLeftEdge"][:]
        self.grid_right_edge[:] = self.enzo.hierarchy_information["GridRightEdge"][:]
        self.grid_levels[:] = self.enzo.hierarchy_information["GridLevels"][:]
        self.grid_procs = self.enzo.hierarchy_information["GridProcs"].copy()
        self.grid_particle_count[:] = self.enzo.hierarchy_information["GridNumberOfParticles"][:]

    def save_data(self, *args, **kwargs):
        pass

class EnzoHierarchy1D(EnzoHierarchy):

    def _fill_arrays(self, ei, si, LE, RE, np):
        self.grid_dimensions[:,:1] = ei
        self.grid_dimensions[:,:1] -= na.array(si, self.float_type)
        self.grid_dimensions += 1
        self.grid_left_edge[:,:2] = LE
        self.grid_right_edge[:,:2] = RE
        self.grid_particle_count.flat[:] = np
        self.grid_left_edge[:,1:2] = 0.0
        self.grid_right_edge[:,1:2] = 1.0
        self.grid_dimensions[:,1:2] = 1

class EnzoHierarchy2D(EnzoHierarchy):

    def _fill_arrays(self, ei, si, LE, RE, np):
        self.grid_dimensions[:,:2] = ei
        self.grid_dimensions[:,:2] -= na.array(si, self.float_type)
        self.grid_dimensions += 1
        self.grid_left_edge[:,:2] = LE
        self.grid_right_edge[:,:2] = RE
        self.grid_particle_count.flat[:] = np
        self.grid_left_edge[:,2] = 0.0
        self.grid_right_edge[:,2] = 1.0
        self.grid_dimensions[:,2] = 1

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
    grid = OrionGrid
    def __init__(self, pf, data_style='orion_native'):
        self.field_info = OrionFieldContainer()
        self.field_indexes = {}
        self.parameter_file = weakref.proxy(pf)
        header_filename = os.path.join(pf.fullplotdir,'Header')
        self.directory = pf.fullpath
        self.data_style = data_style
        #self._setup_classes()

        self.readGlobalHeader(header_filename,self.parameter_file.paranoid_read) # also sets up the grid objects
        self.__cache_endianness(self.levels[-1].grids[-1])
        AMRHierarchy.__init__(self,pf, self.data_style)
        self._setup_data_io()
        self._setup_field_list()
        self._populate_hierarchy()
        
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
                self.levels[-1].grids.append(self.grid(lo,hi,grid_counter,level,gfn, gfo, dims,start,stop,paranoia=paranoid_read,hierarchy=self))
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
        
    def _populate_grid_objects(self):
        mylog.debug("Creating grid objects")
        self.grids = na.concatenate([level.grids for level in self.levels])
        self.grid_levels = na.concatenate([level.ngrids*[level.level] for level in self.levels])
        self.grid_levels = self.grid_levels.reshape((self.num_grids,1))
        grid_dcs = na.concatenate([level.ngrids*[self.dx[level.level]] for level in self.levels],axis=0)
        self.grid_dxs = grid_dcs[:,0].reshape((self.num_grids,1))
        self.grid_dys = grid_dcs[:,1].reshape((self.num_grids,1))
        self.grid_dzs = grid_dcs[:,2].reshape((self.num_grids,1))
        left_edges = []
        right_edges = []
        for level in self.levels:
            left_edges += [g.LeftEdge for g in level.grids]
            right_edges += [g.RightEdge for g in level.grids]
        self.grid_left_edge = na.array(left_edges)
        self.grid_right_edge = na.array(right_edges)
        self.gridReverseTree = [] * self.num_grids
        self.gridReverseTree = [ [] for i in range(self.num_grids)]
        self.gridTree = [ [] for i in range(self.num_grids)]
        mylog.debug("Done creating grid objects")

    def _populate_hierarchy(self):
        self.__setup_grid_tree()
        #self._setup_grid_corners()
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
        #self._add_object_class('grid', "OrionGrid", OrionGridBase, dd)
        self.object_types.sort()

    def _get_grid_children(self, grid):
        mask = na.zeros(self.num_grids, dtype='bool')
        grids, grid_ind = self.get_box_grids(grid.LeftEdge, grid.RightEdge)
        mask[grid_ind] = True
        mask = na.logical_and(mask, (self.grid_levels == (grid.Level+1)).flat)
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

    def _count_grids(self):
        """this is already provided in 

        """
        pass

    def _initialize_grid_arrays(self):
        mylog.debug("Allocating arrays for %s grids", self.num_grids)
        self.grid_dimensions = na.ones((self.num_grids,3), 'int32')
        self.grid_left_edge = na.zeros((self.num_grids,3), self.float_type)
        self.grid_right_edge = na.ones((self.num_grids,3), self.float_type)
        self.grid_levels = na.zeros((self.num_grids,1), 'int32')
        self.grid_particle_count = na.zeros((self.num_grids,1), 'int32')

    def _parse_hierarchy(self):
        pass
    
    def _detect_fields(self):
        pass

    def _setup_unknown_fields(self):
        pass

    def _setup_derived_fields(self):
        pass

    def _initialize_state_variables(self):
        """override to not re-initialize num_grids in AMRHierarchy.__init__

        """
        self._parallel_locking = False
        self._data_file = None
        self._data_mode = None
        self._max_locations = {}
    
class OrionLevel:
    def __init__(self,level,ngrids):
        self.level = level
        self.ngrids = ngrids
        self.grids = []
    
