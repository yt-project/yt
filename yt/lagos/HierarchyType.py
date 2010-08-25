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
import string, re, gc, time, cPickle, pdb
from itertools import chain, izip
try:
    import yt.ramses_reader as ramses_reader
except ImportError:
    mylog.warning("Ramses Reader not imported")

def num_deep_inc(f):
    def wrap(self, *args, **kwargs):
        self.num_deep += 1
        rv = f(self, *args, **kwargs)
        self.num_deep -= 1
        return rv
    return wrap

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

    def get_levels(self):
        for level in range(self.max_level+1):
            yield self.select_grids(level)

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
        fn = self.pf.storage_filename
        if fn is None:
            if os.path.isfile(os.path.join(self.directory,
                                "%s.yt" % self.pf["CurrentTimeIdentifier"])):
                fn = os.path.join(self.directory,"%s.yt" % self.pf["CurrentTimeIdentifier"])
            else:
                fn = os.path.join(self.directory,
                        "%s.yt" % self.parameter_file.basename)
        dir_to_check = os.path.dirname(fn)
        # We have four options:
        #    Writeable, does not exist      : create, open as append
        #    Writeable, does exist          : open as append
        #    Not writeable, does not exist  : do not attempt to open
        #    Not writeable, does exist      : open as read-only
        exists = os.path.isfile(fn)
        if not exists:
            writeable = os.access(dir_to_check, os.W_OK)
        else:
            writeable = os.access(fn, os.W_OK)
        writeable = writeable and not ytcfg.getboolean('lagos','onlydeserialize')
        # We now have our conditional stuff
        self._barrier()
        if not writeable and not exists: return
        if writeable:
            self._data_mode = 'a'
            if not exists: self.__create_data_file(fn)
        else:
            self._data_mode = 'r'

        self.__data_filename = fn
        self._data_file = h5py.File(fn, self._data_mode)

    def __create_data_file(self, fn):
        # Note that this used to be parallel_root_only; it no longer is,
        # because we have better logic to decide who owns the file.
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
        try:
            return self._data_file[full_name][:]
        except TypeError:
            return self._data_file[full_name]

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
        dims = []
        for level in self.levels:
            left_edges += [g.LeftEdge for g in level.grids]
            right_edges += [g.RightEdge for g in level.grids]
            dims += [g.ActiveDimensions for g in level.grids]
        self.grid_left_edge = na.array(left_edges)
        self.grid_right_edge = na.array(right_edges)
        self.grid_dimensions = na.array(dims)
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
    

class GadgetHierarchy(AMRHierarchy):
    grid = GadgetGrid

    def __init__(self, pf, data_style='gadget_hdf5'):
        self.field_info = GadgetFieldContainer()
        self.directory = os.path.dirname(pf.parameter_filename)
        self.data_style = data_style
        self._handle = h5py.File(pf.parameter_filename)
        AMRHierarchy.__init__(self, pf, data_style)
        self._handle.close()

    def _initialize_data_storage(self):
        pass

    def _detect_fields(self):
        #example string:
        #"(S'VEL'\np1\nS'ID'\np2\nS'MASS'\np3\ntp4\n."
        #fields are surrounded with '
        fields_string=self._handle['root'].attrs['fieldnames']
        #splits=fields_string.split("'")
        #pick out the odd fields
        #fields= [splits[j] for j in range(1,len(splits),2)]
        self.field_list = cPickle.loads(fields_string)
    
    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        AMRHierarchy._setup_classes(self, dd)
        self.object_types.sort()

    def _count_grids(self):
        fh = self._handle #shortcut
        #nodes in the hdf5 file are the same as grids
        #in yt
        #the num of levels and total nodes is already saved
        self._levels   = self.pf._get_param('maxlevel')
        self.num_grids = self.pf._get_param('numnodes')
        
    def _parse_hierarchy(self):
        #for every box, define a self.grid(level,edge1,edge2) 
        #with particle counts, dimensions
        f = self._handle #shortcut
        
        root = f['root']
        grids,numnodes = self._walk_nodes(None,root,[])
        dims = [self.pf.max_grid_size for grid in grids]
        LE = [grid.LeftEdge for grid in grids]
        RE = [grid.RightEdge for grid in grids]
        levels = [grid.Level for grid in grids]
        counts = [(grid.N if grid.IsLeaf else 0) for grid in grids]
        self.grids = na.array(grids,dtype='object')
        self.grid_dimensions[:] = na.array(dims, dtype='int64')
        self.grid_left_edge[:] = na.array(LE, dtype='float64')
        self.grid_right_edge[:] = na.array(RE, dtype='float64')
        self.grid_levels.flat[:] = na.array(levels, dtype='int32')
        self.grid_particle_count.flat[:] = na.array(counts, dtype='int32')
            
    def _walk_nodes(self,parent,node,grids,idx=0):
        pi = cPickle.loads
        loc = node.attrs['h5address']
        
        kwargs = {}
        kwargs['Address'] = loc
        kwargs['Parent'] = parent
        kwargs['Axis']  = self.pf._get_param('divideaxis',location=loc)
        kwargs['Level']  = self.pf._get_param('level',location=loc)
        kwargs['LeftEdge'] = self.pf._get_param('leftedge',location=loc) 
        kwargs['RightEdge'] = self.pf._get_param('rightedge',location=loc)
        kwargs['IsLeaf'] = self.pf._get_param('isleaf',location=loc)
        kwargs['N'] = self.pf._get_param('n',location=loc)
        kwargs['NumberOfParticles'] = self.pf._get_param('n',location=loc)
        dx = self.pf._get_param('dx',location=loc)
        dy = self.pf._get_param('dy',location=loc)
        dz = self.pf._get_param('dz',location=loc)
        divdims = na.array([1,1,1])
        if not kwargs['IsLeaf']: 
            divdims[kwargs['Axis']] = 2
        kwargs['ActiveDimensions'] = divdims
        #Active dimensions:
        #This is the number of childnodes, along with dimensiolaity
        #ie, binary tree can be (2,1,1) but octree is (2,2,2)
        
        idx+=1
        #pdb.set_trace()
        children = []
        if not kwargs['IsLeaf']:
            for child in node.values():
                children,idx=self._walk_nodes(node,child,children,idx=idx)
        
        kwargs['Children'] = children
        grid = self.grid(idx,self.pf.parameter_filename,self,**kwargs)
        grids += children
        grids += [grid,]
        return grids,idx

    def _populate_grid_objects(self):
        for g in self.grids:
            g._prepare_grid()
        self.max_level = cPickle.loads(self._handle['root'].attrs['maxlevel'])
    
    def _setup_unknown_fields(self):
        pass

    def _setup_derived_fields(self):
        self.derived_field_list = []

    def _get_grid_children(self, grid):
        #given a grid, use it's address to find subchildren
        pass

class GadgetHierarchyOld(AMRHierarchy):
    #Kept here to compare for the time being
    grid = GadgetGrid

    def __init__(self, pf, data_style):
        self.directory = pf.fullpath
        self.data_style = data_style
        AMRHierarchy.__init__(self, pf, data_style)

    def _count_grids(self):
        # We actually construct our octree here!
        # ...but we do read in our particles, it seems.
        LE = na.zeros(3, dtype='float64')
        RE = na.ones(3, dtype='float64')
        base_grid = ProtoGadgetGrid(0, LE, RE, self.pf.particles)
        self.proto_grids = base_grid.refine(8)
        self.num_grids = len(self.proto_grids)
        self.max_level = max( (g.level for g in self.proto_grids) )

    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        AMRHierarchy._setup_classes(self, dd)
        self.object_types.sort()

    def _parse_hierarchy(self):
        grids = []
        # We need to fill in dims, LE, RE, level, count
        dims, LE, RE, levels, counts = [], [], [], [], []
        self.proto_grids.sort(key = lambda a: a.level)
        for i, pg in enumerate(self.proto_grids):
            g = self.grid(i, self, pg)
            pg.real_grid = g
            grids.append(g)
            dims.append(g.ActiveDimensions)
            LE.append(g.LeftEdge)
            RE.append(g.RightEdge)
            levels.append(g.Level)
            counts.append(g.NumberOfParticles)
        del self.proto_grids
        self.grids = na.array(grids, dtype='object')
        self.grid_dimensions[:] = na.array(dims, dtype='int64')
        self.grid_left_edge[:] = na.array(LE, dtype='float64')
        self.grid_right_edge[:] = na.array(RE, dtype='float64')
        self.grid_levels.flat[:] = na.array(levels, dtype='int32')
        self.grid_particle_count.flat[:] = na.array(counts, dtype='int32')

    def _populate_grid_objects(self):
        # We don't need to do anything here
        for g in self.grids: g._setup_dx()

    def _detect_fields(self):
        self.field_list = ['particle_position_%s' % ax for ax in 'xyz']

    def _setup_unknown_fields(self):
        pass

    def _setup_derived_fields(self):
        self.derived_field_list = []

class ChomboHierarchy(AMRHierarchy):

    grid = ChomboGrid
    
    def __init__(self,pf,data_style='chombo_hdf5'):
        self.data_style = data_style
        self.field_info = ChomboFieldContainer()
        self.field_indexes = {}
        self.parameter_file = weakref.proxy(pf)
        # for now, the hierarchy file is the parameter file!
        self.hierarchy_filename = self.parameter_file.parameter_filename
        self.directory = os.path.dirname(self.hierarchy_filename)
        self._fhandle = h5py.File(self.hierarchy_filename)

        self.float_type = self._fhandle['/level_0']['data:datatype=0'].dtype.name
        self._levels = self._fhandle.listnames()[1:]
        AMRHierarchy.__init__(self,pf,data_style)

        self._fhandle.close()

    def _initialize_data_storage(self):
        pass

    def _detect_fields(self):
        ncomp = int(self._fhandle['/'].attrs['num_components'])
        self.field_list = [c[1] for c in self._fhandle['/'].attrs.listitems()[-ncomp:]]
    
    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        AMRHierarchy._setup_classes(self, dd)
        self.object_types.sort()

    def _count_grids(self):
        self.num_grids = 0
        for lev in self._levels:
            self.num_grids += self._fhandle[lev]['Processors'].len()
        
    def _parse_hierarchy(self):
        f = self._fhandle # shortcut
        
        # this relies on the first Group in the H5 file being
        # 'Chombo_global'
        levels = f.listnames()[1:]
        self.grids = []
        i = 0
        for lev in levels:
            level_number = int(re.match('level_(\d+)',lev).groups()[0])
            boxes = f[lev]['boxes'].value
            dx = f[lev].attrs['dx']
            for level_id, box in enumerate(boxes):
                si = na.array([box['lo_%s' % ax] for ax in 'ijk'])
                ei = na.array([box['hi_%s' % ax] for ax in 'ijk'])
                pg = self.grid(len(self.grids),self,level=level_number,
                               start = si, stop = ei)
                self.grids.append(pg)
                self.grids[-1]._level_id = level_id
                self.grid_left_edge[i] = dx*si.astype(self.float_type)
                self.grid_right_edge[i] = dx*(ei.astype(self.float_type) + 1)
                self.grid_particle_count[i] = 0
                self.grid_dimensions[i] = ei - si + 1
                i += 1
        self.grids = na.array(self.grids, dtype='object')

    def _populate_grid_objects(self):
        for g in self.grids:
            g._prepare_grid()
            g._setup_dx()

        for g in self.grids:
            g.Children = self._get_grid_children(g)
            for g1 in g.Children:
                g1.Parent.append(g)
        self.max_level = self.grid_levels.max()

    def _setup_unknown_fields(self):
        pass

    def _setup_derived_fields(self):
        self.derived_field_list = []

    def _get_grid_children(self, grid):
        mask = na.zeros(self.num_grids, dtype='bool')
        grids, grid_ind = self.get_box_grids(grid.LeftEdge, grid.RightEdge)
        mask[grid_ind] = True
        return [g for g in self.grids[mask] if g.Level == grid.Level + 1]

class TigerHierarchy(AMRHierarchy):

    grid = TigerGrid

    def __init__(self, pf, data_style):
        self.directory = pf.fullpath
        self.data_style = data_style
        AMRHierarchy.__init__(self, pf, data_style)

    def _count_grids(self):
        # Tiger is unigrid
        self.ngdims = [i/j for i,j in
                izip(self.pf.root_size, self.pf.max_grid_size)]
        self.num_grids = na.prod(self.ngdims)
        self.max_level = 0

    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        AMRHierarchy._setup_classes(self, dd)
        self.object_types.sort()

    def _parse_hierarchy(self):
        grids = []
        # We need to fill in dims, LE, RE, level, count
        dims, LE, RE, levels, counts = [], [], [], [], []
        DLE = self.pf["DomainLeftEdge"]
        DRE = self.pf["DomainRightEdge"] 
        DW = DRE - DLE
        gds = DW / self.ngdims
        rd = [self.pf.root_size[i]-self.pf.max_grid_size[i] for i in range(3)]
        glx, gly, glz = na.mgrid[DLE[0]:DRE[0]-gds[0]:self.ngdims[0]*1j,
                                 DLE[1]:DRE[1]-gds[1]:self.ngdims[1]*1j,
                                 DLE[2]:DRE[2]-gds[2]:self.ngdims[2]*1j]
        gdx, gdy, gdz = na.mgrid[0:rd[0]:self.ngdims[0]*1j,
                                 0:rd[1]:self.ngdims[1]*1j,
                                 0:rd[2]:self.ngdims[2]*1j]
        LE, RE, levels, counts = [], [], [], []
        i = 0
        for glei, gldi in izip(izip(glx.flat, gly.flat, glz.flat),
                               izip(gdx.flat, gdy.flat, gdz.flat)):
            gld = na.array(gldi)
            gle = na.array(glei)
            gre = gle + gds
            g = self.grid(i, self, gle, gre, gld, gld+self.pf.max_grid_size)
            grids.append(g)
            dims.append(self.pf.max_grid_size)
            LE.append(g.LeftEdge)
            RE.append(g.RightEdge)
            levels.append(g.Level)
            counts.append(g.NumberOfParticles)
            i += 1
        self.grids = na.array(grids, dtype='object')
        self.grid_dimensions[:] = na.array(dims, dtype='int64')
        self.grid_left_edge[:] = na.array(LE, dtype='float64')
        self.grid_right_edge[:] = na.array(RE, dtype='float64')
        self.grid_levels.flat[:] = na.array(levels, dtype='int32')
        self.grid_particle_count.flat[:] = na.array(counts, dtype='int32')

    def _populate_grid_objects(self):
        # We don't need to do anything here
        for g in self.grids: g._setup_dx()

    def _detect_fields(self):
        self.file_mapping = {"Density" : "rhob",
                             "Temperature" : "temp"}

    @property
    def field_list(self):
        return self.file_mapping.keys()

    def _setup_unknown_fields(self):
        for field in self.field_list:
            add_tiger_field(field, lambda a, b: None)

    def _setup_derived_fields(self):
        self.derived_field_list = []

class FLASHHierarchy(AMRHierarchy):

    grid = FLASHGrid
    _handle = None
    
    def __init__(self,pf,data_style='chombo_hdf5'):
        self.data_style = data_style
        self.field_info = FLASHFieldContainer()
        self.field_indexes = {}
        self.parameter_file = weakref.proxy(pf)
        # for now, the hierarchy file is the parameter file!
        self.hierarchy_filename = self.parameter_file.parameter_filename
        self.directory = os.path.dirname(self.hierarchy_filename)
        self._handle = h5py.File(self.hierarchy_filename)

        self.float_type = na.float64
        AMRHierarchy.__init__(self,pf,data_style)

        self._handle.close()
        self._handle = None

    def _initialize_data_storage(self):
        pass

    def _detect_fields(self):
        ncomp = self._handle["/unknown names"].shape[0]
        self.field_list = [s.strip() for s in self._handle["/unknown names"][:].flat]
    
    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        AMRHierarchy._setup_classes(self, dd)
        self.object_types.sort()

    def _count_grids(self):
        try:
            self.num_grids = self.parameter_file._find_parameter(
                "integer", "globalnumblocks", True, self._handle)
        except KeyError:
            self.num_grids = self._handle["/simulation parameters"][0][0]
        
    def _parse_hierarchy(self):
        f = self._handle # shortcut
        pf = self.parameter_file # shortcut
        
        self.grid_left_edge[:] = f["/bounding box"][:,:,0]
        self.grid_right_edge[:] = f["/bounding box"][:,:,1]
        # Move this to the parameter file
        try:
            nxb = pf._find_parameter("integer", "nxb", True, f)
            nyb = pf._find_parameter("integer", "nyb", True, f)
            nzb = pf._find_parameter("integer", "nzb", True, f)
        except KeyError:
            nxb, nyb, nzb = [int(f["/simulation parameters"]['n%sb' % ax])
                              for ax in 'xyz']
        self.grid_dimensions[:] *= (nxb, nyb, nzb)
        # particle count will need to be fixed somehow:
        #   by getting access to the particle file we can get the number of
        #   particles in each brick.  but how do we handle accessing the
        #   particle file?

        # This will become redundant, as _prepare_grid will reset it to its
        # current value.  Note that FLASH uses 1-based indexing for refinement
        # levels, but we do not, so we reduce the level by 1.
        self.grid_levels.flat[:] = f["/refine level"][:][:] - 1
        g = [self.grid(i+1, self, self.grid_levels[i,0])
                for i in xrange(self.num_grids)]
        self.grids = na.array(g, dtype='object')

    def _populate_grid_objects(self):
        # We only handle 3D data, so offset is 7 (nfaces+1)
        
        offset = 7
        ii = na.argsort(self.grid_levels.flat)
        gid = self._handle["/gid"][:]
        for g in self.grids[ii].flat:
            gi = g.id - g._id_offset
            # FLASH uses 1-indexed group info
            g.Children = [self.grids[i - 1] for i in gid[gi,7:] if i > -1]
            for g1 in g.Children:
                g1.Parent = g
            g._prepare_grid()
            g._setup_dx()
        self.max_level = self.grid_levels.max()

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

class RAMSESHierarchy(AMRHierarchy):

    grid = RAMSESGrid
    _handle = None
    
    def __init__(self,pf,data_style='ramses'):
        self.data_style = data_style
        self.field_info = RAMSESFieldContainer()
        self.parameter_file = weakref.proxy(pf)
        # for now, the hierarchy file is the parameter file!
        self.hierarchy_filename = self.parameter_file.parameter_filename
        self.directory = os.path.dirname(self.hierarchy_filename)
        self.tree_proxy = pf.ramses_tree

        self.float_type = na.float64
        AMRHierarchy.__init__(self,pf,data_style)

    def _initialize_data_storage(self):
        pass

    def _detect_fields(self):
        self.field_list = self.tree_proxy.field_names[:]
    
    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        AMRHierarchy._setup_classes(self, dd)
        self.object_types.sort()

    def _count_grids(self):
        # We have to do all the patch-coalescing here.
        level_info = self.tree_proxy.count_zones()
        num_ogrids = sum(level_info)
        ogrid_left_edge = na.zeros((num_ogrids,3), dtype='float64')
        ogrid_right_edge = na.zeros((num_ogrids,3), dtype='float64')
        ogrid_levels = na.zeros((num_ogrids,1), dtype='int32')
        ogrid_file_locations = na.zeros((num_ogrids,6), dtype='int64')
        ochild_masks = na.zeros((num_ogrids, 8), dtype='int32')
        self.tree_proxy.fill_hierarchy_arrays(
            self.pf["TopGridDimensions"],
            ogrid_left_edge, ogrid_right_edge,
            ogrid_levels, ogrid_file_locations, ochild_masks)
        # Now we can rescale
        mi, ma = ogrid_left_edge.min(), ogrid_right_edge.max()
        DL = self.pf["DomainLeftEdge"]
        DR = self.pf["DomainRightEdge"]
        ogrid_left_edge = (ogrid_left_edge - mi)/(ma - mi) * (DR - DL) + DL
        ogrid_right_edge = (ogrid_right_edge - mi)/(ma - mi) * (DR - DL) + DL
        #import pdb;pdb.set_trace()
        # We now have enough information to run the patch coalescing 
        self.proto_grids = []
        for level in xrange(len(level_info)):
            if level_info[level] == 0: continue
            ggi = (ogrid_levels == level).ravel()
            mylog.info("Re-gridding level %s: %s octree grids", level, ggi.sum())
            nd = self.pf["TopGridDimensions"] * 2**level
            dims = na.ones((ggi.sum(), 3), dtype='int64') * 2
            fl = ogrid_file_locations[ggi,:]
            # Now our initial protosubgrid
            #if level == 6: raise RuntimeError
            # We want grids that cover no more than MAX_EDGE cells in every direction
            MAX_EDGE = 128
            psgs = []
            left_index = na.rint((ogrid_left_edge[ggi,:]) * nd).astype('int64')
            right_index = left_index + 2
            lefts = [na.mgrid[0:nd[i]:MAX_EDGE] for i in range(3)]
            #lefts = zip(*[l.ravel() for l in lefts])
            pbar = get_pbar("Re-gridding ", lefts[0].size)
            min_ind = na.min(left_index, axis=0)
            max_ind = na.max(right_index, axis=0)
            for i,dli in enumerate(lefts[0]):
                pbar.update(i)
                if min_ind[0] > dli + nd[0]: continue
                if max_ind[0] < dli: continue
                idim = min(nd[0] - dli, MAX_EDGE)
                gdi = ((dli  <= right_index[:,0])
                     & (dli + idim >= left_index[:,0]))
                if not na.any(gdi): continue
                for dlj in lefts[1]:
                    if min_ind[1] > dlj + nd[1]: continue
                    if max_ind[1] < dlj: continue
                    idim = min(nd[1] - dlj, MAX_EDGE)
                    gdj = ((dlj  <= right_index[:,1])
                         & (dlj + idim >= left_index[:,1])
                         & (gdi))
                    if not na.any(gdj): continue
                    for dlk in lefts[2]:
                        if min_ind[2] > dlk + nd[2]: continue
                        if max_ind[2] < dlk: continue
                        idim = min(nd[2] - dlk, MAX_EDGE)
                        gdk = ((dlk  <= right_index[:,2])
                             & (dlk + idim >= left_index[:,2])
                             & (gdj))
                        if not na.any(gdk): continue
                        left = na.array([dli, dlj, dlk])
                        domain_left = left.ravel()
                        initial_left = na.zeros(3, dtype='int64') + domain_left
                        idims = na.ones(3, dtype='int64') * na.minimum(nd - domain_left, MAX_EDGE)
                        # We want to find how many grids are inside.
                        dleft_index = left_index[gdk,:]
                        dright_index = right_index[gdk,:]
                        ddims = dims[gdk,:]
                        dfl = fl[gdk,:]
                        psg = ramses_reader.ProtoSubgrid(initial_left, idims,
                                        dleft_index, dright_index, ddims, dfl)
                        #print "Gridding from %s to %s + %s" % (
                        #    initial_left, initial_left, idims)
                        if psg.efficiency <= 0: continue
                        self.num_deep = 0
                        psgs.extend(self._recursive_patch_splitting(
                            psg, idims, initial_left, 
                            dleft_index, dright_index, ddims, dfl))
                        #psgs.extend([psg])
            pbar.finish()
            self.proto_grids.append(psgs)
            sums = na.zeros(3, dtype='int64')
            mylog.info("Final grid count: %s", len(self.proto_grids[level]))
            if len(self.proto_grids[level]) == 1: continue
            for g in self.proto_grids[level]:
                sums += [s.sum() for s in g.sigs]
            assert(na.all(sums == dims.prod(axis=1).sum()))
        self.num_grids = sum(len(l) for l in self.proto_grids)

    num_deep = 0

    @num_deep_inc
    def _recursive_patch_splitting(self, psg, dims, ind,
            left_index, right_index, gdims, fl):
        min_eff = 0.1 # This isn't always respected.
        if self.num_deep > 40:
            # If we've recursed more than 100 times, we give up.
            psg.efficiency = min_eff
            return [psg]
        if psg.efficiency > min_eff or psg.efficiency < 0.0:
            return [psg]
        tt, ax, fp = psg.find_split()
        if (fp % 2) != 0:
            if dims[ax] != fp + 1:
                fp += 1
            else:
                fp -= 1
        #print " " * self.num_deep + "Got ax", ax, "fp", fp
        dims_l = dims.copy()
        dims_l[ax] = fp
        li_l = ind.copy()
        if na.any(dims_l <= 0): return [psg]
        L = ramses_reader.ProtoSubgrid(
                li_l, dims_l, left_index, right_index, gdims, fl)
        #print " " * self.num_deep + "L", tt, L.efficiency
        if L.efficiency > 1.0: raise RuntimeError
        if L.efficiency <= 0.0: L = []
        elif L.efficiency < min_eff:
            L = self._recursive_patch_splitting(L, dims_l, li_l,
                    left_index, right_index, gdims, fl)
        else:
            L = [L]
        dims_r = dims.copy()
        dims_r[ax] -= fp
        li_r = ind.copy()
        li_r[ax] += fp
        if na.any(dims_r <= 0): return [psg]
        R = ramses_reader.ProtoSubgrid(
                li_r, dims_r, left_index, right_index, gdims, fl)
        #print " " * self.num_deep + "R", tt, R.efficiency
        if R.efficiency > 1.0: raise RuntimeError
        if R.efficiency <= 0.0: R = []
        elif R.efficiency < min_eff:
            R = self._recursive_patch_splitting(R, dims_r, li_r,
                    left_index, right_index, gdims, fl)
        else:
            R = [R]
        return L + R
        
    def _parse_hierarchy(self):
        # We have important work to do
        grids = []
        gi = 0
        for level, grid_list in enumerate(self.proto_grids):
            for g in grid_list:
                fl = g.grid_file_locations
                props = g.get_properties()
                self.grid_left_edge[gi,:] = props[0,:] / (2.0**(level+1))
                self.grid_right_edge[gi,:] = props[1,:] / (2.0**(level+1))
                self.grid_dimensions[gi,:] = props[2,:]
                self.grid_levels[gi,:] = level
                grids.append(self.grid(gi, self, level, fl, props[0,:]))
                gi += 1
        self.grids = na.array(grids, dtype='object')

    def _get_grid_parents(self, grid, LE, RE):
        mask = na.zeros(self.num_grids, dtype='bool')
        grids, grid_ind = self.get_box_grids(LE, RE)
        mask[grid_ind] = True
        mask = na.logical_and(mask, (self.grid_levels == (grid.Level-1)).flat)
        return self.grids[mask]

    def _populate_grid_objects(self):
        for gi,g in enumerate(self.grids):
            parents = self._get_grid_parents(g,
                            self.grid_left_edge[gi,:],
                            self.grid_right_edge[gi,:])
            if len(parents) > 0:
                g.Parent.extend(parents.tolist())
                for p in parents: p.Children.append(g)
            g._prepare_grid()
            g._setup_dx()
        self.max_level = self.grid_levels.max()

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

    def _setup_data_io(self):
        self.io = io_registry[self.data_style](self.tree_proxy)

