"""
Wrapping code for Oliver Hahn's RamsesRead++

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD
Author: Oliver Hahn <ohahn@stanford.edu>
Affiliation: KIPAC / Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2010 Matthew Turk.  All Rights Reserved.

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

# Cython wrapping code for Oliver Hahn's RAMSES reader
from cython.operator cimport dereference as deref, preincrement as inc
from libc.stdlib cimport malloc, free, abs

import numpy as np
cimport numpy as np
cimport cython

cdef extern from "<vector>" namespace "std":
    cdef cppclass vector[T]:
        cppclass iterator:
            T operator*()
            iterator operator++()
            bint operator==(iterator)
            bint operator!=(iterator)
        vector()
        void push_back(T&)
        T& operator[](int)
        T& at(int)
        iterator begin()
        iterator end()

cdef extern from "<map>" namespace "std":
    cdef cppclass map[A,B]:
        pass

cdef extern from "string" namespace "std":
    cdef cppclass string:
        string(char *cstr)
        char *c_str()
        string operator*()

cdef extern from "RAMSES_typedefs.h":
    pass

cdef extern from "RAMSES_info.hh" namespace "RAMSES":
    enum codeversion:
        version1
        version2
        version3

    cdef cppclass snapshot:
        string m_filename
        codeversion m_version
        cppclass info_data:
            unsigned ncpu
            unsigned ndim
            unsigned levelmin
            unsigned levelmax
            unsigned ngridmax
            unsigned nstep_coarse
            double boxlen
            double time
            double aexp
            double H0
            double omega_m
            double omega_l
            double omega_k
            double omega_b
            double unit_l
            double unit_d
            double unit_t

        info_data m_header
        vector[double] ind_min
        vector[double] ind_max

        snapshot(string info_filename, codeversion ver)

        unsigned get_snapshot_num()
        unsigned getdomain_bykey( double key )

    #void hilbert3d(vector[double] &x, vector[double] &y, vector[double] &z,
    #               vector[double] &order, unsigned bit_length)

cdef extern from "RAMSES_amr_data.hh" namespace "RAMSES::AMR":
    cdef cppclass vec[real_t]:
        real_t x, y, z
        vec( real_t x_, real_t y_, real_t z_)
        vec ( vec& v)
        vec ( )

    # This class definition is out of date.  I have unrolled the template
    # below.
    cdef cppclass cell_locally_essential[id_t, real_t]:
        id_t m_neighbor[6]
        id_t m_father
        id_t m_son[8]
        real_t m_xg[3]
        id_t m_cpu

        char m_pos

        cell_locally_essential()

        bint is_refined(int ison)

    cdef cppclass RAMSES_cell:
        unsigned m_neighbour[6]
        unsigned m_father
        unsigned m_son[8]
        float m_xg[3]
        unsigned m_cpu

        char m_pos

        RAMSES_cell()

        bint is_refined(int ison)

    cdef cppclass cell_simple[id_t, real_t]:
        id_t m_son[1]
        real_t m_xg[3]
        id_t m_cpu

        char m_pos
        cell_simple()
        bint is_refined(int ison)

    # AMR level implementation

    # This class definition is out of date.  I have unrolled the template
    # below.
    cdef cppclass level[Cell_]:
        unsigned m_ilevel
        vector[Cell_] m_level_cells
        double m_xc[8]
        double m_yc[8]
        double m_zc[8]

        # I skipped the typedefs here

        double m_dx
        unsigned m_nx
        level (unsigned ilevel)

        void register_cell( Cell_ cell )
        vector[Cell_].iterator begin()
        vector[Cell_].iterator end()

        Cell_& operator[]( unsigned i)
        unsigned size()

    cdef cppclass RAMSES_level:
        unsigned m_ilevel
        vector[RAMSES_cell] m_level_cells
        double m_xc[8]
        double m_yc[8]
        double m_zc[8]

        # I skipped the typedefs here

        double m_dx
        unsigned m_nx
        RAMSES_level (unsigned ilevel)

        void register_cell( RAMSES_cell cell )
        vector[RAMSES_cell].iterator begin()
        vector[RAMSES_cell].iterator end()

        RAMSES_cell& operator[]( unsigned i)
        unsigned size()

    # This class definition is out of date.  I have unrolled the template
    # below.
    cdef cppclass tree[Cell_, Level_]:
        cppclass header:
            vector[int] nx
            vector[int] nout
            vector[int] nsteps
            int ncpu
            int ndim
            int nlevelmax
            int ngridmax
            int nboundary
            int ngrid_current
            double boxlen
            vector[double] tout
            vector[double] aout
            vector[double] dtold
            vector[double] dtnew
            vector[double] cosm
            vector[double] timing
            double t
            double mass_sph

        vector[Level_] m_AMR_levels
        vector[unsigned] mheadl, m_numbl, m_taill

        int m_cpu
        int m_minlevel
        int m_maxlevel
        string m_fname
        unsigned m_ncoarse
        header m_header

        # This is from later on in the .hh file ... I don't think we need them
        # typedef proto_iterator<const tree*> const_iterator;
        # typedef proto_iterator<tree *> iterator;

        tree (snapshot &snap, int cpu, int maxlevel, int minlevel = 1)

    cppclass tree_iterator "RAMSES::AMR::RAMSES_tree::iterator":
        tree_iterator operator*()
        RAMSES_cell operator*()
        tree_iterator begin()
        tree_iterator end()
        tree_iterator to_parent()
        tree_iterator get_parent()
        void next()
        bint operator!=(tree_iterator other)
        unsigned get_cell_father()
        bint is_refined(int ison)
        int get_absolute_position()

    cdef cppclass RAMSES_tree:
        # This is, strictly speaking, not a header.  But, I believe it is
        # going to work alright.
        cppclass header:
            vector[int] nx
            vector[int] nout
            vector[int] nsteps
            int ncpu
            int ndim
            int nlevelmax
            int ngridmax
            int nboundary
            int ngrid_current
            double boxlen
            vector[double] tout
            vector[double] aout
            vector[double] dtold
            vector[double] dtnew
            vector[double] cosm
            vector[double] timing
            double t
            double mass_sph

        vector[RAMSES_level] m_AMR_levels
        vector[unsigned] mheadl, m_numbl, m_taill

        int m_cpu
        int m_minlevel
        int m_maxlevel
        string m_fname
        unsigned m_ncoarse
        header m_header

        unsigned size()

        # This is from later on in the .hh file ... I don't think we need them
        # typedef proto_iterator<const tree*> const_iterator;
        # typedef proto_iterator<tree *> iterator;

        RAMSES_tree(snapshot &snap, int cpu, int maxlevel, int minlevel)
        void read()

        tree_iterator begin(int ilevel)
        tree_iterator end(int ilevel)

        tree_iterator begin()
        tree_iterator end()

        vec[double] cell_pos_double "cell_pos<double>" (tree_iterator it, unsigned ind) 
        vec[double] grid_pos_double "grid_pos<double>" (tree_iterator it)
        vec[float] cell_pos_float "cell_pos<float>" (tree_iterator it, unsigned ind) 
        vec[float] grid_pos_float "grid_pos<float>" (tree_iterator it)

cdef extern from "RAMSES_amr_data.hh" namespace "RAMSES::HYDRO":
    enum hydro_var:
        density
        velocity_x
        velocity_y
        velocity_z
        pressure
        metallicit

    char ramses_hydro_variables[][64]

    # There are a number of classes here that we could wrap and utilize.
    # However, I will only wrap the methods I know we need.

    # I have no idea if this will work.
    cdef cppclass TreeTypeIterator[TreeType_]:
        pass

    # This class definition is out of date.  I have unrolled the template
    # below.
    cdef cppclass data[TreeType_, Real_]:
        cppclass header:
            unsigned ncpu
            unsigned nvar
            unsigned ndim
            unsigned nlevelmax
            unsigned nboundary
            double gamma
        string m_fname
        header m_header

        # I don't want to implement proto_data, so we put this here
        Real_& cell_value(TreeTypeIterator[TreeType_] &it, int ind)

        unsigned m_nvars
        vector[string] m_varnames
        map[string, unsigned] m_var_name_map

        data(TreeType_ &AMRtree)
        #_OutputIterator get_var_names[_OutputIterator](_OutputIterator names)
        void read(string varname)

    cdef cppclass RAMSES_hydro_data:
        cppclass header:
            unsigned ncpu
            unsigned nvar
            unsigned ndim
            unsigned nlevelmax
            unsigned nboundary
            double gamma
        string m_fname
        header m_header

        # I don't want to implement proto_data, so we put this here
        double cell_value (tree_iterator &it, int ind)

        unsigned m_nvars
        vector[string] m_varnames
        map[string, unsigned] m_var_name_map

        RAMSES_hydro_data(RAMSES_tree &AMRtree)
        #_OutputIterator get_var_names[_OutputIterator](_OutputIterator names)
        void read(string varname)

cdef class RAMSES_tree_proxy:
    cdef string *snapshot_name
    cdef snapshot *rsnap

    cdef RAMSES_tree** trees
    cdef RAMSES_hydro_data*** hydro_datas

    cdef int *loaded

    cdef public object field_ind
    cdef public object field_names

    # We will store this here so that we have a record, independent of the
    # header, of how many things we have allocated
    cdef int ndomains, nfields
    
    def __cinit__(self, char *fn):
        cdef int idomain, ifield, ii
        cdef RAMSES_tree *local_tree
        cdef RAMSES_hydro_data *local_hydro_data
        self.snapshot_name = new string(fn)
        self.rsnap = new snapshot(deref(self.snapshot_name), version3)
        # We now have to get our field names to fill our array
        self.trees = <RAMSES_tree**>\
            malloc(sizeof(RAMSES_tree*) * self.rsnap.m_header.ncpu)
        self.hydro_datas = <RAMSES_hydro_data ***>\
                       malloc(sizeof(RAMSES_hydro_data**) * self.rsnap.m_header.ncpu)
        self.ndomains = self.rsnap.m_header.ncpu
        #for ii in range(self.ndomains): self.trees[ii] = NULL
        for idomain in range(1, self.rsnap.m_header.ncpu + 1):
            local_tree = new RAMSES_tree(deref(self.rsnap), idomain,
                                         self.rsnap.m_header.levelmax, 0)
            local_tree.read()
            local_hydro_data = new RAMSES_hydro_data(deref(local_tree))
            self.hydro_datas[idomain - 1] = <RAMSES_hydro_data **>\
                malloc(sizeof(RAMSES_hydro_data*) * local_hydro_data.m_nvars)
            del local_hydro_data
            for ii in range(local_hydro_data.m_nvars):
                self.hydro_datas[idomain - 1][ii] = \
                    new RAMSES_hydro_data(deref(local_tree))
            self.trees[idomain - 1] = local_tree
            # We do not delete anything
        # Only once, we read all the field names
        self.nfields = local_hydro_data.m_nvars
        cdef string *field_name
        self.field_names = []
        self.field_ind = {}
        self.loaded = <int *> malloc(sizeof(int) * local_hydro_data.m_nvars)
        for ifield in range(local_hydro_data.m_nvars):
            field_name = &(local_hydro_data.m_varnames[ifield])
            # Does this leak?
            self.field_names.append(field_name.c_str())
            self.field_ind[self.field_names[-1]] = ifield
            self.loaded[ifield] = 0
        # This all needs to be cleaned up in the deallocator

    def __dealloc__(self):
        cdef int idomain, ifield
        for idomain in range(self.ndomains):
            for ifield in range(self.nfields):
                if self.hydro_datas[idomain][ifield] != NULL:
                    del self.hydro_datas[idomain][ifield]
            if self.trees[idomain] != NULL:
                del self.trees[idomain]
            free(self.hydro_datas[idomain])
        free(self.trees)
        free(self.hydro_datas)
        free(self.loaded)
        if self.snapshot_name != NULL: del self.snapshot_name
        if self.rsnap != NULL: del self.rsnap
        
    def count_zones(self):
        # We need to do simulation domains here

        cdef unsigned idomain, ilevel
        cdef RAMSES_tree *local_tree
        cdef RAMSES_hydro_data *local_hydro_data
        cdef RAMSES_level *local_level

        # All the loop-local pointers must be declared up here

        cell_count = []
        for ilevel in range(self.rsnap.m_header.levelmax + 1):
            cell_count.append(0)
        for idomain in range(1, self.rsnap.m_header.ncpu + 1):
            local_tree = new RAMSES_tree(deref(self.rsnap), idomain,
                                         self.rsnap.m_header.levelmax, 0)
            local_tree.read()
            local_hydro_data = new RAMSES_hydro_data(deref(local_tree))
            for ilevel in range(local_tree.m_maxlevel + 1):
                local_level = &local_tree.m_AMR_levels[ilevel]
                cell_count[ilevel] += local_level.size()
            del local_tree, local_hydro_data

        return cell_count

    def ensure_loaded(self, char *varname, int domain_index):
        # this domain_index must be zero-indexed
        cdef int varindex = self.field_ind[varname]
        cdef string *field_name = new string(varname)
        if self.loaded[varindex] == 1: return
        print "READING FROM DISK"
        self.hydro_datas[domain_index][varindex].read(deref(field_name))
        self.loaded[varindex] = 1
        del field_name

    def clear_tree(self, char *varname, int domain_index):
        # this domain_index must be zero-indexed
        # We delete and re-create
        cdef int varindex = self.field_ind[varname]
        cdef string *field_name = new string(varname)
        if self.loaded[varindex] == 0: return
        del self.hydro_datas[domain_index][varindex]
        self.hydro_datas[domain_index - 1][varindex] = \
            new RAMSES_hydro_data(deref(self.trees[domain_index]))
        self.loaded[varindex] = 0
        del field_name

    def get_file_info(self):
        header_info = {}
        header_info["ncpu"] = self.rsnap.m_header.ncpu
        header_info["ndim"] = self.rsnap.m_header.ndim
        header_info["levelmin"] = self.rsnap.m_header.levelmin
        header_info["levelmax"] = self.rsnap.m_header.levelmax
        header_info["ngridmax"] = self.rsnap.m_header.ngridmax
        header_info["nstep_coarse"] = self.rsnap.m_header.nstep_coarse
        header_info["boxlen"] = self.rsnap.m_header.boxlen
        header_info["time"] = self.rsnap.m_header.time
        header_info["aexp"] = self.rsnap.m_header.aexp
        header_info["H0"] = self.rsnap.m_header.H0
        header_info["omega_m"] = self.rsnap.m_header.omega_m
        header_info["omega_l"] = self.rsnap.m_header.omega_l
        header_info["omega_k"] = self.rsnap.m_header.omega_k
        header_info["omega_b"] = self.rsnap.m_header.omega_b
        header_info["unit_l"] = self.rsnap.m_header.unit_l
        header_info["unit_d"] = self.rsnap.m_header.unit_d
        header_info["unit_t"] = self.rsnap.m_header.unit_t
        return header_info

    def fill_hierarchy_arrays(self, 
                              np.ndarray[np.float64_t, ndim=2] left_edges,
                              np.ndarray[np.float64_t, ndim=2] right_edges,
                              np.ndarray[np.int32_t, ndim=2] grid_levels,
                              np.ndarray[np.int64_t, ndim=2] grid_file_locations,
                              np.ndarray[np.int32_t, ndim=2] child_mask):
        # We need to do simulation domains here

        cdef unsigned idomain, ilevel

        # All the loop-local pointers must be declared up here
        cdef RAMSES_tree *local_tree
        cdef RAMSES_hydro_data *local_hydro_data
        cdef unsigned father

        cdef tree_iterator grid_it, grid_end, father_it
        cdef vec[double] gvec
        cdef int grid_ind = 0
        cdef unsigned parent_ind
        cdef bint ci

        cdef double grid_half_width 

        cdef np.int32_t rr
        cell_count = []
        level_cell_counts = {}
        for idomain in range(1, self.rsnap.m_header.ncpu + 1):
            local_tree = new RAMSES_tree(deref(self.rsnap), idomain,
                                         self.rsnap.m_header.levelmax, 0)
            local_tree.read()
            local_hydro_data = new RAMSES_hydro_data(deref(local_tree))
            for ilevel in range(local_tree.m_maxlevel + 1):
                # this gets overwritten for every domain, which is okay
                level_cell_counts[ilevel] = grid_ind 
                grid_half_width = self.rsnap.m_header.boxlen / (2**(ilevel + 1))
                grid_it = local_tree.begin(ilevel)
                grid_end = local_tree.end(ilevel)
                while grid_it != grid_end:
                    gvec = local_tree.grid_pos_double(grid_it)
                    left_edges[grid_ind, 0] = gvec.x - grid_half_width
                    left_edges[grid_ind, 1] = gvec.y - grid_half_width
                    left_edges[grid_ind, 2] = gvec.z - grid_half_width
                    right_edges[grid_ind, 0] = gvec.x + grid_half_width
                    right_edges[grid_ind, 1] = gvec.y + grid_half_width
                    right_edges[grid_ind, 2] = gvec.z + grid_half_width
                    grid_levels[grid_ind, 0] = ilevel
                    # Now the harder part
                    father_it = grid_it.get_parent()
                    grid_file_locations[grid_ind, 0] = idomain
                    grid_file_locations[grid_ind, 1] = grid_ind - level_cell_counts[ilevel]
                    parent_ind = father_it.get_absolute_position()
                    if ilevel > 0:
                        # We calculate the REAL parent index
                        grid_file_locations[grid_ind, 2] = \
                            level_cell_counts[ilevel - 1] + parent_ind
                    else:
                        grid_file_locations[grid_ind, 2] = -1
                    for ci in range(8):
                        rr = <np.int32_t> grid_it.is_refined(ci)
                        child_mask[grid_ind, ci] = rr
                    grid_ind += 1
                    grid_it.next()
            del local_tree, local_hydro_data

    def read_grid(self, char *field, int level, int domain, int grid_id):

        self.ensure_loaded(field, domain - 1)
        cdef int varindex = self.field_ind[field]
        cdef int i

        cdef np.ndarray[np.float64_t, ndim=3] tr = np.empty((2,2,2), dtype='float64')
        cdef tree_iterator grid_it, grid_end
        cdef double* data = <double*> tr.data

        cdef RAMSES_tree *local_tree = self.trees[domain - 1]
        cdef RAMSES_hydro_data *local_hydro_data = self.hydro_datas[domain - 1][varindex]
        grid_it = local_tree.begin(level)
        grid_end = local_tree.end(level)
        
        cdef int cur_id = 0

        while grid_it != grid_end and cur_id < grid_id:
            cur_id += 1
            grid_it.next()
        for i in range(8): 
            data[i] = local_hydro_data.cell_value(grid_it, i)
        return tr
