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

        bool is_refined(int ison)

    cdef cppclass RAMSES_cell:
        unsigned m_neighbor[6]
        unsigned m_father
        unsigned m_son[8]
        float m_xg[3]
        unsigned m_cpu

        char m_pos

        RAMSES_cell()

        bool is_refined(int ison)

    cdef cppclass cell_simple[id_t, real_t]:
        id_t m_son[1]
        real_t m_xg[3]
        id_t m_cpu

        char m_pos
        cell_simple()
        bool is_refined(int ison)

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
        tree_iterator begin()
        tree_iterator end()
        void next()
        bint operator!=(tree_iterator other)

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
        double& cell_value(TreeTypeIterator[RAMSES_tree] &it, int ind)

        unsigned m_nvars
        vector[string] m_varnames
        map[string, unsigned] m_var_name_map

        RAMSES_hydro_data(RAMSES_tree &AMRtree)
        #_OutputIterator get_var_names[_OutputIterator](_OutputIterator names)
        void read(string varname)

def get_file_info(char *fn):
    cdef string *sfn = new string(fn)
    cdef snapshot *rsnap = new snapshot(deref(sfn), version3)
    header_info = {}
    header_info["ncpu"] = rsnap.m_header.ncpu
    header_info["ndim"] = rsnap.m_header.ndim
    header_info["levelmin"] = rsnap.m_header.levelmin
    header_info["levelmax"] = rsnap.m_header.levelmax
    header_info["ngridmax"] = rsnap.m_header.ngridmax
    header_info["nstep_coarse"] = rsnap.m_header.nstep_coarse
    header_info["boxlen"] = rsnap.m_header.boxlen
    header_info["time"] = rsnap.m_header.time
    header_info["aexp"] = rsnap.m_header.aexp
    header_info["H0"] = rsnap.m_header.H0
    header_info["omega_m"] = rsnap.m_header.omega_m
    header_info["omega_l"] = rsnap.m_header.omega_l
    header_info["omega_k"] = rsnap.m_header.omega_k
    header_info["omega_b"] = rsnap.m_header.omega_b
    header_info["unit_l"] = rsnap.m_header.unit_l
    header_info["unit_d"] = rsnap.m_header.unit_d
    header_info["unit_t"] = rsnap.m_header.unit_t

    del rsnap
    del sfn
    return header_info

def fill_buffer(char *fn, int image_levelmin, int image_levelmax,
                unsigned minlvl, unsigned maxlvl,
                double xmin, double xmax, double ymin, double ymax,
                double zmin, double zmax):
    cdef string *sfn = new string(fn)
    cdef snapshot *rsnap = new snapshot(deref(sfn), version3)
    # We need to do simulation domains here
    cdef int nx = <int>( 2**image_levelmax + 0.5 )
    cdef int ny = nx

    cdef np.ndarray[np.int64_t, ndim=2] count
    cdef np.ndarray[np.float64_t, ndim=2] image
    count = np.zeros( (nx, ny), dtype='int64')
    image = np.zeros( (nx, ny), dtype='float64')

    cdef unsigned idomain, ilevel
    cdef int cellsize, i

    # All the loop-local pointers must be declared up here
    cdef RAMSES_tree *local_tree
    cdef RAMSES_hydro_data *local_hydro_data
    cdef string *field_name = new string("density")

    cdef tree_iterator grid_it, grid_end

    print "idomain in range ", rsnap.m_header.ncpu
    for idomain in range(1, rsnap.m_header.ncpu + 1):
        print "Handling domain", idomain
        local_tree = new RAMSES_tree(deref(rsnap), idomain, maxlvl, minlvl)
        local_tree.read()
        local_hydro_data = new RAMSES_hydro_data(deref(local_tree))
        print "Reading"
        local_hydro_data.read(deref(field_name))
        print "Finished reading"
        for ilevel in range(image_levelmin, image_levelmax + 1):
            cellsize = <int>(nx/2**(ilevel+1))
            grid_it = local_tree.begin(ilevel)
            grid_end = local_tree.end(ilevel)
            i = 0
            while grid_it != grid_end:
                i += 1
                grid_it.next()
            print "on level %s we have %s grids" % (ilevel, i)
        del local_tree, local_hydro_data

    del rsnap
    del sfn

def count_zones(char *fn, unsigned minlvl, unsigned maxlvl):
    cdef string *sfn = new string(fn)
    cdef snapshot *rsnap = new snapshot(deref(sfn), version3)
    # We need to do simulation domains here

    cdef unsigned idomain, ilevel
    cdef int i

    # All the loop-local pointers must be declared up here
    cdef RAMSES_tree *local_tree
    cdef RAMSES_hydro_data *local_hydro_data
    cdef string *field_name = new string("density")

    cdef tree_iterator grid_it, grid_end

    print "idomain in range ", rsnap.m_header.ncpu
    for idomain in range(1, rsnap.m_header.ncpu + 1):
        print "Handling domain", idomain
        local_tree = new RAMSES_tree(deref(rsnap), idomain, maxlvl, minlvl)
        local_tree.read()
        local_hydro_data = new RAMSES_hydro_data(deref(local_tree))
        print "Reading"
        local_hydro_data.read(deref(field_name))
        print "Finished reading"
        for ilevel in range(minlvl, maxlvl + 1):
            grid_it = local_tree.begin(ilevel)
            grid_end = local_tree.end(ilevel)
            i = 0
            while grid_it != grid_end:
                i += 1
                grid_it.next()
            print "on level %s we have %s grids" % (ilevel, i)
        del local_tree, local_hydro_data

    del rsnap
    del sfn
