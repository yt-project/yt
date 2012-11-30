"""
Particle operations for Lagrangian Volume

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2011 Matthew Turk.  All Rights Reserved.

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

import numpy as np
import os, sys
cimport numpy as np
cimport cython
from libc.stdlib cimport malloc

from yt.config import ytcfg

cdef import from "particle.h":
    struct particle:
        np.int64_t id
        float pos[6]

cdef import from "io_generic.h":
    ctypedef void (*LPG) (char *filename, particle **p, np.int64_t *num_p)
    void set_load_particles_generic(LPG func)

cdef import from "rockstar.h":
    void rockstar(float *bounds, np.int64_t manual_subs)

cdef import from "config.h":
    void setup_config()

cdef import from "server.h" nogil:
    int server()
    np.int64_t READER_TYPE
    np.int64_t WRITER_TYPE

cdef import from "client.h" nogil:
    void client(np.int64_t in_type)

cdef import from "meta_io.h":
    void read_particles(char *filename)
    void output_halos(np.int64_t id_offset, np.int64_t snap, 
			   np.int64_t chunk, float *bounds)

cdef import from "config_vars.h":
    # Rockstar cleverly puts all of the config variables inside a templated
    # definition of their vaiables.
    char *FILE_FORMAT
    np.float64_t PARTICLE_MASS

    char *MASS_DEFINITION
    np.int64_t MIN_HALO_OUTPUT_SIZE
    np.float64_t FORCE_RES

    np.float64_t SCALE_NOW
    np.float64_t h0
    np.float64_t Ol
    np.float64_t Om

    np.int64_t GADGET_ID_BYTES
    np.float64_t GADGET_MASS_CONVERSION
    np.float64_t GADGET_LENGTH_CONVERSION
    np.int64_t GADGET_SKIP_NON_HALO_PARTICLES
    np.int64_t RESCALE_PARTICLE_MASS

    np.int64_t PARALLEL_IO
    char *PARALLEL_IO_SERVER_ADDRESS
    char *PARALLEL_IO_SERVER_PORT
    np.int64_t PARALLEL_IO_WRITER_PORT
    char *PARALLEL_IO_SERVER_INTERFACE
    char *RUN_ON_SUCCESS

    char *INBASE
    char *FILENAME
    np.int64_t STARTING_SNAP
    np.int64_t NUM_SNAPS
    np.int64_t NUM_BLOCKS
    np.int64_t NUM_READERS
    np.int64_t PRELOAD_PARTICLES
    char *SNAPSHOT_NAMES
    char *LIGHTCONE_ALT_SNAPS
    char *BLOCK_NAMES

    char *OUTBASE
    np.float64_t OVERLAP_LENGTH
    np.int64_t NUM_WRITERS
    np.int64_t FORK_READERS_FROM_WRITERS
    np.int64_t FORK_PROCESSORS_PER_MACHINE

    char *OUTPUT_FORMAT
    np.int64_t DELETE_BINARY_OUTPUT_AFTER_FINISHED
    np.int64_t FULL_PARTICLE_CHUNKS
    char *BGC2_SNAPNAMES

    np.int64_t BOUND_PROPS
    np.int64_t BOUND_OUT_TO_HALO_EDGE
    np.int64_t DO_MERGER_TREE_ONLY
    np.int64_t IGNORE_PARTICLE_IDS
    np.float64_t TRIM_OVERLAP
    np.float64_t ROUND_AFTER_TRIM
    np.int64_t LIGHTCONE
    np.int64_t PERIODIC

    np.float64_t LIGHTCONE_ORIGIN[3]
    np.float64_t LIGHTCONE_ALT_ORIGIN[3]

    np.float64_t LIMIT_CENTER[3]
    np.float64_t LIMIT_RADIUS

    np.int64_t SWAP_ENDIANNESS
    np.int64_t GADGET_VARIANT

    np.float64_t FOF_FRACTION
    np.float64_t FOF_LINKING_LENGTH
    np.float64_t INCLUDE_HOST_POTENTIAL_RATIO
    np.float64_t DOUBLE_COUNT_SUBHALO_MASS_RATIO
    np.int64_t TEMPORAL_HALO_FINDING
    np.int64_t MIN_HALO_PARTICLES
    np.float64_t UNBOUND_THRESHOLD
    np.int64_t ALT_NFW_METRIC

    np.int64_t TOTAL_PARTICLES
    np.float64_t BOX_SIZE
    np.int64_t OUTPUT_HMAD
    np.int64_t OUTPUT_PARTICLES
    np.int64_t OUTPUT_LEVELS
    np.float64_t DUMP_PARTICLES[3]

    np.float64_t AVG_PARTICLE_SPACING
    np.int64_t SINGLE_SNAP

# Forward declare
cdef class RockstarInterface

cdef void rh_read_particles(char *filename, particle **p, np.int64_t *num_p) with gil:
    global SCALE_NOW
    cdef np.float64_t conv[6], left_edge[6]
    cdef np.ndarray[np.int64_t, ndim=1] arri
    cdef np.ndarray[np.float64_t, ndim=1] arr
    cdef unsigned long long pi,fi,i
    pf = rh.tsl.next()
    print 'reading from particle filename %s: %s'%(filename,pf.basename)
    block = int(str(filename).rsplit(".")[-1])
    n = rh.block_ratio

    SCALE_NOW = 1.0/(pf.current_redshift+1.0)
    # Now we want to grab data from only a subset of the grids for each reader.
    all_fields = set(pf.h.derived_field_list + pf.h.field_list)

    # First we need to find out how many this reader is going to read in
    # if the number of readers > 1.
    if NUM_BLOCKS > 1:
        local_parts = 0
        for g in pf.h._get_objs("grids"):
            if g.NumberOfParticles == 0: continue
            if "particle_type" in all_fields:
                iddm = g["particle_type"] == rh.dm_type
            else:
                iddm = Ellipsis
            arri = g["particle_index"].astype("int64")
            arri = arri[iddm] #pick only DM
            local_parts += arri.size
    else:
        local_parts = TOTAL_PARTICLES

    #print "local_parts", local_parts

    p[0] = <particle *> malloc(sizeof(particle) * local_parts)

    conv[0] = conv[1] = conv[2] = pf["mpchcm"]
    conv[3] = conv[4] = conv[5] = 1e-5
    left_edge[0] = pf.domain_left_edge[0]
    left_edge[1] = pf.domain_left_edge[1]
    left_edge[2] = pf.domain_left_edge[2]
    left_edge[3] = left_edge[4] = left_edge[5] = 0.0
    pi = 0
    for g in pf.h._get_objs("grids"):
        if g.NumberOfParticles == 0: continue
        if "particle_type" in all_fields:
            iddm = g["particle_type"] == rh.dm_type
        else:
            iddm = Ellipsis
        arri = g["particle_index"].astype("int64")
        arri = arri[iddm] #pick only DM
        npart = arri.size
        for i in range(npart):
            p[0][i+pi].id = arri[i]
        fi = 0
        for field in ["particle_position_x", "particle_position_y",
                      "particle_position_z",
                      "particle_velocity_x", "particle_velocity_y",
                      "particle_velocity_z"]:
            arr = g[field].astype("float64")
            arr = arr[iddm] #pick DM
            for i in range(npart):
                p[0][i+pi].pos[fi] = (arr[i]-left_edge[fi])*conv[fi]
            fi += 1
        pi += npart
    num_p[0] = local_parts

cdef class RockstarInterface:

    cdef public object data_source
    cdef public object ts
    cdef public object tsl
    cdef int rank
    cdef int size
    cdef public int block_ratio
    cdef public int dm_type
    cdef public int total_particles

    def __cinit__(self, ts):
        self.ts = ts
        self.tsl = ts.__iter__() #timseries generator used by read

    def setup_rockstar(self, char *server_address, char *server_port,
                       int num_snaps, np.int64_t total_particles,
                       int dm_type,
                       np.float64_t particle_mass,
                       int parallel = False, int num_readers = 1,
                       int num_writers = 1,
                       int writing_port = -1, int block_ratio = 1,
                       int periodic = 1, force_res=None,
                       int min_halo_size = 25, outbase = "None"):
        global PARALLEL_IO, PARALLEL_IO_SERVER_ADDRESS, PARALLEL_IO_SERVER_PORT
        global FILENAME, FILE_FORMAT, NUM_SNAPS, STARTING_SNAP, h0, Ol, Om
        global BOX_SIZE, PERIODIC, PARTICLE_MASS, NUM_BLOCKS, NUM_READERS
        global FORK_READERS_FROM_WRITERS, PARALLEL_IO_WRITER_PORT, NUM_WRITERS
        global rh, SCALE_NOW, OUTBASE, MIN_HALO_OUTPUT_SIZE
        global OVERLAP_LENGTH, TOTAL_PARTICLES, FORCE_RES
        if force_res is not None:
            FORCE_RES=np.float64(force_res)
            #print "set force res to ",FORCE_RES
        OVERLAP_LENGTH = 0.0
        if parallel:
            PARALLEL_IO = 1
            PARALLEL_IO_SERVER_ADDRESS = server_address
            PARALLEL_IO_SERVER_PORT = server_port
            if writing_port > 0:
                PARALLEL_IO_WRITER_PORT = writing_port
        else:
            PARALLEL_IO = 0
            PARALLEL_IO_SERVER_ADDRESS = server_address
            PARALLEL_IO_SERVER_PORT = server_port
        FILENAME = "inline.<block>"
        FILE_FORMAT = "GENERIC"
        OUTPUT_FORMAT = "ASCII"
        NUM_SNAPS = num_snaps
        NUM_READERS = num_readers
        NUM_WRITERS = num_writers
        NUM_BLOCKS = num_readers
        MIN_HALO_OUTPUT_SIZE=min_halo_size
        TOTAL_PARTICLES = total_particles
        self.block_ratio = block_ratio
        
        tpf = self.ts[0]
        h0 = tpf.hubble_constant
        Ol = tpf.omega_lambda
        Om = tpf.omega_matter
        SCALE_NOW = 1.0/(tpf.current_redshift+1.0)
        if not outbase =='None'.decode('UTF-8'):
            #output directory. since we can't change the output filenames
            #workaround is to make a new directory
            OUTBASE = outbase 

        PARTICLE_MASS = particle_mass
        PERIODIC = periodic
        BOX_SIZE = (tpf.domain_right_edge[0] -
                    tpf.domain_left_edge[0]) * tpf['mpchcm']
        setup_config()
        rh = self
        rh.dm_type = dm_type
        cdef LPG func = rh_read_particles
        set_load_particles_generic(func)

    def call_rockstar(self):
        read_particles("generic")
        rockstar(NULL, 0)
        output_halos(0, 0, 0, NULL)

    def start_server(self):
        with nogil:
            server()

    def start_reader(self):
        cdef np.int64_t in_type = np.int64(READER_TYPE)
        client(in_type)

    def start_writer(self):
        cdef np.int64_t in_type = np.int64(WRITER_TYPE)
        client(in_type)
