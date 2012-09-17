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

cdef import from "server.h":
    int server()

cdef import from "client.h":
    void client()

cdef import from "meta_io.h":
    void read_particles(char *filename)
    void output_and_free_halos(np.int64_t id_offset, np.int64_t snap, 
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

def print_rockstar_settings():
    # We have to do the config
    print "FILE_FORMAT =", FILE_FORMAT
    print "PARTICLE_MASS =", PARTICLE_MASS

    print "MASS_DEFINITION =", MASS_DEFINITION
    print "MIN_HALO_OUTPUT_SIZE =", MIN_HALO_OUTPUT_SIZE
    print "FORCE_RES =", FORCE_RES

    print "SCALE_NOW =", SCALE_NOW
    print "h0 =", h0
    print "Ol =", Ol
    print "Om =", Om

    print "GADGET_ID_BYTES =", GADGET_ID_BYTES
    print "GADGET_MASS_CONVERSION =", GADGET_MASS_CONVERSION
    print "GADGET_LENGTH_CONVERSION =", GADGET_LENGTH_CONVERSION
    print "GADGET_SKIP_NON_HALO_PARTICLES =", GADGET_SKIP_NON_HALO_PARTICLES
    print "RESCALE_PARTICLE_MASS =", RESCALE_PARTICLE_MASS

    print "PARALLEL_IO =", PARALLEL_IO
    print "PARALLEL_IO_SERVER_ADDRESS =", PARALLEL_IO_SERVER_ADDRESS
    print "PARALLEL_IO_SERVER_PORT =", PARALLEL_IO_SERVER_PORT
    print "PARALLEL_IO_WRITER_PORT =", PARALLEL_IO_WRITER_PORT
    print "PARALLEL_IO_SERVER_INTERFACE =", PARALLEL_IO_SERVER_INTERFACE
    print "RUN_ON_SUCCESS =", RUN_ON_SUCCESS

    print "INBASE =", INBASE
    print "FILENAME =", FILENAME
    print "STARTING_SNAP =", STARTING_SNAP
    print "NUM_SNAPS =", NUM_SNAPS
    print "NUM_BLOCKS =", NUM_BLOCKS
    print "NUM_READERS =", NUM_READERS
    print "PRELOAD_PARTICLES =", PRELOAD_PARTICLES
    print "SNAPSHOT_NAMES =", SNAPSHOT_NAMES
    print "LIGHTCONE_ALT_SNAPS =", LIGHTCONE_ALT_SNAPS
    print "BLOCK_NAMES =", BLOCK_NAMES

    print "OUTBASE =", OUTBASE
    print "OVERLAP_LENGTH =", OVERLAP_LENGTH
    print "NUM_WRITERS =", NUM_WRITERS
    print "FORK_READERS_FROM_WRITERS =", FORK_READERS_FROM_WRITERS
    print "FORK_PROCESSORS_PER_MACHINE =", FORK_PROCESSORS_PER_MACHINE

    print "OUTPUT_FORMAT =", OUTPUT_FORMAT
    print "DELETE_BINARY_OUTPUT_AFTER_FINISHED =", DELETE_BINARY_OUTPUT_AFTER_FINISHED
    print "FULL_PARTICLE_CHUNKS =", FULL_PARTICLE_CHUNKS
    print "BGC2_SNAPNAMES =", BGC2_SNAPNAMES

    print "BOUND_PROPS =", BOUND_PROPS
    print "BOUND_OUT_TO_HALO_EDGE =", BOUND_OUT_TO_HALO_EDGE
    print "DO_MERGER_TREE_ONLY =", DO_MERGER_TREE_ONLY
    print "IGNORE_PARTICLE_IDS =", IGNORE_PARTICLE_IDS
    print "TRIM_OVERLAP =", TRIM_OVERLAP
    print "ROUND_AFTER_TRIM =", ROUND_AFTER_TRIM
    print "LIGHTCONE =", LIGHTCONE
    print "PERIODIC =", PERIODIC

    print "LIGHTCONE_ORIGIN =", LIGHTCONE_ORIGIN[0]
    print "LIGHTCONE_ORIGIN[1] =", LIGHTCONE_ORIGIN[1]
    print "LIGHTCONE_ORIGIN[2] =", LIGHTCONE_ORIGIN[2]
    print "LIGHTCONE_ALT_ORIGIN =", LIGHTCONE_ALT_ORIGIN[0]
    print "LIGHTCONE_ALT_ORIGIN[1] =", LIGHTCONE_ALT_ORIGIN[1]
    print "LIGHTCONE_ALT_ORIGIN[2] =", LIGHTCONE_ALT_ORIGIN[2]

    print "LIMIT_CENTER =", LIMIT_CENTER[0]
    print "LIMIT_CENTER[1] =", LIMIT_CENTER[1]
    print "LIMIT_CENTER[2] =", LIMIT_CENTER[2]
    print "LIMIT_RADIUS =", LIMIT_RADIUS

    print "SWAP_ENDIANNESS =", SWAP_ENDIANNESS
    print "GADGET_VARIANT =", GADGET_VARIANT

    print "FOF_FRACTION =", FOF_FRACTION
    print "FOF_LINKING_LENGTH =", FOF_LINKING_LENGTH
    print "INCLUDE_HOST_POTENTIAL_RATIO =", INCLUDE_HOST_POTENTIAL_RATIO
    print "DOUBLE_COUNT_SUBHALO_MASS_RATIO =", DOUBLE_COUNT_SUBHALO_MASS_RATIO
    print "TEMPORAL_HALO_FINDING =", TEMPORAL_HALO_FINDING
    print "MIN_HALO_PARTICLES =", MIN_HALO_PARTICLES
    print "UNBOUND_THRESHOLD =", UNBOUND_THRESHOLD
    print "ALT_NFW_METRIC =", ALT_NFW_METRIC

    print "TOTAL_PARTICLES =", TOTAL_PARTICLES
    print "BOX_SIZE =", BOX_SIZE
    print "OUTPUT_HMAD =", OUTPUT_HMAD
    print "OUTPUT_PARTICLES =", OUTPUT_PARTICLES
    print "OUTPUT_LEVELS =", OUTPUT_LEVELS
    print "DUMP_PARTICLES =", DUMP_PARTICLES[0]
    print "DUMP_PARTICLES[1] =", DUMP_PARTICLES[1]
    print "DUMP_PARTICLES[2] =", DUMP_PARTICLES[2]

    print "AVG_PARTICLE_SPACING =", AVG_PARTICLE_SPACING
    print "SINGLE_SNAP =", SINGLE_SNAP

cdef class RockstarInterface
cdef void rh_read_particles(char *filename, particle **p, np.int64_t *num_p):
    global SCALE_NOW, TOTAL_PARTICLES
    pf = rh.tsl.next()
    print 'reading from particle filename %s: %s'%(filename,pf.basename)
    cdef np.float64_t conv[6], left_edge[6]
    cdef np.ndarray[np.int64_t, ndim=1] arri
    cdef np.ndarray[np.float64_t, ndim=1] arr
    block = int(str(filename).rsplit(".")[-1])
    

    # Now we want to grab data from only a subset of the grids.
    n = rh.block_ratio
    dd = pf.h.all_data()
    SCALE_NOW = 1.0/(pf.current_redshift+1.0)
    grids = np.array_split(dd._grids, NUM_BLOCKS)[block]
    tnpart = 0
    for g in grids:
        tnpart += np.sum(dd._get_data_from_grid(g, "particle_type")==rh.dm_type)
    p[0] = <particle *> malloc(sizeof(particle) * tnpart)
    #print "Loading indices: size = ", tnpart
    conv[0] = conv[1] = conv[2] = pf["mpchcm"]
    conv[3] = conv[4] = conv[5] = 1e-5
    left_edge[0] = pf.domain_left_edge[0]
    left_edge[1] = pf.domain_left_edge[1]
    left_edge[2] = pf.domain_left_edge[2]
    left_edge[3] = left_edge[4] = left_edge[5] = 0.0
    pi = 0
    for g in grids:
        iddm = dd._get_data_from_grid(g, "particle_type")==rh.dm_type
        arri = dd._get_data_from_grid(g, "particle_index").astype("int64")
        arri = arri[iddm] #pick only DM
        npart = arri.size
        for i in range(npart):
            p[0][i+pi].id = arri[i]
        fi = 0
        for field in ["particle_position_x", "particle_position_y",
                      "particle_position_z",
                      "particle_velocity_x", "particle_velocity_y",
                      "particle_velocity_z"]:
            arr = dd._get_data_from_grid(g, field).astype("float64")
            arr = arr[iddm] #pick DM
            for i in range(npart):
                p[0][i+pi].pos[fi] = (arr[i]-left_edge[fi])*conv[fi]
            fi += 1
        pi += npart
    num_p[0] = tnpart
    TOTAL_PARTICLES = tnpart
    #print 'first particle coordinates'
    #for i in range(3):
    #    print p[0][0].pos[i],
    #print ""
    #print 'last particle coordinates'
    #for i in range(3):
    #    print p[0][tnpart-1].pos[i],
    #print ""

cdef class RockstarInterface:

    cdef public object data_source
    cdef public object ts
    cdef public object tsl
    cdef int rank
    cdef int size
    cdef public int block_ratio
    cdef public int dm_type
    cdef public int total_particles

    def __cinit__(self, ts, data_source):
        self.ts = ts
        self.tsl = ts.__iter__() #timseries generator used by read
        self.data_source = data_source

    def setup_rockstar(self, char *server_address, char *server_port,
                       int num_snaps, np.int64_t total_particles,
                       int dm_type,
                       np.float64_t particle_mass = -1.0,
                       int parallel = False, int num_readers = 1,
                       int num_writers = 1,
                       int writing_port = -1, int block_ratio = 1,
                       int periodic = 1, 
                       int min_halo_size = 25, outbase = "None"):
        global PARALLEL_IO, PARALLEL_IO_SERVER_ADDRESS, PARALLEL_IO_SERVER_PORT
        global FILENAME, FILE_FORMAT, NUM_SNAPS, STARTING_SNAP, h0, Ol, Om
        global BOX_SIZE, PERIODIC, PARTICLE_MASS, NUM_BLOCKS, NUM_READERS
        global FORK_READERS_FROM_WRITERS, PARALLEL_IO_WRITER_PORT, NUM_WRITERS
        global rh, SCALE_NOW, OUTBASE, MIN_HALO_OUTPUT_SIZE
        global OVERLAP_LENGTH
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

        if particle_mass < 0:
            particle_mass = tpf.h.grids[0]["ParticleMassMsun"][0] / h0
        PARTICLE_MASS = particle_mass
        PERIODIC = periodic
        BOX_SIZE = (tpf.domain_right_edge[0] -
                    tpf.domain_left_edge[0]) * tpf['mpchcm']
        setup_config()
        rh = self
        cdef LPG func = rh_read_particles
        set_load_particles_generic(func)

    def call_rockstar(self):
        read_particles("generic")
        rockstar(NULL, 0)
        output_and_free_halos(0, 0, 0, NULL)

    def start_server(self):
        server()

    def start_client(self):
        client()
