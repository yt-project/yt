import numpy as np
import os, sys
cimport numpy as np
cimport cython
#from cpython.mem cimport PyMem_Malloc
from libc.stdlib cimport malloc, free
import sys



# Importing relevant rockstar data types particle, fof halo, halo

cdef import from "particle.h":
    struct particle:
        np.int64_t id
        float pos[6]

cdef import from "fof.h":
    struct fof:
        np.int64_t num_p
        particle *particles

cdef import from "halo.h":
    struct halo:
        np.int64_t id
        float pos[6], corevel[3], bulkvel[3]
        float m, r, child_r, vmax_r, mgrav,    vmax, rvmax, rs, klypin_rs, vrms
        float J[3], energy, spin, alt_m[4], Xoff, Voff, b_to_a, c_to_a, A[3]
        float bullock_spin, kin_to_pot
        np.int64_t num_p, num_child_particles, p_start, desc, flags, n_core
        float min_pos_err, min_vel_err, min_bulkvel_err

# For finding sub halos import finder function and global variable
# rockstar uses to store the results

cdef import from "groupies.h":
    void find_subs(fof *f) nogil
    halo *halos
    np.int64_t num_halos
    void calc_mass_definition() nogil

# For outputing halos, rockstar style

cdef import from "meta_io.h":
    void output_halos(np.int64_t id_offset, np.int64_t snap, np.int64_t chunk, float *bounds) nogil

# For setting up the configuration of rockstar

cdef import from "config.h":
    void setup_config() nogil

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



cdef class RockstarGroupiesInterface:
    
    cdef public object pf
    cdef public object fof

    # For future use/consistency
    def __cinit__(self,pf):
        self.pf = pf

    def setup_rockstar(self,
                        particle_mass,
                        int periodic = 1, force_res=None,
                        int min_halo_size = 25, outbase = "None",
                        callbacks = None):
        global FILENAME, FILE_FORMAT, NUM_SNAPS, STARTING_SNAP, h0, Ol, Om
        global BOX_SIZE, PERIODIC, PARTICLE_MASS, NUM_BLOCKS, NUM_READERS
        global FORK_READERS_FROM_WRITERS, PARALLEL_IO_WRITER_PORT, NUM_WRITERS
        global rh, SCALE_NOW, OUTBASE, MIN_HALO_OUTPUT_SIZE
        global OVERLAP_LENGTH, TOTAL_PARTICLES, FORCE_RES
        

        if force_res is not None:
            FORCE_RES=np.float64(force_res)

        OVERLAP_LENGTH = 0.0
        
        FILENAME = "inline.<block>"
        FILE_FORMAT = "GENERIC"
        OUTPUT_FORMAT = "ASCII"
        MIN_HALO_OUTPUT_SIZE=min_halo_size
        
        pf = self.pf

        h0 = pf.hubble_constant
        Ol = pf.omega_lambda
        Om = pf.omega_matter
        
        SCALE_NOW = 1.0/(pf.current_redshift+1.0)
        
        if not outbase =='None'.decode('UTF-8'):
            #output directory. since we can't change the output filenames
            #workaround is to make a new directory
            OUTBASE = outbase 


        PARTICLE_MASS = particle_mass.in_units('Msun/h')
        PERIODIC = periodic
        BOX_SIZE = pf.domain_width.in_units('Mpccm/h')[0]

        # Set up the configuration options
        setup_config()

        # Needs to be called so rockstar can use the particle mass parameter
        # to calculate virial quantities properly
        calc_mass_definition()

    def output_halos(self):
        output_halos(0, 0, 0, NULL) 

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def make_rockstar_fof(self, np.ndarray[np.int64_t, ndim=1] pid,
                                np.ndarray[np.float64_t, ndim=2] pos,
                                np.ndarray[np.float64_t, ndim=2] vel,
                                np.ndarray[np.int64_t, ndim=1] fof_tags,
                                np.int64_t nfof,
                                np.int64_t npart_max):

        # Define fof object

        # Find number of particles
        cdef np.int64_t i, j
        cdef np.int64_t num_particles = pid.shape[0]

        # Allocate space for correct number of particles
        cdef particle* particles = <particle*> malloc(npart_max * sizeof(particle))
        cdef fof fof_obj
        fof_obj.particles = particles

        cdef np.int64_t last_fof_tag = 1
        cdef np.int64_t k = 0
        for i in range(num_particles):
            if fof_tags[i] < 0:
                continue
            if fof_tags[i] != last_fof_tag:
                last_fof_tag = fof_tags[i]
                if k > 16:
                    fof_obj.num_p = k
                    find_subs(&fof_obj)
                k = 0
            particles[k].id = pid[i]

            # fill in locations & velocities
            for j in range(3):
                particles[k].pos[j] = pos[i,j]
                particles[k].pos[j+3] = vel[i,j]
            k += 1
        free(particles)



