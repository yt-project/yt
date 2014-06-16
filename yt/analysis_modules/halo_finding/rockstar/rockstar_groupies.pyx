import numpy as np
import os, sys
cimport numpy as np
cimport cython
#from cpython.mem cimport PyMem_Malloc
from libc.stdlib cimport malloc, free
import sys

ctypedef fused anyfloat:
    np.float32_t
    np.float64_t

# Importing relevant rockstar data types particle, fof halo, halo

cdef import from "particle.h":
    struct particle:
        np.int64_t id
        float pos[6]

cdef import from "rockstar.h":
    particle *global_particles "p"
    void rockstar_cleanup()

cdef import from "fof.h":
    struct fof:
        np.int64_t num_p
        particle *particles

cdef import from "halo.h":
    struct halo:
        np.int64_t id
        float pos[6]
        float corevel[3]
        float bulkvel[3]
        float m, r, child_r, vmax_r, mgrav,    vmax, rvmax, rs, klypin_rs, vrms
        float J[3]
        float energy, spin
        float alt_m[4]
        float Xoff, Voff, b_to_a, c_to_a
        float A[3]
        float bullock_spin, kin_to_pot
        np.int64_t num_p, num_child_particles, p_start, desc, flags, n_core
        float min_pos_err, min_vel_err, min_bulkvel_err

ctypedef packed struct haloflat:
    np.int64_t id
    float pos_x, pos_y, pos_z, pos_v, pos_u, pos_w
    float corevel_x, corevel_y, corevel_z
    float bulkvel_x, bulkvel_y, bulkvel_z
    float m, r, child_r, vmax_r, mgrav,    vmax, rvmax, rs, klypin_rs, vrms
    float J1, J2, J3
    float energy, spin
    float alt_m1, alt_m2, alt_m3, alt_m4
    float Xoff, Voff, b_to_a, c_to_a
    float A1, A2, A3
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
    void free_particle_copies() nogil
    void alloc_particle_copies(np.int64_t total_copies) nogil
    void free_halos() nogil

# For outputing halos, rockstar style

cdef import from "meta_io.h":
    void output_halos(np.int64_t id_offset, np.int64_t snap, np.int64_t chunk, float *bounds) nogil

# For setting up the configuration of rockstar

cdef import from "config.h":
    void setup_config() nogil
    void output_config(char *fn) nogil

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
    
    cdef public object ds
    cdef public object fof

    # For future use/consistency
    def __cinit__(self,ds):
        self.ds = ds

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
        
        ds = self.ds

        h0 = ds.hubble_constant
        Ol = ds.omega_lambda
        Om = ds.omega_matter
        
        SCALE_NOW = 1.0/(ds.current_redshift+1.0)
        
        if not outbase =='None'.decode('UTF-8'):
            #output directory. since we can't change the output filenames
            #workaround is to make a new directory
            OUTBASE = outbase 


        PARTICLE_MASS = particle_mass.in_units('Msun/h')
        PERIODIC = periodic
        BOX_SIZE = ds.domain_width.in_units('Mpccm/h')[0]

        # Set up the configuration options
        setup_config()

        # Needs to be called so rockstar can use the particle mass parameter
        # to calculate virial quantities properly
        calc_mass_definition()

    def output_halos(self):
        output_halos(0, 0, 0, NULL) 

    def return_halos(self):
        cdef haloflat[:] haloview = <haloflat[:num_halos]> (<haloflat*> halos)
        rv = np.asarray(haloview).copy()
        rockstar_cleanup()
        free_halos()
        return rv

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def make_rockstar_fof(self, np.ndarray[np.int64_t, ndim=1] pind,
                                np.ndarray[np.int64_t, ndim=1] fof_tags,
                                np.ndarray[anyfloat, ndim=2] pos,
                                np.ndarray[anyfloat, ndim=2] vel):

        # Define fof object

        # Find number of particles
        cdef np.int64_t i, j, k, ind, offset
        cdef np.int64_t num_particles = pind.shape[0]
        global global_particles

        # Allocate space for correct number of particles
        cdef fof fof_obj

        cdef np.int64_t max_count = 0
        cdef np.int64_t next_tag, local_tag, last_fof_tag = -1
        fof_obj.num_p = 0
        j = 0
        # We're going to do one iteration to get the most frequent value.
        for i in range(pind.shape[0]):
            ind = pind[i]
            local_tag = fof_tags[ind]
            # Don't count the null group
            if local_tag == -1: continue
            if local_tag != last_fof_tag:
                if j > max_count:
                    max_count = j
                last_fof_tag = local_tag
                j = 1
            else:
                j += 1
        if j > max_count:
            max_count = j
        #print >> sys.stderr, "Most frequent occurrance: %s" % max_count
        fof_obj.particles = <particle*> malloc(max_count * sizeof(particle))
        j = 0
        cdef int counter = 0, ndone = 0
        cdef np.ndarray[np.int64_t, ndim=1] pcounts 
        pcounts = np.zeros(np.unique(fof_tags).size, dtype="int64")
        cdef np.int64_t frac = <np.int64_t> (pcounts.shape[0] / 20.0)
        free_halos()
        for i in range(pind.shape[0]):
            ind = pind[i]
            local_tag = fof_tags[ind]
            # Skip this one -- it means no group.
            if local_tag == -1:
                continue
            if i == pind.shape[0] - 1:
                next_tag = local_tag + 1
            else:
                next_tag = fof_tags[pind[i+1]]
            for k in range(3):
                fof_obj.particles[j].pos[k] = pos[ind,k]
                fof_obj.particles[j].pos[k+3] = vel[ind,k]
            fof_obj.particles[j].id = j
            fof_obj.num_p += 1
            j += 1
            # Now we check if we're the last one
            if local_tag != next_tag:
                pcounts[ndone] = fof_obj.num_p
                counter += 1
                ndone += 1
                if counter == frac:
                    print >> sys.stderr, "R*-ing % 5.1f%% done (%0.3f -> %0.3f)" % (
                        (100.0 * ndone)/pcounts.size,
                        fof_obj.particles[0].pos[2],
                        halos[num_halos - 1].pos[2])
                    counter = 0
                global_particles = &fof_obj.particles[0]
                find_subs(&fof_obj)
                # Now we reset
                fof_obj.num_p = j = 0
        free(fof_obj.particles)
        return pcounts
