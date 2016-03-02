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
        float m, r, child_r, vmax_r, mgrav, vmax, rvmax, rs, klypin_rs, vrms
        float J[3]
        float energy, spin
        float alt_m[4]
        float Xoff, Voff, b_to_a, c_to_a
        float A[3]
        float b_to_a2, c_to_a2
        float A2[3]
        float bullock_spin, kin_to_pot, m_pe_b, m_pe_d
        np.int64_t num_p, num_child_particles, p_start, desc, flags, n_core
        float min_pos_err, min_vel_err, min_bulkvel_err, _pad

ctypedef struct haloflat:
    np.int64_t id
    float pos_x, pos_y, pos_z, vel_x, vel_y, vel_z
    float corevel_x, corevel_y, corevel_z
    float bulkvel_x, bulkvel_y, bulkvel_z
    float m, r, child_r, vmax_r, mgrav, vmax, rvmax, rs, klypin_rs, vrms
    float Jx, Jy, Jz
    float energy, spin
    float alt_m1, alt_m2, alt_m3, alt_m4
    float Xoff, Voff, b_to_a, c_to_a
    float Ax, Ay, Az
    float b_to_a2, c_to_a2, A2x, A2y, A2z
    float bullock_spin, kin_to_pot, m_pe_b, m_pe_d
    np.int64_t num_p, num_child_particles, p_start, desc, flags, n_core
    float min_pos_err, min_vel_err, min_bulkvel_err, _pad

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
    float max_halo_radius(halo *h) nogil

# global in groupies.c
cdef extern double particle_thresh_dens[5]

# For outputing halos, rockstar style

cdef import from "meta_io.h":
    void output_halos(np.int64_t id_offset, np.int64_t snap, np.int64_t chunk, float *bounds) nogil

# For setting up the configuration of rockstar

cdef import from "config.h":
    void setup_config() nogil
    void output_config(char *fn) nogil

cdef import from "distance.h":
    void init_cosmology() nogil

cdef import from "config_vars.h":
    # Rockstar cleverly puts all of the config variables inside a templated
    # definition of their vaiables.
    char *FILE_FORMAT
    np.float64_t PARTICLE_MASS

    char *MASS_DEFINITION
    char *MASS_DEFINITION2
    char *MASS_DEFINITION3
    char *MASS_DEFINITION4
    char *MASS_DEFINITION5
    np.int64_t STRICT_SO_MASSES
    np.int64_t MIN_HALO_OUTPUT_SIZE
    np.float64_t FORCE_RES
    np.float64_t FORCE_RES_PHYS_MAX

    np.float64_t SCALE_NOW
    np.float64_t h0
    np.float64_t Ol
    np.float64_t Om
    np.float64_t W0
    np.float64_t WA

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
    np.int64_t RESTART_SNAP
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

    np.int64_t SHAPE_ITERATIONS
    np.int64_t WEIGHTED_SHAPES
    np.int64_t BOUND_PROPS
    np.int64_t BOUND_OUT_TO_HALO_EDGE
    np.int64_t DO_MERGER_TREE_ONLY
    np.int64_t IGNORE_PARTICLE_IDS
    np.float64_t EXACT_LL_CALC
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
    np.int64_t ART_VARIANT

    np.float64_t FOF_FRACTION
    np.float64_t FOF_LINKING_LENGTH
    np.float64_t INITIAL_METRIC_SCALING
    np.float64_t INCLUDE_HOST_POTENTIAL_RATIO
    np.int64_t TEMPORAL_HALO_FINDING
    np.int64_t MIN_HALO_PARTICLES
    np.float64_t UNBOUND_THRESHOLD
    np.int64_t ALT_NFW_METRIC
    np.int64_t EXTRA_PROFILING

    np.int64_t TOTAL_PARTICLES
    np.float64_t BOX_SIZE
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
                       int periodic = 1, force_res = None,
                       int min_halo_size = 25, outbase = "None",
                       write_config = False,  exact_ll_calc = False,
                       lightcone = False, lightcone_origin = [0,0,0],
                       callbacks = None, unbound_threshold=None):
        global FILENAME, FILE_FORMAT, NUM_SNAPS, STARTING_SNAP, h0, Ol, Om
        global BOX_SIZE, PERIODIC, PARTICLE_MASS, NUM_BLOCKS, NUM_READERS
        global FORK_READERS_FROM_WRITERS, PARALLEL_IO_WRITER_PORT, NUM_WRITERS
        global rh, SCALE_NOW, OUTBASE, MIN_HALO_OUTPUT_SIZE
        global OVERLAP_LENGTH, TOTAL_PARTICLES, FORCE_RES
        global OUTPUT_FORMAT, EXTRA_PROFILING
        global STRICT_SO_MASSES, EXACT_LL_CALC
        global LIGHTCONE, LIGHTCONE_ORIGIN

        if force_res is not None:
            FORCE_RES=np.float64(force_res)

        OVERLAP_LENGTH = 0.0
        
        # Set to 0.0 if you plan on calculating spherical overdensity masses.
        # Otherwise filtering of halos in rockstar meta_io.c _should_print
        # will filter the wrong halos when halo mass is re-calculated before
        # output_halos
        global UNBOUND_THRESHOLD
        if unbound_threshold is not None:
            UNBOUND_THRESHOLD = unbound_threshold
        
        FILENAME = "inline.<block>"
        FILE_FORMAT = "GENERIC"
        OUTPUT_FORMAT = "BOTH"
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

        if exact_ll_calc: EXACT_LL_CALC = 1
        STRICT_SO_MASSES = 1    # presumably unused in our code path
        EXTRA_PROFILING = 0

        if lightcone:
            LIGHTCONE = 1
            LIGHTCONE_ORIGIN[0] = lightcone_origin[0]
            LIGHTCONE_ORIGIN[1] = lightcone_origin[1]
            LIGHTCONE_ORIGIN[2] = lightcone_origin[2]

        # Set up the configuration options
        setup_config()

        # Needs to be called so rockstar can use the particle mass parameter
        # to calculate virial quantities properly
        init_cosmology()
        calc_mass_definition()

        if write_config: output_config(NULL)

    def particle_thresh_dens(self):
        cdef np.ndarray d = np.array([particle_thresh_dens[0],
                                      particle_thresh_dens[1],
                                      particle_thresh_dens[2],
                                      particle_thresh_dens[3],
                                      particle_thresh_dens[4]],
                                     dtype=np.float64)
        return d

    def assign_masses(self, h, np.ndarray[np.float32_t, ndim=1] r, float force_res, \
                      double pmass, np.ndarray[np.float64_t, ndim=1] dens_thresh,
                      early_termination=False):
        """
        Assign spherical overdensity masses to halos.  r must be sorted

        Parameters
        ----------
        h: struct haloflat
            Assign masses to this halo
        r: np.ndarray
            Sorted array of particle radii
        force_res: float
            Force resolution, below which density is smoothed.
        dens_thresh: np.ndarray
            Thresholds for spherical overdensity mass calculation
        early_termination: bool
            Specifies whether or not to terminate mass calculation when
            first particle density is below the lowest density threshold.
            If False, may lead to overestimate of SO masses for subhalos,
            but gives a better comparison to plain rockstar masses with
            STRICT_SO=1. Default: False
        Returns
        -------
        None
        """
        cdef double total_mass = 0.0
        cdef double m = 0.0
        cdef double alt_m1 = 0.0
        cdef double alt_m2 = 0.0
        cdef double alt_m3 = 0.0
        cdef double alt_m4 = 0.0
        cdef double rr
        cdef double cur_dens
        cdef int min_ind = np.argmin(dens_thresh)
        cdef int eterm = early_termination
        for rr in r:
            if rr < force_res: rr = force_res
            total_mass += pmass
            cur_dens = total_mass/(rr*rr*rr)
            if cur_dens > dens_thresh[0]: m = total_mass
            if cur_dens > dens_thresh[1]: alt_m1 = total_mass
            if cur_dens > dens_thresh[2]: alt_m2 = total_mass
            if cur_dens > dens_thresh[3]: alt_m3 = total_mass
            if cur_dens > dens_thresh[4]: alt_m4 = total_mass
            if eterm and cur_dens <= dens_thresh[min_ind]:
                break
        h['m'] = m
        h['alt_m1'] = alt_m1
        h['alt_m2'] = alt_m2
        h['alt_m3'] = alt_m3
        h['alt_m4'] = alt_m4
        # if cur_dens > dens_thresh[1]:
            # This is usually a subhalo problem, and we don't know who is a subhalo
            # print >> sys.stderr, "r too small in assign_masses, m200b will be wrong!"
            # print >> sys.stderr, "edge_dens/dens_thresh[1] %.3f" % (cur_dens/dens_thresh[1])

    def max_halo_radius(self, int i):
        return max_halo_radius(&halos[i])

    def output_halos(self, np.int64_t idoffset, np.ndarray[np.float32_t, ndim=2] bbox):
        cdef float bounds[6]
        if idoffset is None: idoffset = 0
        if bbox is None:
            output_halos(idoffset, 0, 0, NULL) 
        else:
            for i in range(3):
                bounds[i] = bbox[i,0]
                bounds[i+3] = bbox[i,1]
            output_halos(idoffset, 0, 0, bounds) 

    def output_config(self):
        output_config(NULL) 

    def return_halos(self):
        cdef haloflat[:] haloview = <haloflat[:num_halos]> (<haloflat*> halos)
        return np.asarray(haloview)

    def finish(self):
        rockstar_cleanup()
        free_halos()

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def make_rockstar_fof(self, np.ndarray[np.int64_t, ndim=1] pind,
                                np.ndarray[np.int64_t, ndim=1] fof_tags,
                                np.ndarray[anyfloat, ndim=2] pos,
                                np.ndarray[anyfloat, ndim=2] vel):

        verbose = False
        # Define fof object

        # Find number of particles
        cdef np.int64_t i, j, k, ind
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
        #print >> sys.stderr, "Most frequent occurrence: %s" % max_count
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
                if verbose and counter == frac:
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
        global_particles = NULL
        return pcounts
