cimport cython
import numpy as np
cimport numpy as np
import sys

from yt.geometry.selection_routines cimport \
    SelectorObject, AlwaysSelector, OctreeSubsetSelector
from yt.utilities.lib.fp_utils cimport imax
from yt.geometry.oct_container cimport \
    SparseOctreeContainer
from yt.geometry.oct_visitors cimport Oct
from yt.geometry.particle_deposit cimport \
    ParticleDepositOperation
from libc.stdlib cimport malloc, free
from libc.string cimport memcpy
import data_structures
from yt.utilities.lib.misc_utilities import OnceIndirect

cdef extern from "platform_dep.h":
    ctypedef int int32_t
    ctypedef long long int64_t
    void *alloca(int)

cdef extern from "cosmology.h":
    ctypedef struct CosmologyParameters "CosmologyParameters" :
        pass
    CosmologyParameters *cosmology_allocate()
    void cosmology_free(CosmologyParameters *c)
    void cosmology_set_fixed(CosmologyParameters *c)

    void cosmology_set_OmegaM(CosmologyParameters *c, double value)
    void cosmology_set_OmegaB(CosmologyParameters *c, double value)
    void cosmology_set_OmegaL(CosmologyParameters *c, double value)
    void cosmology_set_h(CosmologyParameters *c, double value)
    void cosmology_set_DeltaDC(CosmologyParameters *c, double value)

    double abox_from_auni(CosmologyParameters *c, double a) nogil
    double tcode_from_auni(CosmologyParameters *c, double a) nogil
    double tphys_from_auni(CosmologyParameters *c, double a) nogil

    double auni_from_abox(CosmologyParameters *c, double v) nogil
    double auni_from_tcode(CosmologyParameters *c, double v) nogil
    double auni_from_tphys(CosmologyParameters *c, double v) nogil

    double abox_from_tcode(CosmologyParameters *c, double tcode) nogil
    double tcode_from_abox(CosmologyParameters *c, double abox) nogil

    double tphys_from_abox(CosmologyParameters *c, double abox) nogil
    double tphys_from_tcode(CosmologyParameters *c, double tcode) nogil

cdef extern from "artio.h":
    ctypedef struct artio_fileset_handle "artio_fileset" :
        pass
    ctypedef struct artio_selection "artio_selection" :
        pass
    ctypedef struct artio_context :
        pass
    cdef extern artio_context *artio_context_global

    # open modes
    cdef int ARTIO_OPEN_HEADER "ARTIO_OPEN_HEADER"
    cdef int ARTIO_OPEN_GRID "ARTIO_OPEN_GRID"
    cdef int ARTIO_OPEN_PARTICLES "ARTIO_OPEN_PARTICLES"

    # parameter constants
    cdef int ARTIO_TYPE_STRING "ARTIO_TYPE_STRING"
    cdef int ARTIO_TYPE_CHAR "ARTIO_TYPE_CHAR"
    cdef int ARTIO_TYPE_INT "ARTIO_TYPE_INT"
    cdef int ARTIO_TYPE_FLOAT "ARTIO_TYPE_FLOAT"
    cdef int ARTIO_TYPE_DOUBLE "ARTIO_TYPE_DOUBLE"
    cdef int ARTIO_TYPE_LONG "ARTIO_TYPE_LONG"

    cdef int ARTIO_MAX_STRING_LENGTH "ARTIO_MAX_STRING_LENGTH"

    cdef int ARTIO_PARAMETER_EXHAUSTED "ARTIO_PARAMETER_EXHAUSTED"

    # grid read options
    cdef int ARTIO_READ_LEAFS "ARTIO_READ_LEAFS"
    cdef int ARTIO_READ_REFINED "ARTIO_READ_REFINED"
    cdef int ARTIO_READ_ALL "ARTIO_READ_ALL"
    cdef int ARTIO_READ_REFINED_NOT_ROOT "ARTIO_READ_REFINED_NOT_ROOT"
    cdef int ARTIO_RETURN_CELLS "ARTIO_RETURN_CELLS"
    cdef int ARTIO_RETURN_OCTS "ARTIO_RETURN_OCTS"

    # errors
    cdef int ARTIO_SUCCESS "ARTIO_SUCCESS"
    cdef int ARTIO_ERR_MEMORY_ALLOCATION "ARTIO_ERR_MEMORY_ALLOCATION"

    artio_fileset_handle *artio_fileset_open(char *file_prefix, int type, artio_context *context )
    int artio_fileset_close( artio_fileset_handle *handle )
    int artio_fileset_open_particle( artio_fileset_handle *handle )
    int artio_fileset_open_grid(artio_fileset_handle *handle)
    int artio_fileset_close_grid(artio_fileset_handle *handle)

    int artio_fileset_has_grid( artio_fileset_handle *handle )
    int artio_fileset_has_particles( artio_fileset_handle *handle )

    # selection functions
    artio_selection *artio_selection_allocate( artio_fileset_handle *handle )
    artio_selection *artio_select_all( artio_fileset_handle *handle )
    artio_selection *artio_select_volume( artio_fileset_handle *handle, double lpos[3], double rpos[3] )
    int artio_selection_add_root_cell( artio_selection *selection, int coords[3] )
    int artio_selection_destroy( artio_selection *selection )
    int artio_selection_iterator( artio_selection *selection,
            int64_t max_range_size, int64_t *start, int64_t *end )
    int64_t artio_selection_size( artio_selection *selection )
    void artio_selection_print( artio_selection *selection )

    # parameter functions
    int artio_parameter_iterate( artio_fileset_handle *handle, char *key, int *type, int *length )
    int artio_parameter_get_int_array(artio_fileset_handle *handle, char * key, int length, int32_t *values)
    int artio_parameter_get_float_array(artio_fileset_handle *handle, char * key, int length, float *values)
    int artio_parameter_get_long_array(artio_fileset_handle *handle, char * key, int length, int64_t *values)
    int artio_parameter_get_double_array(artio_fileset_handle *handle, char * key, int length, double *values)
    int artio_parameter_get_string_array(artio_fileset_handle *handle, char * key, int length, char **values )

    # grid functions
    int artio_grid_cache_sfc_range(artio_fileset_handle *handle, int64_t start, int64_t end)
    int artio_grid_clear_sfc_cache( artio_fileset_handle *handle )

    int artio_grid_read_root_cell_begin(artio_fileset_handle *handle, int64_t sfc,
        double *pos, float *variables,
        int *num_tree_levels, int *num_octs_per_level)
    int artio_grid_read_root_cell_end(artio_fileset_handle *handle)

    int artio_grid_read_level_begin(artio_fileset_handle *handle, int level )
    int artio_grid_read_level_end(artio_fileset_handle *handle)

    int artio_grid_read_oct(artio_fileset_handle *handle, double *pos,
            float *variables, int *refined)

    int artio_grid_count_octs_in_sfc_range(artio_fileset_handle *handle,
            int64_t start, int64_t end, int64_t *num_octs)

    #particle functions
    int artio_fileset_open_particles(artio_fileset_handle *handle)
    int artio_particle_read_root_cell_begin(artio_fileset_handle *handle, int64_t sfc,
                        int * num_particle_per_species)
    int artio_particle_read_root_cell_end(artio_fileset_handle *handle)
    int artio_particle_read_particle(artio_fileset_handle *handle, int64_t *pid, int *subspecies,
                        double *primary_variables, float *secondary_variables)
    int artio_particle_cache_sfc_range(artio_fileset_handle *handle, int64_t sfc_start, int64_t sfc_end)
    int artio_particle_clear_sfc_cache(artio_fileset_handle *handle)
    int artio_particle_read_species_begin(artio_fileset_handle *handle, int species)
    int artio_particle_read_species_end(artio_fileset_handle *handle)


cdef extern from "artio_internal.h":
    np.int64_t artio_sfc_index( artio_fileset_handle *handle, int coords[3] ) nogil
    void artio_sfc_coords( artio_fileset_handle *handle, int64_t index, int coords[3] ) nogil

cdef void check_artio_status(int status, char *fname="[unknown]"):
    if status != ARTIO_SUCCESS:
        import traceback
        traceback.print_stack()
        callername = sys._getframe().f_code.co_name
        nline = sys._getframe().f_lineno
        raise RuntimeError('failure with status', status, 'in function',fname,'from caller', callername, nline)

cdef class artio_fileset :
    cdef public object parameters
    cdef artio_fileset_handle *handle

    # cosmology parameter for time unit conversion
    cdef CosmologyParameters *cosmology
    cdef float tcode_to_years

    # common attributes
    cdef public int num_grid
    cdef int64_t num_root_cells
    cdef int64_t sfc_min, sfc_max

    # grid attributes
    cdef public int has_grid
    cdef public int min_level, max_level
    cdef public int num_grid_variables
    cdef int *num_octs_per_level
    cdef float *grid_variables

    # particle attributes
    cdef public int has_particles
    cdef public int num_species
    cdef int *particle_position_index
    cdef int *num_particles_per_species
    cdef double *primary_variables
    cdef float *secondary_variables

    def __init__(self, char *file_prefix) :
        cdef int artio_type = ARTIO_OPEN_HEADER
        cdef int64_t num_root

        self.handle = artio_fileset_open( file_prefix, artio_type, artio_context_global )
        if not self.handle :
            raise RuntimeError

        self.read_parameters()

        self.num_root_cells = self.parameters['num_root_cells'][0]
        self.num_grid = 1
        num_root = self.num_root_cells
        while num_root > 1 :
            self.num_grid <<= 1
            num_root >>= 3

        self.sfc_min = 0
        self.sfc_max = self.num_root_cells-1

        # initialize cosmology module
        if "abox" in self.parameters:
            self.cosmology = cosmology_allocate()
            cosmology_set_OmegaM(self.cosmology, self.parameters['OmegaM'][0])
            cosmology_set_OmegaL(self.cosmology, self.parameters['OmegaL'][0])
            cosmology_set_OmegaB(self.cosmology, self.parameters['OmegaB'][0])
            cosmology_set_h(self.cosmology, self.parameters['hubble'][0])
            cosmology_set_DeltaDC(self.cosmology, self.parameters['DeltaDC'][0])
            cosmology_set_fixed(self.cosmology)
        else:
            self.cosmology = NULL
            self.tcode_to_years = self.parameters['time_unit'][0]/(365.25*86400)

        # grid detection
        self.min_level = 0
        self.max_level = self.parameters['grid_max_level'][0]
        self.num_grid_variables = self.parameters['num_grid_variables'][0]

        self.num_octs_per_level = <int *>malloc(self.max_level*sizeof(int))
        self.grid_variables = <float *>malloc(8*self.num_grid_variables*sizeof(float))
        if (not self.num_octs_per_level) or (not self.grid_variables) :
            raise MemoryError

        if artio_fileset_has_grid(self.handle):
            status = artio_fileset_open_grid(self.handle)
            check_artio_status(status)
            self.has_grid = 1
        else:
            self.has_grid = 0

        # particle detection
        if ( artio_fileset_has_particles(self.handle) ):
            status = artio_fileset_open_particles(self.handle)
            check_artio_status(status)
            self.has_particles = 1
	    
            for v in ["num_particle_species","num_primary_variables","num_secondary_variables"]:
                if v not in self.parameters:
                    raise RuntimeError("Unable to locate particle header information in artio header: key=", v)

            self.num_species = self.parameters['num_particle_species'][0]
            self.particle_position_index = <int *>malloc(3*sizeof(int)*self.num_species)
            if not self.particle_position_index :
                raise MemoryError
            for ispec in range(self.num_species) :
                species_labels = "species_%02d_primary_variable_labels"% (ispec,)
                if species_labels not in self.parameters:
                    raise RuntimeError("Unable to locate variable labels for species",ispec)

                labels = self.parameters[species_labels]
                try :
                    self.particle_position_index[3*ispec+0] = labels.index('POSITION_X')
                    self.particle_position_index[3*ispec+1] = labels.index('POSITION_Y')
                    self.particle_position_index[3*ispec+2] = labels.index('POSITION_Z')
                except ValueError :
                    raise RuntimeError("Unable to locate position information for particle species", ispec)

            self.num_particles_per_species =  <int *>malloc(sizeof(int)*self.num_species)
            self.primary_variables = <double *>malloc(sizeof(double)*max(self.parameters['num_primary_variables']))
            self.secondary_variables = <float *>malloc(sizeof(float)*max(self.parameters['num_secondary_variables']))
            if (not self.num_particles_per_species) or (not self.primary_variables) or (not self.secondary_variables) :
                raise MemoryError
        else:
            self.has_particles = 0

    def __dealloc__(self) :
        if self.num_octs_per_level : free(self.num_octs_per_level)
        if self.grid_variables : free(self.grid_variables)

        if self.particle_position_index : free(self.particle_position_index)
        if self.num_particles_per_species : free(self.num_particles_per_species)
        if self.primary_variables : free(self.primary_variables)
        if self.secondary_variables : free(self.secondary_variables)

        if self.cosmology : cosmology_free(self.cosmology)

        if self.handle : artio_fileset_close(self.handle)

    def read_parameters(self) :
        from sys import version
        cdef char key[64]
        cdef int type
        cdef int length
        cdef char ** char_values
        cdef int32_t *int_values
        cdef int64_t *long_values
        cdef float *float_values
        cdef double *double_values

        self.parameters = {}

        while artio_parameter_iterate( self.handle, key, &type, &length ) == ARTIO_SUCCESS :
	    
            if type == ARTIO_TYPE_STRING :
                char_values = <char **>malloc(length*sizeof(char *))
                for i in range(length) :
                    char_values[i] = <char *>malloc( ARTIO_MAX_STRING_LENGTH*sizeof(char) )
                artio_parameter_get_string_array( self.handle, key, length, char_values )
                parameter = [ char_values[i] for i in range(length) ]
                for i in range(length) :
                    free(char_values[i])
                free(char_values)
                if version[0] == '3':
                    for i in range(0,len(parameter)):
                        parameter[i] = parameter[i].decode('utf-8')
            elif type == ARTIO_TYPE_INT :
                int_values = <int32_t *>malloc(length*sizeof(int32_t))
                artio_parameter_get_int_array( self.handle, key, length, int_values )
                parameter = [ int_values[i] for i in range(length) ]
                free(int_values)
            elif type == ARTIO_TYPE_LONG :
                long_values = <int64_t *>malloc(length*sizeof(int64_t))
                artio_parameter_get_long_array( self.handle, key, length, long_values )
                parameter = [ long_values[i] for i in range(length) ]
                free(long_values)
            elif type == ARTIO_TYPE_FLOAT :
                float_values = <float *>malloc(length*sizeof(float))
                artio_parameter_get_float_array( self.handle, key, length, float_values )
                parameter = [ float_values[i] for i in range(length) ]
                free(float_values)
            elif type == ARTIO_TYPE_DOUBLE :
                double_values = <double *>malloc(length*sizeof(double))
                artio_parameter_get_double_array( self.handle, key, length, double_values )
                parameter = [ double_values[i] for i in range(length) ]
                free(double_values)
            else :
                raise RuntimeError("ARTIO file corruption detected: invalid type!")

            if version[0] == '3':
                self.parameters[key.decode('utf-8')] = parameter
            else:
                self.parameters[key] = parameter

    def abox_from_auni(self, np.float64_t a):
        if self.cosmology:
            return abox_from_auni(self.cosmology, a)
        else:
            raise RuntimeError("abox_from_auni called for non-cosmological ARTIO fileset!")

    def tcode_from_auni(self, np.float64_t a):
        if self.cosmology:
            return tcode_from_auni(self.cosmology, a)
        else:
            raise RuntimeError("tcode_from_auni called for non-cosmological ARTIO fileset!")

    def tphys_from_auni(self, np.float64_t a):
        if self.cosmology:
            return tphys_from_auni(self.cosmology, a)
        else:
            raise RuntimeError("tphys_from_auni called for non-cosmological ARTIO fileset!")

    def auni_from_abox(self, np.float64_t v):
        if self.cosmology:
            return auni_from_abox(self.cosmology, v)
        else:
            raise RuntimeError("auni_from_abox called for non-cosmological ARTIO fileset!")

    def auni_from_tcode(self, np.float64_t v):
        if self.cosmology:
            return auni_from_tcode(self.cosmology, v)
        else:
            raise RuntimeError("auni_from_tcode called for non-cosmological ARTIO fileset!")

    @cython.wraparound(False)
    @cython.boundscheck(False)
    def auni_from_tcode_array(self, np.ndarray[np.float64_t] tcode):
        cdef int i, nauni
        cdef np.ndarray[np.float64_t, ndim=1] auni
        cdef CosmologyParameters *cosmology = self.cosmology
        if not cosmology:
            raise RuntimeError("auni_from_tcode_array called for non-cosmological ARTIO fileset!")
        auni = np.empty_like(tcode)
        nauni = auni.shape[0]
        with nogil:
            for i in range(nauni):
                auni[i] = auni_from_tcode(self.cosmology, tcode[i])
        return auni

    def auni_from_tphys(self, np.float64_t v):
        if self.cosmology:
            return auni_from_tphys(self.cosmology, v)
        else:
            raise RuntimeError("auni_from_tphys called for non-cosmological ARTIO fileset!")

    def abox_from_tcode(self, np.float64_t abox):
        if self.cosmology:
            return abox_from_tcode(self.cosmology, abox)
        else:
            raise RuntimeError("abox_from_tcode called for non-cosmological ARTIO fileset!")

    def tcode_from_abox(self, np.float64_t abox):
        if self.cosmology:
            return tcode_from_abox(self.cosmology, abox)
        else:
            raise RuntimeError("tcode_from_abox called for non-cosmological ARTIO fileset!")

    def tphys_from_abox(self, np.float64_t abox):
        if self.cosmology:
            return tphys_from_abox(self.cosmology, abox)
        else:
            raise RuntimeError("tphys_from_abox called for non-cosmological ARTIO fileset!")

    def tphys_from_tcode(self, np.float64_t tcode):
        if self.cosmology:
            return tphys_from_tcode(self.cosmology, tcode)
        else:
            return self.tcode_to_years*tcode

    @cython.wraparound(False)
    @cython.boundscheck(False)
    def tphys_from_tcode_array(self, np.ndarray[np.float64_t] tcode):
        cdef int i, ntphys
        cdef np.ndarray[np.float64_t, ndim=1] tphys
        cdef CosmologyParameters *cosmology = self.cosmology
        tphys = np.empty_like(tcode)
        ntphys = tcode.shape[0]

        if cosmology:
            tphys = np.empty_like(tcode)
            ntphys = tcode.shape[0]
            with nogil:
                for i in range(ntphys):
                    tphys[i] = tphys_from_tcode(cosmology, tcode[i])
            return tphys
        else:
            return tcode*self.tcode_to_years

#    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def read_particle_chunk(self, SelectorObject selector, int64_t sfc_start, int64_t sfc_end, fields) :
        cdef int i
        cdef int status
        cdef int subspecies
        cdef int64_t pid

        cdef int num_fields = len(fields)
        cdef np.float64_t pos[3]

        # since RuntimeErrors are not fatal, ensure no artio_particles* functions
        # called if fileset lacks particles
        if not self.has_particles: return

        data = {}
        accessed_species = np.zeros( self.num_species, dtype="int")
        selected_mass = [ None for i in range(self.num_species)]
        selected_pid = [ None for i in range(self.num_species)]
        selected_species = [ None for i in range(self.num_species)]
        selected_primary = [ [] for i in range(self.num_species)]
        selected_secondary = [ [] for i in range(self.num_species)]

        for species,field in fields :
            if species < 0 or species > self.num_species :
                raise RuntimeError("Invalid species provided to read_particle_chunk")
            accessed_species[species] = 1

            if self.parameters["num_primary_variables"][species] > 0 and \
                    field in self.parameters["species_%02u_primary_variable_labels"%(species,)] :
                selected_primary[species].append((self.parameters["species_%02u_primary_variable_labels"%(species,)].index(field),(species,field)))
                data[(species,field)] = np.empty(0,dtype="float64")
            elif self.parameters["num_secondary_variables"][species] > 0 and \
                    field in self.parameters["species_%02u_secondary_variable_labels"%(species,)] :
                selected_secondary[species].append((self.parameters["species_%02u_secondary_variable_labels"%(species,)].index(field),(species,field)))
                data[(species,field)] = np.empty(0,dtype="float64")
            elif field == "MASS" :
                selected_mass[species] = (species,field)
                data[(species,field)] = np.empty(0,dtype="float64")
            elif field == "PID" :
                selected_pid[species] = (species,field)
                data[(species,field)] = np.empty(0,dtype="int64")
            elif field == "SPECIES" :
                selected_species[species] = (species,field)
                data[(species,field)] = np.empty(0,dtype="int8")
            else :
                raise RuntimeError("invalid field name provided to read_particle_chunk")

        # cache the range
        status = artio_particle_cache_sfc_range( self.handle, self.sfc_min, self.sfc_max )
        check_artio_status(status)

        for sfc in range( sfc_start, sfc_end+1 ) :
            status = artio_particle_read_root_cell_begin( self.handle, sfc,
                    self.num_particles_per_species )
            check_artio_status(status)

            for ispec in range(self.num_species) :
                if accessed_species[ispec] :
                    status = artio_particle_read_species_begin(self.handle, ispec)
                    check_artio_status(status)

                    for particle in range( self.num_particles_per_species[ispec] ) :
                        status = artio_particle_read_particle(self.handle,
                                &pid, &subspecies, self.primary_variables,
                                self.secondary_variables)
                        check_artio_status(status)

                        for i in range(3) :
                            pos[i] = self.primary_variables[self.particle_position_index[3*ispec+i]]

                        if selector.select_point(pos) :
                            # loop over primary variables
                            for i,field in selected_primary[ispec] :
                                count = len(data[field])
                                data[field].resize(count+1)
                                data[field][count] = self.primary_variables[i]

                            # loop over secondary variables
                            for i,field in selected_secondary[ispec] :
                                count = len(data[field])
                                data[field].resize(count+1)
                                data[field][count] = self.secondary_variables[i]

                            # add particle id
                            if selected_pid[ispec] :
                                count = len(data[selected_pid[ispec]])
                                data[selected_pid[ispec]].resize(count+1)
                                data[selected_pid[ispec]][count] = pid

                            # add mass if requested
                            if selected_mass[ispec] :
                                count = len(data[selected_mass[ispec]])
                                data[selected_mass[ispec]].resize(count+1)
                                data[selected_mass[ispec]][count] = self.parameters["particle_species_mass"][ispec]

                    status = artio_particle_read_species_end( self.handle )
                    check_artio_status(status)

            status = artio_particle_read_root_cell_end( self.handle )
            check_artio_status(status)

        return data

    #@cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def read_grid_chunk(self, SelectorObject selector, int64_t sfc_start, int64_t sfc_end, fields ):
        cdef int i
        cdef int level
        cdef int num_oct_levels
        cdef int refined[8]
        cdef int status
        cdef int64_t count
        cdef int64_t max_octs
        cdef double dpos[3]
        cdef np.float64_t left[3]
        cdef np.float64_t right[3]
        cdef np.float64_t dds[3]

        cdef int *field_order
        cdef int num_fields  = len(fields)
        field_order = <int*>malloc(sizeof(int)*num_fields)

        # translate fields from ARTIO names to indices
        var_labels = self.parameters['grid_variable_labels']
        for i, f in enumerate(fields):
            if f not in var_labels:
                raise RuntimeError("Field",f,"is not known to ARTIO")
            field_order[i] = var_labels.index(f)

        status = artio_grid_cache_sfc_range( self.handle, self.sfc_min, self.sfc_max )
        check_artio_status(status)

        # determine max number of cells we could hit (optimize later)
        #status = artio_grid_count_octs_in_sfc_range( self.handle,
        #        sfc_start, sfc_end, &max_octs )
        #check_artio_status(status)
        #max_cells = sfc_end-sfc_start+1 + max_octs*8

        # allocate space for _fcoords, _icoords, _fwidth, _ires
        #fcoords = np.empty((max_cells, 3), dtype="float64")
        #ires = np.empty(max_cells, dtype="int64")
        fcoords = np.empty((0, 3), dtype="float64")
        ires = np.empty(0, dtype="int64")

        #data = [ np.empty(max_cells, dtype="float32") for i in range(num_fields) ]
        data = [ np.empty(0,dtype="float64") for i in range(num_fields)]

        count = 0
        for sfc in range( sfc_start, sfc_end+1 ) :
            status = artio_grid_read_root_cell_begin( self.handle, sfc,
                    dpos, self.grid_variables, &num_oct_levels, self.num_octs_per_level )
            check_artio_status(status)

            if num_oct_levels == 0 :
                for i in range(num_fields) :
                    data[i].resize(count+1)
                    data[i][count] = self.grid_variables[field_order[i]]
                fcoords.resize((count+1,3))
                for i in range(3) :
                    fcoords[count][i] = dpos[i]
                ires.resize(count+1)
                ires[count] = 0
                count += 1

            for level in range(1,num_oct_levels+1) :
                status = artio_grid_read_level_begin( self.handle, level )
                check_artio_status(status)

                for i in range(3) :
                    dds[i] = 2.**-level

                for oct in range(self.num_octs_per_level[level-1]) :
                    status = artio_grid_read_oct( self.handle, dpos, self.grid_variables, refined )
                    check_artio_status(status)

                    for child in range(8) :
                        if not refined[child] :
                            for i in range(3) :
                                left[i] = (dpos[i]-dds[i]) if (child & (i<<1)) else dpos[i]
                                right[i] = left[i] + dds[i]

                            if selector.select_bbox(left,right) :
                                fcoords.resize((count+1, 3))
                                for i in range(3) :
                                    fcoords[count][i] = left[i]+0.5*dds[i]
                                ires.resize(count+1)
                                ires[count] = level
                                for i in range(num_fields) :
                                    data[i].resize(count+1)
                                    data[i][count] = self.grid_variables[self.num_grid_variables*child+field_order[i]]
                                count += 1
                status = artio_grid_read_level_end( self.handle )
                check_artio_status(status)

            status = artio_grid_read_root_cell_end( self.handle )
            check_artio_status(status)

        free(field_order)

        #fcoords.resize((count,3))
        #ires.resize(count)
        #
        #for i in range(num_fields) :
        #    data[i].resize(count)

        return (fcoords, ires, data)

    def root_sfc_ranges_all(self, int max_range_size = 1024) :
        cdef int64_t sfc_start, sfc_end
        cdef artio_selection *selection

        selection = artio_select_all( self.handle )
        if selection == NULL :
            raise RuntimeError
        sfc_ranges = []
        while artio_selection_iterator(selection, max_range_size,
                &sfc_start, &sfc_end) == ARTIO_SUCCESS :
            sfc_ranges.append([sfc_start, sfc_end])
        artio_selection_destroy(selection)
        return sfc_ranges

    def root_sfc_ranges(self, SelectorObject selector,
                        int max_range_size = 1024):
        cdef int coords[3]
        cdef int64_t sfc_start, sfc_end
        cdef np.float64_t left[3]
        cdef np.float64_t right[3]
        cdef np.float64_t dds[3]
        cdef artio_selection *selection
        cdef int i, j, k

        sfc_ranges=[]
        selection = artio_selection_allocate(self.handle)
        for i in range(self.num_grid) :
            # stupid cython
            coords[0] = i
            left[0] = coords[0]
            right[0] = left[0] + 1.0
            for j in range(self.num_grid) :
                coords[1] = j
                left[1] = coords[1]
                right[1] = left[1] + 1.0
                for k in range(self.num_grid) :
                    coords[2] = k
                    left[2] = coords[2]
                    right[2] = left[2] + 1.0
                    if selector.select_bbox(left,right) :
                        status = artio_selection_add_root_cell(selection, coords)
                        check_artio_status(status)

        while artio_selection_iterator(selection, max_range_size,
                &sfc_start, &sfc_end) == ARTIO_SUCCESS :
            sfc_ranges.append([sfc_start, sfc_end])

        artio_selection_destroy(selection)
        return sfc_ranges

###################################################
def artio_is_valid( char *file_prefix ) :
    cdef artio_fileset_handle *handle = artio_fileset_open( file_prefix,
            ARTIO_OPEN_HEADER, artio_context_global )
    if handle == NULL :
        return False
    else :
        artio_fileset_close(handle)
    return True

cdef class ARTIOSFCRangeHandler:
    cdef public np.int64_t sfc_start
    cdef public np.int64_t sfc_end
    cdef public artio_fileset artio_handle
    cdef public object root_mesh_handler
    cdef public object oct_count
    cdef public object octree_handler
    cdef artio_fileset_handle *handle
    cdef np.float64_t DLE[3]
    cdef np.float64_t DRE[3]
    cdef np.float64_t dds[3]
    cdef np.int64_t dims[3]
    cdef public np.int64_t total_octs
    cdef np.int64_t *doct_count
    cdef np.int64_t **pcount
    cdef float **root_mesh_data
    cdef np.int64_t nvars[2]
    cdef int cache_root_mesh

    def __init__(self, domain_dimensions, # cells
                 domain_left_edge,
                 domain_right_edge,
                 artio_fileset artio_handle,
                 sfc_start, sfc_end, int cache_root_mesh = 0):
        cdef int i
        cdef np.int64_t sfc
        self.sfc_start = sfc_start
        self.sfc_end = sfc_end
        self.artio_handle = artio_handle
        self.root_mesh_handler = None
        self.octree_handler = None
        self.handle = artio_handle.handle
        self.oct_count = None
        self.root_mesh_data = NULL
        self.pcount = NULL
        self.cache_root_mesh = cache_root_mesh

        if artio_handle.has_particles:
            self.pcount = <np.int64_t **> malloc(sizeof(np.int64_t*)
                * artio_handle.num_species)
            self.nvars[0] = artio_handle.num_species
            for i in range(artio_handle.num_species):
                self.pcount[i] = <np.int64_t*> malloc(sizeof(np.int64_t)
                    * (self.sfc_end - self.sfc_start + 1))
                for sfc in range(self.sfc_end - self.sfc_start + 1):
                    self.pcount[i][sfc] = 0
        else:
            self.nvars[0] = 0

        if artio_handle.has_grid:
            self.nvars[1] = artio_handle.num_grid_variables
        else:
            self.nvars[1] = 0

        for i in range(3):
            self.dims[i] = domain_dimensions[i]
            self.DLE[i] = domain_left_edge[i]
            self.DRE[i] = domain_right_edge[i]
            self.dds[i] = (self.DRE[i] - self.DLE[i])/self.dims[i]

    def __dealloc__(self):
        cdef int i
        if self.artio_handle.has_particles:
            for i in range(self.nvars[0]):
                free(self.pcount[i])
            free(self.pcount)
        if self.artio_handle.has_grid:
            if self.root_mesh_data != NULL:
                for i in range(self.nvars[1]):
                    free(self.root_mesh_data[i])
                free(self.root_mesh_data)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def construct_mesh(self):
        cdef int status, level, ngv
        cdef np.int64_t sfc, oc, i
        cdef double dpos[3]
        cdef int num_oct_levels
        cdef int max_level = self.artio_handle.max_level
        cdef int *num_octs_per_level = <int *>malloc(
            (max_level + 1)*sizeof(int))
        cdef int num_species = self.artio_handle.num_species
        cdef int *num_particles_per_species
        cdef ARTIOOctreeContainer octree
        ngv = self.nvars[1]
        cdef float *grid_variables = <float *>malloc(
            ngv * sizeof(float))
        self.octree_handler = octree = ARTIOOctreeContainer(self)
        if self.cache_root_mesh == 1:
            self.root_mesh_data = <float **>malloc(sizeof(float *) * ngv)
            for i in range(ngv):
                self.root_mesh_data[i] = <float *>malloc(sizeof(float) * \
                    (self.sfc_end - self.sfc_start + 1))
        # We want to pre-allocate an array of root pointers.  In the future,
        # this will be pre-determined by the ARTIO library.  However, because
        # realloc plays havoc with our tree searching, we can't utilize an
        # expanding array at the present time.
        octree.allocate_domains([], self.sfc_end - self.sfc_start + 1)
        cdef np.ndarray[np.int64_t, ndim=1] oct_count
        oct_count = np.zeros(self.sfc_end - self.sfc_start + 1, dtype="int64")
        status = artio_grid_cache_sfc_range(self.handle, self.sfc_start,
                                            self.sfc_end)
        check_artio_status(status)
        for sfc in range(self.sfc_start, self.sfc_end + 1):
            status = artio_grid_read_root_cell_begin( self.handle,
                sfc, dpos, grid_variables, &num_oct_levels,
                num_octs_per_level)
            check_artio_status(status)
            for i in range(ngv * self.cache_root_mesh):
                self.root_mesh_data[i][sfc - self.sfc_start] = \
                    grid_variables[i]
            if num_oct_levels > 0:
                oc = 0
                for level in range(num_oct_levels):
                    oc += num_octs_per_level[level]
                self.total_octs += oc
                oct_count[sfc - self.sfc_start] = oc
                octree.initialize_local_mesh(oc, num_oct_levels,
                    num_octs_per_level, sfc)
            status = artio_grid_read_root_cell_end(self.handle)
            check_artio_status(status)
        status = artio_grid_clear_sfc_cache(self.handle)
        check_artio_status(status)
        if self.artio_handle.has_particles:
            num_particles_per_species =  <int *>malloc(
                    sizeof(int)*num_species)

            # Now particles
            status = artio_particle_cache_sfc_range(self.handle, self.sfc_start,
                                            self.sfc_end)
            check_artio_status(status)
            for sfc in range(self.sfc_start, self.sfc_end + 1):
                # Now particles
                status = artio_particle_read_root_cell_begin(self.handle,
                        sfc, num_particles_per_species)
                check_artio_status(status)

                for i in range(num_species):
                    self.pcount[i][sfc - self.sfc_start] = \
                        num_particles_per_species[i]

                status = artio_particle_read_root_cell_end(self.handle)
                check_artio_status(status)

            status = artio_particle_clear_sfc_cache(self.handle)
            check_artio_status(status)

            free(num_particles_per_species)

        free(grid_variables)
        free(num_octs_per_level)
        self.oct_count = oct_count
        self.doct_count = <np.int64_t *> oct_count.data
        self.root_mesh_handler = ARTIORootMeshContainer(self)

    def free_mesh(self):
        self.octree_handler = None
        self.root_mesh_handler = None
        self.doct_count = NULL
        self.oct_count = None

def get_coords(artio_fileset handle, np.int64_t s):
    cdef int coords[3]
    artio_sfc_coords(handle.handle, s, coords)
    return (coords[0], coords[1], coords[2])

cdef struct particle_var_pointers:
    # The number of particles we have filled
    np.int64_t count
    # Number of primary variables and pointers to their indices
    int n_p
    int p_ind[16] # Max of 16 vars
    # Number of secondary variables and pointers to their indices
    int n_s
    int s_ind[16] # Max of 16 vars
    # Pointers to the bools and data arrays for mass, pid and species
    int n_mass
    np.float64_t *mass
    int n_pid
    np.int64_t *pid
    int n_species
    np.int8_t *species
    # Pointers to the pointers to primary and secondary vars
    np.float64_t *pvars[16]
    np.float64_t *svars[16]

cdef class ARTIOOctreeContainer(SparseOctreeContainer):
    # This is a transitory, created-on-demand OctreeContainer.  It should not
    # be considered to be long-lasting, and during its creation it will read
    # the index file.  This means that when created it will then be able to
    # provide coordinates, but then on subsequent IO accesses it will pass over
    # the file again, despite knowing the indexing system already.  Because of
    # this, we will avoid creating it as long as possible.

    cdef public artio_fileset artio_handle
    cdef ARTIOSFCRangeHandler range_handler
    cdef np.int64_t level_indices[32]
    cdef np.int64_t sfc_start
    cdef np.int64_t sfc_end

    def __init__(self, ARTIOSFCRangeHandler range_handler):
        self.range_handler = range_handler
        self.sfc_start = range_handler.sfc_start
        self.sfc_end = range_handler.sfc_end
        self.artio_handle = range_handler.artio_handle
        # Note the final argument is partial_coverage, which indicates whether
        # or not an Oct can be partially refined.
        dims, DLE, DRE = [], [], []
        for i in range(32):
            self.level_indices[i] = 0
        for i in range(3):
            # range_handler has dims in cells, which is the same as the number
            # of possible octs.  This is because we have a forest of octrees.
            dims.append(range_handler.dims[i])
            DLE.append(range_handler.DLE[i])
            DRE.append(range_handler.DRE[i])
        super(ARTIOOctreeContainer, self).__init__(dims, DLE, DRE)
        self.artio_handle = range_handler.artio_handle
        self.level_offset = 1
        self.domains = NULL
        self.root_nodes = NULL

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void initialize_local_mesh(self, np.int64_t oct_count,
                              int num_oct_levels, int *num_octs_per_level,
                              np.int64_t sfc):
        # We actually will not be initializing the root mesh here, we will be
        # initializing the entire mesh between sfc_start and sfc_end.
        cdef np.int64_t oct_ind, tot_octs, ipos, nadded
        cdef int i, status, level, num_root, num_octs
        cdef int num_level_octs
        cdef artio_fileset_handle *handle = self.artio_handle.handle
        cdef int coords[3]
        cdef int max_level = self.artio_handle.max_level
        cdef double dpos[3]
        cdef np.float64_t f64pos[3]
        cdef np.float64_t dds[3]
        cdef np.ndarray[np.float64_t, ndim=2] pos
        # NOTE: We do not cache any SFC ranges here, as we should only ever be
        # called from within a pre-cached operation in the SFC handler.

        # We only allow one root oct.
        self.append_domain(oct_count)
        self.domains[self.num_domains - 1].con_id = sfc

        oct_ind = -1
        ipos = 0
        for level in range(num_oct_levels):
            oct_ind = imax(oct_ind, num_octs_per_level[level])
            self.level_indices[level] = ipos
            ipos += num_octs_per_level[level]
        pos = np.empty((oct_ind, 3), dtype="float64")

        # Now we initialize
        # Note that we also assume we have already started reading the level.
        ipos = 0
        for level in range(num_oct_levels):
            status = artio_grid_read_level_begin(handle, level + 1)
            check_artio_status(status)
            for oct_ind in range(num_octs_per_level[level]):
                status = artio_grid_read_oct(handle, dpos, NULL, NULL)
                for i in range(3):
                    pos[oct_ind, i] = dpos[i]
                check_artio_status(status)
            status = artio_grid_read_level_end(handle)
            check_artio_status(status)
            nadded = self.add(self.num_domains, level, pos[:num_octs_per_level[level],:])
            if nadded != num_octs_per_level[level]:
                raise RuntimeError

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fill_sfc(self,
             np.ndarray[np.uint8_t, ndim=1] levels,
             np.ndarray[np.uint8_t, ndim=1] cell_inds,
             np.ndarray[np.int64_t, ndim=1] file_inds,
             np.ndarray[np.int64_t, ndim=1] domain_counts,
             field_indices, dest_fields):
        cdef np.ndarray[np.float32_t, ndim=2] source
        cdef np.ndarray[np.float64_t, ndim=1] dest
        cdef int n, status, i, di, num_oct_levels, nf, ngv, max_level
        cdef int level, j, oct_ind, si
        cdef np.int64_t sfc, ipos
        cdef np.float64_t val
        cdef artio_fileset_handle *handle = self.artio_handle.handle
        cdef double dpos[3]
        # We duplicate some of the grid_variables stuff here so that we can
        # potentially release the GIL
        nf = len(field_indices)
        ngv = self.artio_handle.num_grid_variables
        max_level = self.artio_handle.max_level
        cdef int *num_octs_per_level = <int *>malloc(
            (max_level + 1)*sizeof(int))
        cdef float *grid_variables = <float *>malloc(
            8 * ngv * sizeof(float))
        cdef int* field_ind = <int*> malloc(
            nf * sizeof(int))
        cdef np.float32_t **field_vals = <np.float32_t**> malloc(
            nf * sizeof(np.float32_t*))
        source_arrays = []
        ipos = -1
        for i in range(self.num_domains):
            ipos = imax(ipos, self.domains[i].n)
        for i in range(nf):
            field_ind[i] = field_indices[i]
            # Note that we subtract one, because we're not using the root mesh.
            source = np.zeros((ipos, 8), dtype="float32")
            source_arrays.append(source)
            field_vals[i] = <np.float32_t*> source.data
        cdef np.int64_t level_position[32]
        cdef np.int64_t lp
        # First we need to walk the mesh in the file.  Then we fill in the dest
        # location based on the file index.
        # A few ways this could be improved:
        #   * Create a new visitor function that actually queried the data,
        #     rather than our somewhat hokey double-loop over SFC arrays.
        #   * Go to pointers for the dest arrays.
        #   * Enable preloading during mesh initialization
        #   * Calculate domain indices on the fly rather than with a
        #     double-loop to calculate domain_counts
        # The cons should be in order
        cdef np.int64_t sfc_start, sfc_end
        sfc_start = self.domains[0].con_id
        sfc_end = self.domains[self.num_domains - 1].con_id
        status = artio_grid_cache_sfc_range(handle, sfc_start, sfc_end)
        check_artio_status(status)
        cdef np.int64_t offset = 0
        for si in range(self.num_domains):
            sfc = self.domains[si].con_id
            status = artio_grid_read_root_cell_begin( handle, sfc,
                    dpos, NULL, &num_oct_levels, num_octs_per_level)
            check_artio_status(status)
            lp = 0
            for level in range(num_oct_levels):
                status = artio_grid_read_level_begin(handle, level + 1)
                check_artio_status(status)
                level_position[level] = lp
                for oct_ind in range(num_octs_per_level[level]):
                    status = artio_grid_read_oct(handle, dpos, grid_variables, NULL)
                    check_artio_status(status)
                    for j in range(8):
                        for i in range(nf):
                            field_vals[i][(oct_ind + lp)*8+j] = \
                                grid_variables[ngv*j+field_ind[i]]
                status = artio_grid_read_level_end(handle)
                check_artio_status(status)
                lp += num_octs_per_level[level]
            status = artio_grid_read_root_cell_end(handle)
            check_artio_status(status)
            # Now we have all our sources.
            for j in range(nf):
                dest = dest_fields[j]
                source = source_arrays[j]
                for i in range(domain_counts[si]):
                    level = levels[i + offset]
                    oct_ind = file_inds[i + offset] + level_position[level]
                    dest[i + offset] = source[oct_ind, cell_inds[i + offset]]
            # Now, we offset by the actual number filled here.
            offset += domain_counts[si]
        status = artio_grid_clear_sfc_cache(handle)
        check_artio_status(status)
        free(field_ind)
        free(field_vals)
        free(grid_variables)
        free(num_octs_per_level)

    def fill_sfc_particles(self, fields):
        # This handles not getting particles for refined sfc values.
        rv = read_sfc_particles(self.artio_handle, self.sfc_start,
                                self.sfc_end, 0, fields,
                                self.range_handler.doct_count,
                                self.range_handler.pcount)
        return rv

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef read_sfc_particles(artio_fileset artio_handle,
                        np.int64_t sfc_start, np.int64_t sfc_end,
                        int read_unrefined, fields,
                        np.int64_t *doct_count,
                        np.int64_t **pcount):
    cdef int status, ispec, subspecies
    cdef np.int64_t sfc, particle, pid, ind, vind
    cdef int num_species = artio_handle.num_species
    cdef artio_fileset_handle *handle = artio_handle.handle
    cdef int num_oct_levels
    cdef int *num_particles_per_species =  <int *>malloc(
        sizeof(int)*num_species)
    cdef int *accessed_species =  <int *>malloc(
        sizeof(int)*num_species)
    cdef int *total_particles = <int *>malloc(
        sizeof(int)*num_species)
    cdef particle_var_pointers *vpoints = <particle_var_pointers *> malloc(
        sizeof(particle_var_pointers)*num_species)
    cdef double *primary_variables
    cdef double dpos[3]
    cdef float *secondary_variables
    cdef np.int64_t tp
    cdef int max_level = artio_handle.max_level
    cdef int *num_octs_per_level = <int *>malloc(
        (max_level + 1)*sizeof(int))

    cdef np.ndarray[np.int8_t, ndim=1] npi8arr
    cdef np.ndarray[np.int64_t, ndim=1] npi64arr
    cdef np.ndarray[np.float64_t, ndim=1] npf64arr

    if not artio_handle.has_particles:
        raise RuntimeError("Attempted to read non-existent particles in ARTIO")

    # Now we set up our field pointers
    params = artio_handle.parameters

    npri_vars = params["num_primary_variables"]
    nsec_vars = params["num_secondary_variables"]
    primary_variables = <double *>malloc(sizeof(double) * max(npri_vars))
    secondary_variables = <float *>malloc(sizeof(float) * max(nsec_vars))

    cdef particle_var_pointers *vp

    for ispec in range(num_species):
        total_particles[ispec] = 0
        accessed_species[ispec] = 0
        # Initialize our vpoints array
        vpoints[ispec].count = 0
        vpoints[ispec].n_mass = 0
        vpoints[ispec].n_pid = 0
        vpoints[ispec].n_species = 0
        vpoints[ispec].n_p = 0
        vpoints[ispec].n_s = 0

    # Pass through once.  We want every single particle.
    tp = 0
    cdef np.int64_t c
    for sfc in range(sfc_start, sfc_end + 1):
        c = doct_count[sfc - sfc_start]
        if read_unrefined == 1 and c > 0: continue
        if read_unrefined == 0 and c == 0: continue

        for ispec in range(num_species):
            total_particles[ispec] += pcount[ispec][sfc - sfc_start]

    # Now we allocate our final fields, which will be filled
    #for ispec in range(num_species):
    #    print "In SFC %s to %s reading %s of species %s" % (
    #        sfc_start, sfc_end + 1, total_particles[ispec], ispec)
    data = {}
    for species, field in sorted(fields):
        accessed_species[species] = 1
        pri_vars = params.get(
            "species_%02u_primary_variable_labels" % (species,), [])
        sec_vars = params.get(
            "species_%02u_secondary_variable_labels" % (species,), [])
        tp = total_particles[species]
        vp = &vpoints[species]
        if field == "PID":
            vp.n_pid = 1
            data[(species, field)] = np.zeros(tp, dtype="int64")
            npi64arr = data[(species, field)]
            vp.pid = <np.int64_t*> npi64arr.data
        elif field == "SPECIES":
            vp.n_species = 1
            data[(species, field)] = np.zeros(tp, dtype="int8")
            npi8arr = data[(species, field)]
            # We fill this *now*
            npi8arr += species
            vp.species = <np.int8_t*> npi8arr.data
        elif npri_vars[species] > 0 and field in pri_vars :
            data[(species, field)] = np.zeros(tp, dtype="float64")
            npf64arr = data[(species, field)]
            vp.p_ind[vp.n_p] = pri_vars.index(field)
            vp.pvars[vp.n_p] = <np.float64_t *> npf64arr.data
            vp.n_p += 1
        elif nsec_vars[species] > 0 and field in sec_vars :
            data[(species, field)] = np.zeros(tp, dtype="float64")
            npf64arr = data[(species, field)]
            vp.s_ind[vp.n_s] = sec_vars.index(field)
            vp.svars[vp.n_s] = <np.float64_t *> npf64arr.data
            vp.n_s += 1
        elif field == "MASS":
            vp.n_mass = 1
            data[(species, field)] = np.zeros(tp, dtype="float64")
            npf64arr = data[(species, field)]
            # We fill this *now*
            npf64arr += params["particle_species_mass"][species]
            vp.mass = <np.float64_t*> npf64arr.data

    status = artio_particle_cache_sfc_range( handle,
            sfc_start, sfc_end )
    check_artio_status(status)

    for sfc in range(sfc_start, sfc_end + 1):
        c = doct_count[sfc - sfc_start]
        check_artio_status(status)
        if read_unrefined == 1 and c > 0: continue
        if read_unrefined == 0 and c == 0: continue
        c = 0
        for ispec in range(num_species) :
            if accessed_species[ispec] == 0: continue
            c += pcount[ispec][sfc - sfc_start]
        if c == 0: continue
        status = artio_particle_read_root_cell_begin( handle, sfc,
                num_particles_per_species )
        check_artio_status(status)
        for ispec in range(num_species) :
            if accessed_species[ispec] == 0: continue
            status = artio_particle_read_species_begin(handle, ispec);
            check_artio_status(status)
            vp = &vpoints[ispec]

            for particle in range(num_particles_per_species[ispec]) :
                status = artio_particle_read_particle(handle,
                        &pid, &subspecies, primary_variables,
                        secondary_variables)
                check_artio_status(status)
                ind = vp.count

                for i in range(vp.n_p):
                    vind = vp.p_ind[i]
                    vp.pvars[i][ind] = primary_variables[vind]

                for i in range(vp.n_s):
                    vind = vp.s_ind[i]
                    vp.svars[i][ind] = secondary_variables[vind]

                if vp.n_pid:
                    vp.pid[ind] = pid

                vp.count += 1

            status = artio_particle_read_species_end( handle )
            check_artio_status(status)

        status = artio_particle_read_root_cell_end( handle )
        check_artio_status(status)

    status = artio_particle_clear_sfc_cache(handle)
    check_artio_status(status)

    free(num_octs_per_level)
    free(num_particles_per_species)
    free(total_particles)
    free(accessed_species)
    free(vpoints)
    free(primary_variables)
    free(secondary_variables)
    return data

cdef class ARTIORootMeshContainer:
    cdef public artio_fileset artio_handle
    cdef np.float64_t DLE[3]
    cdef np.float64_t DRE[3]
    cdef np.float64_t dds[3]
    cdef np.int64_t dims[3]
    cdef artio_fileset_handle *handle
    cdef np.uint64_t sfc_start
    cdef np.uint64_t sfc_end
    cdef public object _last_mask
    cdef public np.int64_t _last_selector_id
    cdef np.int64_t _last_mask_sum
    cdef ARTIOSFCRangeHandler range_handler
    cdef np.uint8_t *sfc_mask
    cdef np.int64_t nsfc

    def __init__(self, ARTIOSFCRangeHandler range_handler):
        cdef int i
        cdef np.int64_t sfci
        for i in range(3):
            self.DLE[i] = range_handler.DLE[i]
            self.DRE[i] = range_handler.DRE[i]
            self.dims[i] = range_handler.dims[i]
            self.dds[i] = range_handler.dds[i]
        self.handle = range_handler.handle
        self.artio_handle = range_handler.artio_handle
        self._last_mask = None
        self._last_selector_id = -1
        self.sfc_start = range_handler.sfc_start
        self.sfc_end = range_handler.sfc_end
        self.range_handler = range_handler
        # We assume that the number of octs has been created and filled
        # already.  We no longer care about ANY of the SFCs that have octs
        # inside them -- this goes for every operation that this object
        # performs.
        self.sfc_mask = <np.uint8_t *>malloc(sizeof(np.uint8_t) *
          self.sfc_end - self.sfc_start + 1)
        self.nsfc = 0
        for sfci in range(self.sfc_end - self.sfc_start + 1):
            if self.range_handler.oct_count[sfci] > 0:
                self.sfc_mask[sfci] = 0
            else:
                self.sfc_mask[sfci] = 1
                self.nsfc += 1

    def __dealloc__(self):
        free(self.sfc_mask)

    @cython.cdivision(True)
    cdef np.int64_t pos_to_sfc(self, np.float64_t pos[3]) nogil:
        # Calculate the index
        cdef int coords[3]
        cdef int i
        cdef np.int64_t sfc
        for i in range(3):
            coords[i] = <int>((pos[i] - self.DLE[i])/self.dds[i])
        sfc = artio_sfc_index(self.handle, coords)
        return sfc

    @cython.cdivision(True)
    cdef void sfc_to_pos(self, np.int64_t sfc, np.float64_t pos[3]) nogil:
        cdef int coords[3]
        cdef int i
        artio_sfc_coords(self.handle, sfc, coords)
        for i in range(3):
            pos[i] = self.DLE[i] + (coords[i] + 0.5) * self.dds[i]

    cdef np.int64_t count_cells(self, SelectorObject selector):
        # We visit each cell if it is not refined and determine whether it is
        # included or not.
        cdef np.int64_t sfc
        cdef np.float64_t pos[3]
        cdef np.float64_t right_edge[3]
        cdef int num_cells = 0
        cdef int i
        return self.mask(selector).sum()

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def icoords(self, SelectorObject selector, np.int64_t num_cells = -1,
                int domain_id = -1):
        # Note that num_octs does not have to equal sfc_end - sfc_start + 1.
        cdef np.int64_t sfc, sfci = -1
        cdef int acoords[3]
        cdef int i
        cdef np.ndarray[np.uint8_t, ndim=1, cast=True] mask
        mask = self.mask(selector)
        num_cells = self._last_mask_sum
        cdef np.ndarray[np.int64_t, ndim=2] coords
        coords = np.empty((num_cells, 3), dtype="int64")
        cdef int filled = 0
        for sfc in range(self.sfc_start, self.sfc_end + 1):
            if self.sfc_mask[sfc - self.sfc_start] == 0: continue
            sfci += 1
            if mask[sfci] == 0: continue
            # Note that we do *no* checks on refinement here.  In fact, this
            # entire setup should not need to touch the disk except if the
            # artio sfc calculators need to.
            artio_sfc_coords(self.handle, sfc, acoords)
            for i in range(3):
                coords[filled, i] = acoords[i]
            filled += 1
        return coords

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fcoords(self, SelectorObject selector, np.int64_t num_cells = -1,
                int domain_id = -1):
        # Note that num_cells does not have to equal sfc_end - sfc_start + 1.
        cdef np.int64_t sfc, sfci = -1
        cdef np.float64_t pos[3]
        cdef int acoords[3]
        cdef int i
        cdef np.ndarray[np.uint8_t, ndim=1, cast=True] mask
        mask = self.mask(selector)
        num_cells = self._last_mask_sum
        cdef np.ndarray[np.float64_t, ndim=2] coords
        coords = np.empty((num_cells, 3), dtype="float64")
        cdef int filled = 0
        for sfc in range(self.sfc_start, self.sfc_end + 1):
            if self.sfc_mask[sfc - self.sfc_start] == 0: continue
            sfci += 1
            if mask[sfci] == 0: continue
            # Note that we do *no* checks on refinement here.  In fact, this
            # entire setup should not need to touch the disk except if the
            # artio sfc calculators need to.
            self.sfc_to_pos(sfc, pos)
            for i in range(3):
                coords[filled, i] = pos[i]
            filled += 1
        return coords

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fwidth(self, SelectorObject selector, np.int64_t num_cells = -1,
                int domain_id = -1):
        cdef int i
        cdef np.ndarray[np.uint8_t, ndim=1, cast=True] mask
        mask = self.mask(selector)
        num_cells = self._last_mask_sum
        cdef np.ndarray[np.float64_t, ndim=2] width
        width = np.zeros((num_cells, 3), dtype="float64")
        for i in range(3):
            width[:,i] = self.dds[i]
        return width

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def ires(self, SelectorObject selector, np.int64_t num_cells = -1,
                int domain_id = -1):
        cdef np.ndarray[np.uint8_t, ndim=1, cast=True] mask
        mask = self.mask(selector)
        num_cells = self._last_mask_sum
        cdef np.ndarray[np.int64_t, ndim=1] res
        res = np.zeros(num_cells, dtype="int64")
        return res

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def selector_fill(self, SelectorObject selector,
                      np.ndarray source,
                      np.ndarray dest = None,
                      np.int64_t offset = 0, int dims = 1,
                      int domain_id = -1):
        # This is where we use the selector to transplant from one to the
        # other.  Note that we *do* apply the selector here.
        cdef np.int64_t num_cells = -1
        cdef np.int64_t ind
        cdef np.int64_t sfc, sfci = -1
        cdef np.float64_t pos[3]
        cdef np.float64_t dpos[3]
        cdef int dim, status, filled = 0
        cdef int num_oct_levels, level
        cdef int max_level = self.artio_handle.max_level
        cdef int *num_octs_per_level = <int *>malloc(
            (max_level + 1)*sizeof(int))
        cdef char *sdata = <char*> source.data
        cdef char *ddata
        cdef int ss = source.dtype.itemsize
        cdef np.ndarray[np.uint8_t, ndim=1, cast=True] mask
        mask = self.mask(selector)
        if dest is None:
            # Note that RAMSES can have partial refinement inside an Oct.  This
            # means we actually do want the number of Octs, not the number of
            # cells.
            num_cells = self._last_mask_sum
            if dims > 1:
                dest = np.zeros((num_cells, dims), dtype=source.dtype,
                    order='C')
            else:
                dest = np.zeros(num_cells, dtype=source.dtype, order='C')
        ddata = (<char*>dest.data) + offset*ss*dims
        ind = 0
        for sfc in range(self.sfc_start, self.sfc_end + 1):
            if self.sfc_mask[sfc - self.sfc_start] == 0: continue
            sfci += 1
            if mask[sfci] == 0: continue
            memcpy(ddata, sdata + ind, dims * ss)
            ddata += dims * ss
            filled += 1
            ind += ss * dims
        if num_cells >= 0:
            return dest
        return filled

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def mask(self, SelectorObject selector, np.int64_t num_cells = -1,
             int domain_id = -1):
        # We take a domain_id here to avoid subclassing
        cdef int i
        cdef np.float64_t pos[3]
        cdef np.int64_t sfc, sfci = -1
        if self._last_selector_id == hash(selector):
            return self._last_mask
        cdef np.ndarray[np.uint8_t, ndim=1] mask
        mask = np.zeros((self.nsfc), dtype="uint8")
        self._last_mask_sum = 0
        for sfc in range(self.sfc_start, self.sfc_end + 1):
            if self.sfc_mask[sfc - self.sfc_start] == 0: continue
            sfci += 1
            self.sfc_to_pos(sfc, pos)
            if selector.select_cell(pos, self.dds) == 0: continue
            mask[sfci] = 1
            self._last_mask_sum += 1
        self._last_mask = mask.astype("bool")
        self._last_selector_id = hash(selector)
        return self._last_mask


    def fill_sfc_particles(self, fields):
        rv = read_sfc_particles(self.artio_handle,
                                self.sfc_start, self.sfc_end,
                                1, fields, self.range_handler.doct_count,
                                self.range_handler.pcount)
        return rv

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fill_sfc(self, SelectorObject selector, field_indices):
        cdef np.ndarray[np.float64_t, ndim=1] dest
        cdef int n, status, i, di, num_oct_levels, nf, ngv, max_level
        cdef np.int64_t sfc, num_cells, sfci = -1
        cdef np.float64_t val
        cdef double dpos[3]
        max_level = self.artio_handle.max_level
        cdef int *num_octs_per_level = <int *>malloc(
            (max_level + 1)*sizeof(int))
        # We duplicate some of the grid_variables stuff here so that we can
        # potentially release the GIL
        nf = len(field_indices)
        ngv = self.artio_handle.num_grid_variables
        cdef float *grid_variables = <float *>malloc(
            ngv * sizeof(float))
        cdef np.ndarray[np.uint8_t, ndim=1, cast=True] mask
        mask = self.mask(selector, -1)
        num_cells = self._last_mask_sum
        tr = []
        for i in range(nf):
            tr.append(np.zeros(num_cells, dtype="float64"))
        cdef int* field_ind = <int*> malloc(
            nf * sizeof(int))
        cdef np.float64_t **field_vals = <np.float64_t**> malloc(
            nf * sizeof(np.float64_t*))
        for i in range(nf):
            field_ind[i] = field_indices[i]
            # This zeros should be an empty once we handle the root grid
            dest = tr[i]
            field_vals[i] = <np.float64_t*> dest.data
        # First we need to walk the mesh in the file.  Then we fill in the dest
        # location based on the file index.
        cdef int filled = 0
        cdef float **mesh_data = self.range_handler.root_mesh_data
        if mesh_data == NULL:
            status = artio_grid_cache_sfc_range(self.handle, self.sfc_start,
                                                self.sfc_end)
            check_artio_status(status)
            for sfc in range(self.sfc_start, self.sfc_end + 1):
                if self.sfc_mask[sfc - self.sfc_start] == 0: continue
                sfci += 1
                if mask[sfci] == 0: continue
                status = artio_grid_read_root_cell_begin( self.handle,
                    sfc, dpos, grid_variables, &num_oct_levels,
                    num_octs_per_level)
                check_artio_status(status)
                for i in range(nf):
                    field_vals[i][filled] = grid_variables[field_ind[i]]
                filled += 1
                status = artio_grid_read_root_cell_end(self.handle)
                check_artio_status(status)
            status = artio_grid_clear_sfc_cache(self.handle)
            check_artio_status(status)
        else:
            for sfc in range(self.sfc_start, self.sfc_end + 1):
                if self.sfc_mask[sfc - self.sfc_start] == 0: continue
                sfci += 1
                if mask[sfci] == 0: continue
                for i in range(nf):
                    field_vals[i][filled] = mesh_data[field_ind[i]][
                        sfc - self.sfc_start]
                filled += 1
        # Now we have all our sources.
        free(field_ind)
        free(field_vals)
        free(grid_variables)
        free(num_octs_per_level)
        return tr

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def deposit(self, ParticleDepositOperation pdeposit,
                SelectorObject selector,
                np.ndarray[np.float64_t, ndim=2] positions,
                fields):
        # This implements the necessary calls to enable particle deposition to
        # occur as needed.
        cdef int nf, i, j
        cdef np.int64_t sfc, sfci
        if fields is None:
            fields = []
        nf = len(fields)
        cdef np.float64_t[::cython.view.indirect, ::1] field_pointers 
        if nf > 0: field_pointers = OnceIndirect(fields)
        cdef np.float64_t[:] field_vals = np.empty(nf, dtype="float64")
        cdef np.ndarray[np.uint8_t, ndim=1, cast=True] mask
        mask = self.mask(selector, -1)
        cdef np.ndarray[np.int64_t, ndim=1] domain_ind
        domain_ind = np.zeros(self.sfc_end - self.sfc_start + 1,
                              dtype="int64") - 1
        j = 0
        sfci = -1
        for sfc in range(self.sfc_start, self.sfc_end + 1):
            if self.sfc_mask[sfc - self.sfc_start] == 0: continue
            sfci += 1
            if mask[sfci] == 0:
                continue
            domain_ind[sfc - self.sfc_start] = j
            j += 1
        cdef np.float64_t pos[3]
        cdef np.float64_t left_edge[3]
        cdef int coords[3]
        cdef int dims[3]
        dims[0] = dims[1] = dims[2] = 1
        cdef np.int64_t offset, moff
        cdef np.int64_t numpart = positions.shape[0]
        for i in range(positions.shape[0]):
            for j in range(nf):
                field_vals[j] = field_pointers[j][i]
            for j in range(3):
                pos[j] = positions[i, j]
                coords[j] = <int>((pos[j] - self.DLE[j])/self.dds[j])
            sfc = artio_sfc_index(self.artio_handle.handle, coords)
            if sfc < self.sfc_start or sfc > self.sfc_end: continue
            offset = domain_ind[sfc - self.sfc_start]
            if offset < 0: continue
            # Check that we found the oct ...
            for j in range(3):
                left_edge[j] = coords[j] * self.dds[j] + self.DLE[j]
            pdeposit.process(dims, left_edge, self.dds,
                         offset, pos, field_vals, sfc)
            if pdeposit.update_values == 1:
                for j in range(nf):
                    field_pointers[j][i] = field_vals[j]

cdef class SFCRangeSelector(SelectorObject):

    cdef SelectorObject base_selector
    cdef ARTIOSFCRangeHandler range_handler
    cdef ARTIORootMeshContainer mesh_container
    cdef np.int64_t sfc_start, sfc_end

    def __init__(self, dobj):
        self.base_selector = dobj.base_selector
        self.mesh_container = dobj.oct_handler
        self.range_handler = self.mesh_container.range_handler
        self.sfc_start = self.mesh_container.sfc_start
        self.sfc_end = self.mesh_container.sfc_end

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def select_grids(self,
                     np.ndarray[np.float64_t, ndim=2] left_edges,
                     np.ndarray[np.float64_t, ndim=2] right_edges,
                     np.ndarray[np.int32_t, ndim=2] levels):
        raise RuntimeError

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        return self.select_point(pos)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_point(self, np.float64_t pos[3]) nogil:
        cdef np.int64_t sfc = self.mesh_container.pos_to_sfc(pos)
        if sfc > self.sfc_end: return 0
        cdef np.int64_t oc = self.range_handler.doct_count[
            sfc - self.sfc_start]
        if oc > 0: return 0
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        return self.base_selector.select_bbox(left_edge, right_edge)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_grid(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3], np.int32_t level,
                         Oct *o = NULL) nogil:
        # Because visitors now use select_grid, we should be explicitly
        # checking this.
        return self.base_selector.select_grid(left_edge, right_edge, level, o)

    def _hash_vals(self):
        return (hash(self.base_selector), self.sfc_start, self.sfc_end)

sfc_subset_selector = AlwaysSelector
#sfc_subset_selector = SFCRangeSelector

