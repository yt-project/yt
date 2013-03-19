"""

"""
cimport cython
import numpy as np
cimport numpy as np
import sys 

from yt.geometry.selection_routines cimport SelectorObject
from libc.stdint cimport int32_t, int64_t
from libc.stdlib cimport malloc, free
import  data_structures  

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
    int artio_particle_read_species_begin(artio_fileset_handle *handle, int species)
    int artio_particle_read_species_end(artio_fileset_handle *handle) 
   

cdef check_artio_status(int status, char *fname="[unknown]"):
    if status!=ARTIO_SUCCESS :
        callername = sys._getframe().f_code.co_name
        nline = sys._getframe().f_lineno
        raise RuntimeError('failure with status', status, 'in function',fname,'from caller', callername, nline)

cdef class artio_fileset :
    cdef public object parameters 
    cdef artio_fileset_handle *handle

    # common attributes
    cdef public int num_grid
    cdef int64_t num_root_cells
    cdef int64_t sfc_min, sfc_max

    # grid attributes
    cdef public int min_level, max_level
    cdef public int num_grid_variables
    cdef int *num_octs_per_level
    cdef float *grid_variables

    # particle attributes
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

        # grid detection
        self.min_level = 0
        self.max_level = self.parameters['grid_max_level'][0]
        self.num_grid_variables = self.parameters['num_grid_variables'][0]

        self.num_octs_per_level = <int *>malloc(self.max_level*sizeof(int))
        self.grid_variables = <float *>malloc(8*self.num_grid_variables*sizeof(float))
        if (not self.num_octs_per_level) or (not self.grid_variables) :
            raise MemoryError

        status = artio_fileset_open_grid( self.handle )
        check_artio_status(status)

        # particle detection
        self.num_species = self.parameters['num_particle_species'][0]
        self.particle_position_index = <int *>malloc(3*sizeof(int)*self.num_species)
        if not self.particle_position_index :
            raise MemoryError
        for ispec in range(self.num_species) :
            labels = self.parameters["species_%02d_primary_variable_labels"% (ispec,)]
            try :
                self.particle_position_index[3*ispec+0] = labels.index('POSITION_X')
                self.particle_position_index[3*ispec+1] = labels.index('POSITION_Y')
                self.particle_position_index[3*ispec+2] = labels.index('POSITION_Z')
            except ValueError :
                raise RuntimeError("Unable to locate position information for particle species", ispec )

        self.num_particles_per_species =  <int *>malloc(sizeof(int)*self.num_species) 
        self.primary_variables = <double *>malloc(sizeof(double)*max(self.parameters['num_primary_variables']))  
        self.secondary_variables = <float *>malloc(sizeof(float)*max(self.parameters['num_secondary_variables']))  
        if (not self.num_particles_per_species) or (not self.primary_variables) or (not self.secondary_variables) :
            raise MemoryError

        status = artio_fileset_open_particles( self.handle )
        check_artio_status(status)
   
    # this should possibly be __dealloc__ 
    def __del__(self) :
        if self.num_octs_per_level : free(self.num_octs_per_level)
        if self.grid_variables : free(self.grid_variables)

        if self.particle_position_index : free(self.particle_position_index)
        if self.num_particles_per_species : free(self.num_particles_per_species)
        if self.primary_variables : free(self.primary_variables)
        if self.secondary_variables : free(self.secondary_variables)

        if self.handle : artio_fileset_close(self.handle)
  
    def read_parameters(self) :
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

            self.parameters[key] = parameter

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
        cdef np.float64_t dds[3]
        cdef int eterm[3]
        for i in range(3) : dds[i] = 0

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
                data[(species,field)] = np.empty(0,dtype="float32")
            elif field == "MASS" :
                selected_mass[species] = (species,field)
                data[(species,field)] = np.empty(0,dtype="float32")
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
                    status = artio_particle_read_species_begin(self.handle, ispec);
                    check_artio_status(status)
 
                    for particle in range( self.num_particles_per_species[ispec] ) :
                        status = artio_particle_read_particle(self.handle,
                                &pid, &subspecies, self.primary_variables,
                                self.secondary_variables)
                        check_artio_status(status)

                        for i in range(3) :
                            pos[i] = self.primary_variables[self.particle_position_index[3*ispec+i]] 

                        if selector.select_cell(pos,dds,eterm) :
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
                                data[selected_mass[ispec]][count] = self.parameters["particle_species_mass"]
                        
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
        cdef np.float64_t pos[3]
        cdef np.float64_t dds[3]
        cdef int eterm[3]
        cdef int *field_order
        cdef int num_fields  = len(fields)
        field_order = <int*>malloc(sizeof(int)*num_fields)

        # translate fields from ARTIO names to indices
        var_labels = self.parameters['grid_variable_labels']
        for i, f in enumerate(fields):
            if f not in var_labels:
                raise RuntimeError("Field",f,"is not known to ARTIO")
            field_order[i] = var_labels.index(f)

        # dhr - cache the entire domain (replace later)
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
        data = [ np.empty(0,dtype="float32") for i in range(num_fields)]

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
                                pos[i] = dpos[i] + dds[i]*(0.5 if (child & (1<<i)) else -0.5)

                            if selector.select_cell( pos, dds, eterm ) :
                                fcoords.resize((count+1, 3))
                                for i in range(3) :
                                    fcoords[count][i] = pos[i]
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

    def root_sfc_ranges_all(self) :
        cdef int max_range_size = 1024
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

    def root_sfc_ranges(self, SelectorObject selector) :
        cdef int max_range_size = 1024
        cdef int coords[3]
        cdef int64_t sfc_start, sfc_end
        cdef np.float64_t pos[3]
        cdef np.float64_t dds[3]
        cdef int eterm[3]
        cdef artio_selection *selection
        cdef int i, j, k
        for i in range(3): dds[i] = 1.0

        sfc_ranges=[]
        selection = artio_selection_allocate(self.handle)
        for i in range(self.num_grid) :
            # stupid cython
            coords[0] = i
            pos[0] = coords[0] + 0.5
            for j in range(self.num_grid) :
                coords[1] = j
                pos[1] = coords[1] + 0.5
                for k in range(self.num_grid) :
                    coords[2] = k 
                    pos[2] = coords[2] + 0.5
                    if selector.select_cell(pos, dds, eterm) :
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
