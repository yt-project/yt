"""

"""
cimport cython
import numpy as np
cimport numpy as np
import sys 

from libc.stdint cimport int32_t, int64_t
from libc.stdlib cimport malloc, free
import  data_structures  
from yt.geometry.oct_container cimport \
    OctreeContainer, \
    ARTIOOctreeContainer, \
    Oct

cdef extern from "stdlib.h":
    void *alloca(int)

cdef extern from "artio.h":
    ctypedef struct artio_fileset_handle "artio_fileset" :
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

    # parameter functions
    int artio_parameter_iterate( artio_fileset_handle *handle, char *key, int *type, int *length )
    int artio_parameter_get_int_array(artio_fileset_handle *handle, char * key, int length, int32_t *values)
    int artio_parameter_get_float_array(artio_fileset_handle *handle, char * key, int length, float *values)
    int artio_parameter_get_long_array(artio_fileset_handle *handle, char * key, int length, int64_t *values)
    int artio_parameter_get_double_array(artio_fileset_handle *handle, char * key, int length, double *values)
    int artio_parameter_get_string_array(artio_fileset_handle *handle, char * key, int length, char **values, int max_length)

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
        print 'failure with status', status, 'in function',fname,'from caller', callername, nline 
        sys.exit(1)

cdef class artio_fileset :
    cdef public object parameters 
    cdef ARTIOOctreeContainer oct_handler
    cdef artio_fileset_handle *handle
    cdef int64_t num_root_cells
    cdef int64_t sfc_min, sfc_max
    cdef public int num_grid

    # grid attributes
    cdef int min_level, max_level
    cdef int num_grid_variables
    cdef int num_particle_variables

    cdef int cpu
    cdef int domain_id

    def __init__(self, char *file_prefix) :
        cdef int artio_type = ARTIO_OPEN_HEADER
        cdef int64_t num_root

        self.handle = artio_fileset_open( file_prefix, artio_type, artio_context_global ) 
        self.read_parameters()
        print 'print parameters in caller.pyx',self.parameters
        print 'done reading header parameters'

        self.num_root_cells = self.parameters['num_root_cells'][0]
        self.num_grid = 1
        num_root = self.num_root_cells
        while num_root > 1 :
            self.num_grid <<= 1
            num_root >>= 3
 
        # dhr - add grid detection code 
        status = artio_fileset_open_grid( self.handle )
        check_artio_status(status)
  

        # grid stuff
        self.oct_handler=None

        self.min_level = 0
        self.max_level = self.parameters['grid_max_level'][0]

        #snl FIX: the sfc values used should come from "subset" and describe the domain for chunking
        # note the root level method may force chunking to be done on 0-level ytocts 
        self.sfc_min = 0
        self.sfc_max = self.parameters['grid_file_sfc_index'][1]-1
        self.num_grid_variables = self.parameters['num_grid_variables'][0]

        # these should be fixed to something meaningful
        self.cpu = 0
        self.domain_id = 0

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
                    char_values[i] = <char *>malloc( 128*sizeof(char) )
                artio_parameter_get_string_array( self.handle, key, length, char_values, 128 ) 
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
                print "ERROR: invalid type!"

            self.parameters[key] = parameter

    def count_refined_octs(self) :
        cdef int64_t num_total_octs = 0

        # this only works if this domain includes all cells!
        if self.parameters.has_key("num_octs_per_level") :
            return self.parameters["num_octs_per_level"].sum()

        status = artio_grid_count_octs_in_sfc_range( self.handle, 
                self.sfc_min, self.sfc_max, &num_total_octs )
        check_artio_status(status) 

        # add octs for root cells
        num_total_octs += (self.sfc_max-self.sfc_min+1)/8
 
        return num_total_octs

    def grid_pos_fill(self, ARTIOOctreeContainer oct_handler) :
        ''' adds all refined octs and a new array of ghost octs for  
        the "fully refined octs" at level=-1 in ART or 0 in yt convention 
        so that root level can consist of children
        '''
        cdef int64_t sfc
        cdef int level, ix, iy, iz
        cdef int num_oct_levels
        cdef int *num_octs_per_level
        cdef double dpos[3]
        cdef int refined[8]
        cdef int oct_count
        cdef int ioct

        cdef Oct **next_level_parents=NULL, **cur_level_parents=NULL, **tmp_parents
        cdef int num_next_parents, num_cur_parents, tmp_size
        cdef int next_parents_size, cur_parents_size, cur_parent

        print 'start filling oct positions'
        self.oct_handler = oct_handler
        oct_count = 0

        for iz in range(self.num_grid/2) :
            dpos[2] = iz*2+1
            for iy in range(self.num_grid/2) :
                dpos[1] = iy*2+1
                for ix in range(self.num_grid/2) :
                    dpos[0]=ix*2+1
                    oct_handler.add_oct(self.cpu+1, NULL, 0, dpos)
                    oct_count += 1

        print "done filling root oct positions"

        status = artio_grid_cache_sfc_range( self.handle, self.sfc_min, self.sfc_max )
        check_artio_status(status) 

        num_octs_per_level = <int *>malloc(self.max_level*sizeof(int))
        if not num_octs_per_level :
            raise MemoryError
        next_level_parents = <Oct **>malloc(sizeof(Oct *))
        if not next_level_parents :
            raise MemoryError
        next_parents_size = 1
        cur_level_parents = <Oct **>malloc(sizeof(Oct *))
        if not cur_level_parents :
            raise MemoryError
        cur_parents_size = 1

        for sfc in range( self.sfc_min, self.sfc_max+1 ) :
            status = artio_grid_read_root_cell_begin( self.handle, sfc, 
                dpos, NULL, &num_oct_levels, num_octs_per_level )
            check_artio_status(status) 

            next_level_parents[0] = oct_handler.get_root_oct(dpos)
            num_next_parents = 1

            for level in range(1,num_oct_levels+1) :
                tmp_parents = cur_level_parents
                tmp_size = cur_parents_size

                cur_level_parents = next_level_parents
                cur_parents_size = next_parents_size
                num_cur_parents = num_next_parents

                cur_parent = 0
                num_next_parents = 0
                next_parents_size = tmp_size
                next_level_parents = tmp_parents

                if level < num_oct_levels and \
                         next_parents_size < num_octs_per_level[level] :
                    free(next_level_parents)
                    next_parents_size = num_octs_per_level[level]
                    next_level_parents = <Oct **>malloc(next_parents_size*sizeof(Oct *))
                    if not next_level_parents :
                        raise MemoryError

                status = artio_grid_read_level_begin( self.handle, level )
                check_artio_status(status) 

                for ioct in range(num_octs_per_level[level-1]) :
                    status = artio_grid_read_oct( self.handle, dpos, NULL, refined )
                    check_artio_status(status) 

                    new_oct = oct_handler.add_oct(self.cpu+1, 
                            cur_level_parents[cur_parent],
                            level, dpos)
                    oct_count += 1
                    cur_parent += 1

                    if level < num_oct_levels :
                        for i in range(8) :
                            if refined[i] :
                                next_level_parents[num_next_parents] = new_oct
                                num_next_parents += 1

                status = artio_grid_read_level_end( self.handle )
                check_artio_status(status) 

            status = artio_grid_read_root_cell_end( self.handle )
            check_artio_status(status) 
       
        status = artio_grid_clear_sfc_cache( self.handle )
        check_artio_status(status)

        free(num_octs_per_level)
        free(next_level_parents)
        free(cur_level_parents)

        print 'done filling oct positions', oct_count


    def particle_var_fill(self, accessed_species, source, selector, fields, static_fields) :
        cdef double *primary_variables
        cdef float *secondary_variables
        cdef int *num_particles_per_species
        cdef int status
        cdef np.ndarray[np.float32_t, ndim=1] arr
        # This relies on the fields being contiguous
        cdef np.float32_t **fpoint
        cdef int nf = len(fields)
        fpoint = <np.float32_t**>malloc(sizeof(np.float32_t*)*nf)
        forder = <int*>malloc(sizeof(int)*nf)
        cdef int i, j, level
        cdef int subspecies
        cdef int64_t pid = 0l
#        cdef np.float64_t dds[3], pos[3]
        dds = []
        pos = []
        eterm = []
        dds[0] = dds[1] = dds[2] = 0.0
#        cdef int eterm[3]
        cdef int **mask
        mask = <int**>malloc(sizeof(int*)*len(accessed_species))

        status = artio_particle_cache_sfc_range( self.handle, self.sfc_min, self.sfc_max )
        check_artio_status(status)

        count_mask = []
        count = []
	
        for species in range(len(accessed_species)) :
             mask[species] = <int*>malloc(self.parameters['particle_species_num'][species] * sizeof(int))
             count_mask[species] = 0
             count[species] = 0             

        for sfc in range( self.sfc_min, self.sfc_max+1 ) :
                status = artio_particl`e_read_root_cell_begin( self.handle, sfc,
                    num_particles_per_species )
                check_artio_status(status)

                for species in range(len(accessed_species) ) :
                        status = artio_particle_read_species_begin(self.handle, species)
                        check_artio_status(status)

                        for particle in range( num_particles_per_species[species] ) :
                               status = artio_particle_read_particle(self.handle,
                                        &pid, &subspecies, primary_variables,
                                        secondary_variables)
                               check_artio_status(status)
                               pos[0] = primary_variables[0]
                               pos[1] = primary_variables[1]
                               pos[2] = primary_variables[2]

                               mask[species][count[species]] = selector.select_cell(pos, dds, eterm)
                               count_mask[species] += mask[species][count_mask[species]]
                               count[species] += 1
	
                        status = artio_particle_read_species_end( self.handle )
                        check_artio_status(status)

                status = artio_particle_read_root_cell_end( self.handle )
                check_artio_status(status)

        var_labels = self.parameters['species_00_primary_variable_labels']
        for i, f in enumerate(fields):
            # It might be better to do this check in the Python code
            if f not in var_labels:
		if f not in static_fields :
                    print "This field is not known to ARTIO: ", f
                    raise RuntimeError
            j = var_labels.index(f)
            source[f] = np.empty(count_mask.sum())    
            arr = source[f]
            fpoint[i] = <np.float32_t *>arr.data
            forder[i] = j
            primary_variables = <double *>malloc(8*len(var_labels)*sizeof(float))

        for species in range(len(accessed_species)) :
             count_mask[species] = 0
             count[species] = 0
	
        count_all_particles = 0
        for sfc in range( self.sfc_min, self.sfc_max+1 ) :
                status = artio_particle_read_root_cell_begin( self.handle, sfc,
                    num_particles_per_species )
                check_artio_status(status)	

                for species in range(len(accessed_species) ) :
                    status = artio_particle_read_species_begin(self.handle, species);
                    check_artio_status(status)

                    for particle in range( num_particles_per_species[species] ) :
                        if mask[species][count[species]] == 1 :

                             status = artio_particle_read_particle(self.handle,
                                        &pid, &subspecies, primary_variables,
                                        secondary_variables)
                             check_artio_status(status)

                             for i in range(nf):
                                 fpoint[i][count_all_particles] = primary_variables[forder[i]]
                             count_all_particles += 1
                        count[species] += 1

                    status = artio_particle_read_species_end( self.handle )
                    check_artio_status(status)

                status = artio_particle_read_root_cell_end( self.handle )
                check_artio_status(status)
	
	for i, f in enumerate(fields):
            # It might be better to do this check in the Python code
            if f not in var_labels:
                if f in static_fields :
			


        free(primary_variables)
        free(fpoint)
        free(forder)

        print 'done filling particle variables', count_all_particles



    #@cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def grid_var_fill(self, source, fields):
        cdef int num_oct_levels
        cdef int *num_octs_per_level
        cdef float *variables
        cdef int status
        cdef int root_oct, child, order
        cdef int ix, iy, iz
        cdef int cx, cy, cz
        cdef double dpos[3]
        cdef np.ndarray[np.float32_t, ndim=1] arr
        # This relies on the fields being contiguous
        cdef np.float32_t **fpoint
        cdef int nf = len(fields)
        fpoint = <np.float32_t**>malloc(sizeof(np.float32_t*)*nf)
        forder = <int*>malloc(sizeof(int)*nf)
        cdef int i, j, level

        # translate fields from ARTIO names to indices
        var_labels = self.parameters['grid_variable_labels']
        for i, f in enumerate(fields):
            # It might be better to do this check in the Python code
            if f not in var_labels:
                print "This field is not known to ARTIO:", f
                raise RuntimeError
            j = var_labels.index(f)
            arr = source[f]
            fpoint[i] = <np.float32_t *>arr.data
            forder[i] = j

        status = artio_grid_cache_sfc_range( self.handle, self.sfc_min, self.sfc_max )
        check_artio_status(status) 

        num_octs_per_level = <int *>malloc(self.max_level*sizeof(int))
        variables = <float *>malloc(8*self.num_grid_variables*sizeof(float))

        count = self.num_root_cells

        for sfc in range( self.sfc_min, self.sfc_max+1 ) :
            status = artio_grid_read_root_cell_begin( self.handle, sfc, 
                    dpos, variables, &num_oct_levels, num_octs_per_level )
            check_artio_status(status) 

            ix = <int>(dpos[0]-0.5) / 2
            iy = <int>(dpos[1]-0.5) / 2
            iz = <int>(dpos[2]-0.5) / 2

            cx = 0 if dpos[0] < (2*ix + 1) else 1
            cy = 0 if dpos[1] < (2*iy + 1) else 1
            cz = 0 if dpos[2] < (2*iz + 1) else 1
            
            root_oct = ix+(self.num_grid/2)*(iy+(self.num_grid/2)*iz)
            child = cx+2*(cy+2*cz)
            order = 8*root_oct + child

            assert( root_oct < self.num_root_cells / 8 )
            assert( child >= 0 and child < 8 )
            assert( order >= 0 and order < self.num_root_cells )

            for i in range(nf):
                fpoint[i][order] = variables[forder[i]]
 
            for level in range(1,num_oct_levels) :
                status = artio_grid_read_level_begin( self.handle, level )
                check_artio_status(status) 

                for oct in range(num_octs_per_level[level-1]) :
                    status = artio_grid_read_oct( self.handle, NULL, variables, NULL )
                    check_artio_status(status) 

                    for child in range(8) :
                        for i in range(nf):
                            fpoint[i][count] = variables[self.num_grid_variables*child+forder[i]]
                        count += 1
 
                status = artio_grid_read_level_end( self.handle )
                check_artio_status(status) 

            status = artio_grid_read_root_cell_end( self.handle )
            check_artio_status(status) 
        
        status = artio_grid_clear_sfc_cache( self.handle )
        check_artio_status(status)

        free(num_octs_per_level) 
        free(variables)
        free(fpoint)
        free(forder)

        print 'done filling oct variables', count

###################################################
def artio_is_valid( char *file_prefix ) :
    cdef artio_fileset_handle *handle = artio_fileset_open( file_prefix, 
            ARTIO_OPEN_HEADER, artio_context_global )
    if handle == NULL :
        return False
    else :
        artio_fileset_close(handle) 
    return True
