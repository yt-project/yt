"""

"""
import numpy as np
cimport numpy as np
import sys 

from libc.stdint cimport int32_t, int64_t
from libc.stdlib cimport malloc, free
import  data_structures  

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

def check_artio_status(status, fname="[unknown]"):
    if status!=ARTIO_SUCCESS :
        callername = sys._getframe().f_code.co_name
        nline = sys._getframe().f_lineno
        print 'failure with status', status, 'in function',fname,'from caller', callername, nline 
        sys.exit(1)

cdef class artio_fileset :
    cdef public object parameters 
    cdef object oct_handler
    cdef artio_fileset_handle *handle
    cdef int64_t num_root_cells
    cdef int64_t sfc_min, sfc_max
    cdef public int num_grid

    # grid attributes
    cdef int min_level, max_level
    cdef int num_grid_variables

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

    def grid_pos_fill(self, oct_handler) :
        ''' adds all refined octs and a new array of ghost octs for  
        the "fully refined octs" at level=-1 in ART or 0 in yt convention 
        so that root level can consist of children
        '''
        cdef int64_t sfc
        cdef int level
        cdef int num_oct_levels
        cdef int *num_octs_per_level
        cdef double dpos[3]
        cdef int oct_count

        pos = np.empty((1,3), dtype='float64')

        print 'start filling oct positions'
        self.oct_handler = oct_handler

        oct_count = 0
        level = 0
        for iz in range(self.num_grid/2) :
            pos[0,2] = iz*2+1
            for iy in range(self.num_grid/2) :
                pos[0,1] = iy*2+1
                for ix in range(self.num_grid/2) :
                    pos[0,0]=ix*2+1
                    self.oct_handler.add(self.cpu+1, level, self.num_grid, pos, self.domain_id) 
                    oct_count += 1

        # Now do real ART octs
        print 'start filling oct positions children'
        status = artio_grid_cache_sfc_range( self.handle, self.sfc_min, self.sfc_max )
        check_artio_status(status) 

        num_octs_per_level = <int *>malloc(self.max_level*sizeof(int))

        for sfc in range( self.sfc_min, self.sfc_max+1 ) :
            status = artio_grid_read_root_cell_begin( self.handle, sfc, 
                dpos, NULL, &num_oct_levels, num_octs_per_level )
            check_artio_status(status) 

            for level in range(num_oct_levels) :
                status = artio_grid_read_level_begin( self.handle, level+1 )
                check_artio_status(status) 

                for oct in range(num_octs_per_level[level]) :
                    status = artio_grid_read_oct( self.handle, dpos, NULL, NULL )
                    check_artio_status(status) 
 
                    pos[0,0] = dpos[0]
                    pos[0,1] = dpos[1]
                    pos[0,2] = dpos[2]
                    oct_count += self.oct_handler.add(self.cpu+1, level+1, self.num_grid, pos, self.domain_id )

                status = artio_grid_read_level_end( self.handle )
                check_artio_status(status) 

            status = artio_grid_read_root_cell_end( self.handle )
            check_artio_status(status) 
        
        status = artio_grid_clear_sfc_cache( self.handle )
        check_artio_status(status)

        free(num_octs_per_level) 

        print 'done filling oct positions', oct_count

    def grid_var_fill(self, source, fields):
        cdef int num_oct_levels
        cdef int *num_octs_per_level
        cdef float *variables
        cdef int status
        cdef int root_oct, child, order
        cdef int ix, iy, iz
        cdef int cx, cy, cz
        cdef int count
        cdef double dpos[3]

        print "Field list:", fields
 
        # translate fields from ARTIO names to indices
        if not all(f in self.parameters['grid_variable_labels'] for f in fields) :
            print "Asked for a variable that doesn't exist!"
            sys.exit(1)
        field_order = dict([(f,self.parameters['grid_variable_labels'].index(f)) for f in fields])

        print "field order: ", field_order
 
        status = artio_grid_cache_sfc_range( self.handle, self.sfc_min, self.sfc_max )
        check_artio_status(status) 

        num_octs_per_level = <int *>malloc(self.max_level*sizeof(int))
        variables = <float *>malloc(8*self.num_grid_variables*sizeof(float))

        #art_order = {}
        #for ix in range(2) :
        #    for iy in range(2) :
        #        for iz in range(2) :
        #        #    art_order[ix+2*(iy+2*iz)] = iz+2*(iy+2*ix)
        #            art_order[iz+2*(iy+2*ix)] = ix+2*(iy+2*iz)
        #
        #for i in range(8) :
        #    print i, art_order[i]

        count = self.num_root_cells
        seen = [False for i in range(self.num_root_cells)]

        for sfc in range( self.sfc_min, self.sfc_max+1 ) :
            status = artio_grid_read_root_cell_begin( self.handle, sfc, 
                    dpos, variables, &num_oct_levels, num_octs_per_level )
            check_artio_status(status) 

            ix = (int)(dpos[0]-0.5) / 2
            iy = (int)(dpos[1]-0.5) / 2
            iz = (int)(dpos[2]-0.5) / 2

            cx = 0 if dpos[0] < (2*ix + 1) else 1
            cy = 0 if dpos[1] < (2*iy + 1) else 1
            cz = 0 if dpos[2] < (2*iz + 1) else 1
            
            root_oct = iz+(self.num_grid/2)*(iy+(self.num_grid/2)*ix)
            child = cz+2*(cy+2*cx)
            order = 8*root_oct + child

            assert( root_oct < self.num_root_cells / 8 )
            assert( child >= 0 and child < 8 )
            assert( order >= 0 and order < self.num_root_cells )

            assert( not seen[order] )
            seen[order] = True

            for f in fields :
                source[f][order] = variables[field_order[f]]
 
            for level in range(num_oct_levels) :
                status = artio_grid_read_level_begin( self.handle, level+1 )
                check_artio_status(status) 

                for oct in range(num_octs_per_level[level]) :
                    status = artio_grid_read_oct( self.handle, NULL, variables, NULL )
                    check_artio_status(status) 

                    for child in range(8) :
                        for f in fields :
                            source[f][count] = variables[self.num_grid_variables*child+field_order[f]]
                        count += 1
 
                status = artio_grid_read_level_end( self.handle )
                check_artio_status(status) 

            status = artio_grid_read_root_cell_end( self.handle )
            check_artio_status(status) 
        
        status = artio_grid_clear_sfc_cache( self.handle )
        check_artio_status(status)

        free(num_octs_per_level) 
        free(variables)

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
