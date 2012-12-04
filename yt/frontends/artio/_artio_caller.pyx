"""

"""

from libc.stdint cimport int32_t, int64_t
from libc.stdlib cimport malloc, free
#from  data_structures  import cell_pos_callback
import  data_structures  

cdef struct artio_context_struct: 
    int comm
ctypedef artio_context_struct * artio_context
cdef artio_context artio_context_global


cdef extern from "sfc.h":
    ctypedef void sfc_coords( int index, int coords[3], int )

cdef extern from "artio.h":
    ctypedef struct artio_file_struct "artio_file_struct" :
        pass
    ctypedef artio_file_struct *artio_file "artio_file"

    ctypedef struct artio_context_struct "artio_context_struct" :
        pass
    ctypedef artio_context_struct *artio_context "artio_context"


    cdef artio_context artio_context_global
    

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
    

    # errors
    cdef int ARTIO_SUCCESS "ARTIO_SUCCESS"
    cdef int ARTIO_ERR_MEMORY_ALLOCATION "ARTIO_ERR_MEMORY_ALLOCATION"

    artio_file artio_fileset_open(char *file_prefix, int type, artio_context context )
    int artio_fileset_close( artio_file handle )
    int artio_fileset_open_particle( artio_file handle )
    int artio_fileset_open_grid(artio_file handle) 
    int artio_fileset_close_grid(artio_file handle) 

    # parameter functions
    int artio_parameter_iterate( artio_file handle, char *key, int *type, int *length )
    int artio_parameter_get_int_array(artio_file handle, char * key, int length, int32_t *values)
    int artio_parameter_get_float_array(artio_file handle, char * key, int length, float *values)
    int artio_parameter_get_long_array(artio_file handle, char * key, int length, int64_t *values)
    int artio_parameter_get_double_array(artio_file handle, char * key, int length, double *values)
    int artio_parameter_get_string_array(artio_file handle, char * key, int length, char **values, int max_length)
    int artio_grid_read_root_cell_end(artio_file handle)
    int artio_grid_read_root_nocts(artio_file handle, int64_t sfc,\
                                  float *variables, int32_t *num_oct_levels, int32_t *num_octs_per_level)

    ctypedef void (* GridCallBackPos)(float * variables, int level, int refined, int64_t sfc_index, double pos[3])
    int artio_grid_read_sfc_range_pos(artio_file handle,\
                int64_t sfc1, int64_t sfc2,\
                int min_level_to_read, int max_level_to_read,\
                int options,\
                GridCallBackPos callback)
    ctypedef void (* GridCallBack)(float * variables, int level, int refined,int64_t sfc_index)
    int artio_grid_read_sfc_range(artio_file handle,\
                int64_t sfc1, int64_t sfc2,\
                int min_level_to_read, int max_level_to_read,\
                int options,\
                GridCallBack callback)
    
########## wrappers calling c
cdef extern from "artio.c":
    artio_file artio_fileset_open(char * file_prefix, int type, artio_context context) 


cdef class read_parameters_artio : 
    cdef public object parameters 
    cdef artio_file handle

    def __init__(self, char *file_prefix, int artio_type) :
        self.handle = artio_fileset_open( file_prefix, artio_type, artio_context_global ) 
        self.read_parameters()
        artio_fileset_close(self.handle)  #snl why didn't Doug close?

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

#        print self.parameters
        
#    def read_parameters(self) : 
#        cdef char key[64]
#        cdef int type
#        cdef int length
#        while artio_parameter_iterate( self.handle, key, &type, &length ) == ARTIO_SUCCESS :
#            print 'hi!!'

cdef class artio_fileset :
    cdef public object parameters 
    cdef artio_file handle
    def __init__(self, char *file_prefix) :
        cdef int artio_type = ARTIO_OPEN_HEADER
        self.handle = artio_fileset_open( file_prefix, artio_type, artio_context_global ) 
        d = read_parameters_artio(file_prefix, artio_type)
        self.parameters = {}
        self.parameters = d.parameters
        print self.parameters
        print 'done reading header parameters'

cdef class artio_fileset_grid :
    cdef public object parameters #is self.parameters public or is parameters public?
    cdef artio_file handle

    def __init__(self, char *file_prefix) :
        cdef int artio_type = ARTIO_OPEN_GRID
        self.handle = artio_fileset_open( file_prefix, artio_type, artio_context_global ) 
        artio_fileset_open_grid( self.handle ) 
        d = read_parameters_artio(file_prefix, artio_type)
        self.parameters = {}
        self.parameters = d.parameters
#        print self.parameters
        print 'done reading grid parameters'

###### callback for positions #############
def count_octs(char *file_prefix,\
                   int64_t sfc1, int64_t sfc2,\
                   int min_level_to_read, int max_level_to_read,\
                   int num_grid_variables
               ) :
    #max_level_to_read is currently unused, but could be useful
    cdef artio_file handle
    handle = artio_fileset_open( file_prefix, ARTIO_OPEN_GRID, artio_context_global ) 
    artio_fileset_open_grid( handle ) 
    num_total_octs = 0 
    cdef float * variables  
    cdef int32_t * num_oct_levels 
    cdef int32_t * num_octs_per_level 
    length = num_grid_variables * 8
    variables = <float *>malloc(length*sizeof(float))
    num_oct_levels = <int32_t *>malloc(1*sizeof(int32_t))
    length = max_level_to_read
    num_octs_per_level = <int32_t *>malloc(length*sizeof(int32_t))
    
    for sfc in xrange(sfc1,sfc2):
        artio_grid_read_root_nocts(handle,\
                                  sfc,\
                                  variables, num_oct_levels,\
                                  num_octs_per_level)
        count_level_octs = {}          
        count_level_octs = [ num_octs_per_level[i] for i in xrange(min(max_level_to_read,1)) ]
        num_sfc_octs = sum(count_level_octs)
        num_total_octs += num_sfc_octs

    artio_grid_read_root_cell_end(handle)
    return num_total_octs

###### callback for positions #############
def grid_pos_fill(char *file_prefix,\
                int64_t sfc1, int64_t sfc2,\
                int min_level_to_read, int max_level_to_read) :
    cdef artio_file handle
    handle = artio_fileset_open( file_prefix, ARTIO_OPEN_GRID, artio_context_global ) 
    artio_fileset_open_grid( handle ) 
    status = artio_grid_read_sfc_range_pos(handle,\
                sfc1, sfc2,\
                min_level_to_read, max_level_to_read,\
                ARTIO_READ_LEAFS,\
                wrap_cell_pos_callback)
    if status != ARTIO_SUCCESS :
        print "exiting sfc range read with error status", status
        if status == ARTIO_ERR_MEMORY_ALLOCATION :
            print "cannot allocate enough memory in one sfc range,",\
                "try loading smaller pieces"

cdef void wrap_cell_pos_callback(float *variables, int level, int refined, int64_t sfc_index, double *pos):
    position = {}
    position = [ pos[i] for i in range(3) ]
    data_structures.cell_pos_callback(level, refined, sfc_index, position)

###### callback for variables #############
def grid_var_fill(char *file_prefix,\
                int64_t sfc1, int64_t sfc2,\
                int min_level_to_read, int max_level_to_read) :
    cdef artio_file handle
    handle = artio_fileset_open( file_prefix, ARTIO_OPEN_GRID, artio_context_global ) 
    artio_fileset_open_grid( handle ) 
    artio_grid_read_sfc_range(handle,\
                sfc1, sfc2,\
                min_level_to_read, max_level_to_read,\
                ARTIO_READ_LEAFS,\
                wrap_cell_var_callback)
def cell_var_callback(level, refined, sfc_index, cell_var):
    print "variable callback success! ",level, refined, sfc_index 
cdef void wrap_cell_var_callback(float *variables, int level, int refined, int64_t sfc_index):
    cell_var={}
#     cdef int num_grid_variables=1 # really from read_parameters_artio(file_prefix, artio_type)
#        for label in parameters.grid_variable_labels
#            cell_var = [ variables[i] for i in range(num_grid_variables) ]
#            if(label == 'density'): 
#                cell_var_callback(cell_var, level, refined, sfc_index)
    cell_var_callback(level, refined, sfc_index, cell_var)
 
def artio_is_valid( char *file_prefix ) :
    cdef artio_file handle = artio_fileset_open( file_prefix, 
            ARTIO_OPEN_HEADER, artio_context_global )
    if handle == NULL :
        return False;
    else :
        artio_fileset_close(handle) 
    return True
