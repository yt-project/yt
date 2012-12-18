"""

"""
import numpy as np
cimport numpy as np
import sys 

from libc.stdint cimport int32_t, int64_t
from libc.stdlib cimport malloc, free
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
    ctypedef struct artio_grid_file_struct "artio_grid_file_struct" :
        pass
    ctypedef artio_grid_file_struct *artio_grid_file "artio_grid_file"


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
    int artio_grid_cache_sfc_range(artio_file handle, int64_t start, int64_t end)

    ctypedef void (* GridCallBackPos)(float * variables, int level, int refined, int64_t sfc_index, double pos[3], void *)
    int artio_grid_read_sfc_range_pos(artio_file handle,\
                int64_t sfc1, int64_t sfc2,\
                int min_level_to_read, int max_level_to_read,\
                int options,\
                GridCallBackPos callback, void *user_data)

    ctypedef void (* GridCallBackBuffer)(float * variables, int level, int refined, int64_t sfc_index, void *)
    int artio_grid_read_sfc_range_buffer(artio_file handle,\
                int64_t sfc1, int64_t sfc2,\
                int min_level_to_read, int max_level_to_read,\
                int options,\
                GridCallBackBuffer callback, void *user_data)

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

def check_artio_status(status, fname):
    callername = sys._getframe().f_code.co_name
    nline = sys._getframe().f_lineno
    if status!=ARTIO_SUCCESS :
        print 'failure with status', status, 'in function',fname,'from caller', callername, nline 
        sys.exit(1)
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
    cdef public object parameters
    cdef artio_file handle
    cdef public file_prefix

    def __init__(self, char *file_prefix) :
        cdef int artio_type = ARTIO_OPEN_GRID
        self.handle = artio_fileset_open( file_prefix, artio_type, artio_context_global ) 
        artio_fileset_open_grid( self.handle ) 
        d = read_parameters_artio(file_prefix, artio_type)
        self.parameters = {}
        self.parameters = d.parameters
        self.file_prefix = file_prefix
        print 'done reading grid parameters'

#snl: subclass some of this?
class artio_grid_routines(object) : 
    def __init__(self, param_handle) :
        self.oct_handler=None
        self.source = None
        self.level_count = None

        self.min_level_to_read = 0
        self.max_level_to_read = param_handle.parameters['grid_max_level'][0]
        self.sfc1 = 0
        self.sfc2 = param_handle.parameters['grid_file_sfc_index'][1]-1
        self.num_grid_variables = param_handle.parameters['num_grid_variables'][0]
        self.param_handle=param_handle
        self.ng=1 
        self.cpu=0
        self.domain_id=0

        self.grid_variable_labels=param_handle.parameters['grid_variable_labels']
 
        # dictionary from artio to yt ... should go elsewhere
        self.label_artio_yt = {}
        for label in self.grid_variable_labels :
            if label == 'HVAR_GAS_DENSITY' :
                self.label_artio_yt[label] = 'Density'
            else :
                self.label_artio_yt[label] = label
        print self.label_artio_yt


        self.grid_variable_labels        
        self.label_index = {}
        self.matched_fields = []
        # not sure about file handling
        # self.handle = <object> handle #<------seg faults
        # and you never bother to close all of these handles
        
    def count_octs(self) :
        cdef int min_level_to_read = self.min_level_to_read
        cdef int max_level_to_read = self.max_level_to_read
        cdef int64_t sfc1 = self.sfc1
        cdef int64_t sfc2 = self.sfc2
        cdef num_grid_variables = self.num_grid_variables
        
        cdef float * variables  
        cdef int32_t * num_oct_levels 
        cdef int32_t * num_octs_per_level 
        
        length = num_grid_variables * 8
        variables = <float *>malloc(length*sizeof(float))
        num_oct_levels = <int32_t *>malloc(1*sizeof(int32_t))
        length = max_level_to_read
        num_octs_per_level = <int32_t *>malloc(length*sizeof(int32_t))
        
        cdef artio_file handle
        handle = artio_fileset_open( self.param_handle.file_prefix, 
                                     ARTIO_OPEN_GRID, artio_context_global ) 
        artio_fileset_open_grid( handle ) 
        
        cdef int64_t num_total_octs =0
        n_levels = max_level_to_read - min_level_to_read + 1
        level_count = np.zeros(n_levels, dtype='int64')

        status = artio_grid_cache_sfc_range(handle, sfc1, sfc2)
        check_artio_status(status, artio_grid_routines.__name__)
        for sfc in xrange(sfc1,sfc2):
            status = artio_grid_read_root_nocts(handle, sfc,
                                                variables, num_oct_levels,
                                                num_octs_per_level)
            check_artio_status(status, artio_grid_routines.__name__)
            noct_levels = num_oct_levels[0]
            count_level_octs = {}          
            count_level_octs = [ num_octs_per_level[i] for i in xrange(noct_levels) ]
            for level in xrange(noct_levels) : 
                level_count[level] += count_level_octs[level] 
            num_sfc_octs = sum(count_level_octs)
            num_total_octs += num_sfc_octs
            status = artio_grid_read_root_cell_end(handle)
            check_artio_status(status, artio_grid_routines.__name__)
            
            #    check_artio_status(-1, count_octs.__name__)
            # dont close file until the end of the object... add __del__: artio_fileset_close_grid(handle)
        self.level_count = level_count
        return num_total_octs

    def grid_pos_fill(self, oct_handler) :
        self.oct_handler = oct_handler
        cdef artio_file handle
        if self.oct_handler == None :
            print 'oct_handler is not assigned!'
            sys.exit(1)
        handle = artio_fileset_open( self.param_handle.file_prefix, 
                                     ARTIO_OPEN_GRID, artio_context_global ) 
        status = artio_grid_read_sfc_range_pos(handle,\
                    self.sfc1, self.sfc2,\
                    self.min_level_to_read, self.max_level_to_read,\
                    ARTIO_READ_REFINED,\
                    wrap_cell_pos_callback, <void*>self) 
        check_artio_status(status, artio_grid_routines.__name__)
        print 'done filling oct positions'
    def cell_pos_callback(self, level, refined, sfc_index, pos):
        self.oct_handler.add(self.cpu + 1, level - self.min_level_to_read, self.ng, pos, self.domain_id)

    def grid_var_fill(self, source, fields):
        self.source = source
        i=-1
        for artlabel in self.grid_variable_labels :
            label = self.label_artio_yt[artlabel]
            i=i+1
            for field in fields : 
                if field == label :
                    print 'match, in fields?', field,label, fields
                    print '!!!!!!!!!!!!!!!'
                    self.label_index[field]=i
                    self.matched_fields.append(field)
        print 'matched fields:',self.matched_fields
        print 'art index of matched fields',self.label_index
        self.count=0
        cdef artio_file handle
        if len(self.label_index) > 0 :
            handle = artio_fileset_open( self.param_handle.file_prefix, 
                                         ARTIO_OPEN_GRID, artio_context_global ) 
            status = artio_grid_read_sfc_range_buffer(
                handle, self.sfc1, self.sfc2,\
                    self.min_level_to_read, self.max_level_to_read,\
                    ARTIO_READ_REFINED,\
                    wrap_cell_var_callback,\
                    <void*>self
                ) #only octs!
            check_artio_status(status, artio_grid_routines.__name__)
        print 'done buffering variables'
    def cell_var_callback(self, level, refined, ichild, cell_var):
        for field in self.matched_fields : 
            self.source[field][self.count] = cell_var[self.label_index[field]] 
        self.count=self.count+1
 

###### callbacks #############
cdef void wrap_cell_pos_callback(float *variables, int level, int refined, 
                                 int64_t sfc_index, double *pos, void *user_data):
    position = np.empty((1, 3), dtype='float64')
    position[0,0] = pos[0] 
    position[0,1] = pos[1] 
    position[0,2] = pos[2] 
    artioroutines = <object>user_data
    artioroutines.cell_pos_callback(level, refined, sfc_index, position)

cdef void wrap_cell_var_callback(float *variables, int level, int refined, 
                                 int64_t ichild, void *user_data):
    artioroutines = <object>user_data
    cell_var={}
    cell_var = [ variables[i] for i in range(artioroutines.num_grid_variables) ]
    artioroutines.cell_var_callback(level, refined, ichild, cell_var)


def artio_is_valid( char *file_prefix ) :
    cdef artio_file handle = artio_fileset_open( file_prefix, 
            ARTIO_OPEN_HEADER, artio_context_global )
    if handle == NULL :
        return False;
    else :
        artio_fileset_close(handle) 
    return True



