"""

"""

from libc.stdint cimport int32_t, int64_t
from libc.stdlib cimport malloc, free

cdef struct artio_context_struct: 
    int comm
ctypedef artio_context_struct * artio_context
cdef artio_context artio_context_global


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

    # errors
    cdef int ARTIO_SUCCESS "ARTIO_SUCCESS"

    artio_file artio_fileset_open(char *file_prefix, int type, artio_context context )
    int artio_fileset_close( artio_file handle )
    int artio_fileset_open_grid( artio_file )
    int artio_fileset_open_particle( artio_file )

    # parameter functions
    int artio_parameter_iterate( artio_file handle, char *key, int *type, int *length )
    int artio_parameter_get_int_array(artio_file handle, char * key, int length, int32_t *values)
    int artio_parameter_get_float_array(artio_file handle, char * key, int length, float *values)
    int artio_parameter_get_long_array(artio_file handle, char * key, int length, int64_t *values)
    int artio_parameter_get_double_array(artio_file handle, char * key, int length, double *values)
    int artio_parameter_get_string_array(artio_file handle, char * key, int length, char **values, int max_length)

########## wrappers calling c
cdef extern from "artio.c":
    artio_file artio_fileset_open(char * file_prefix, int type, artio_context context) 

cpdef read_header(char * file_prefix, int type, rheader) :
    cdef artio_file junk
    cdef artio_context context=NULL
    print "file_prefix=",file_prefix,"type=",type
    junk = artio_fileset_open(file_prefix, type, context)


########## wrappers called by python
cdef class artio_fileset :
    cdef artio_file handle
    cdef public object parameters
 
    def __init__(self, char *file_prefix) :
        self.handle = artio_fileset_open( file_prefix, ARTIO_OPEN_HEADER, artio_context_global )
        self.read_parameters()

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

        print self.parameters
            

def artio_is_valid( char *file_prefix ) :
    cdef artio_file handle = artio_fileset_open( file_prefix, 
            ARTIO_OPEN_HEADER, artio_context_global )
    if handle == NULL :
        return False;
    else :
        artio_fileset_close(handle) 
    return True
