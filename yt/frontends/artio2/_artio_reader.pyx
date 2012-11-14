# warning: this may not be portable
from cpython.cobject cimport PyCObject_FromVoidPtr, PyCObject_AsVoidPtr

from libc.stdint cimport int32_t, int64_t

cdef extern from "artio.h" :
    ctypedef struct artio_file_struct "artio_file_struct" :
        pass
    ctypedef artio_file_struct *artio_file "artio_file"

    ctypedef struct artio_context_struct "artio_context_struct" :
        pass
    ctypedef artio_context_struct *artio_context "artio_context"
    cdef artio_context artio_context_global

    cdef int ARTIO_OPEN_HEADER "ARTIO_OPEN_HEADER"
    cdef int ARTIO_OPEN_GRID "ARTIO_OPEN_GRID"
    cdef int ARTIO_OPEN_PARTICLES "ARTIO_OPEN_PARTICLES" 

    cdef int ARTIO_SUCCESS "ARTIO_SUCCESS"

    artio_file artio_fileset_open(char *file_prefix, int type, artio_context context )
    int artio_fileset_close( artio_file handle )
    int artio_fileset_open_grid( artio_file )
    int artio_fileset_open_particle( artio_file )

    # parameter functions
    int artio_parameter_iterate( artio_file handle, char *key, int *type, int *length )
    int artio_parameter_get_array_length(artio_file handle, char * key, int *length)

    int artio_parameter_get_int(artio_file handle, char * key, int32_t * value)
    int artio_parameter_get_int_array(artio_file handle, char * key, int length,
            int32_t *values)

class artio_fileset :
    def __init__(self, char *file_prefix) :
        self.handle = PyCObject_FromVoidPtr(artio_fileset_open( file_prefix,
                ARTIO_OPEN_HEADER, artio_context_global ),NULL)

def artio_is_valid( char *file_prefix ) :
    cdef artio_file local_handle = artio_fileset_open( file_prefix, 
            ARTIO_OPEN_HEADER, artio_context_global )
    if local_handle == NULL :
        return False;
    else :
        artio_fileset_close(local_handle) 
    return True
