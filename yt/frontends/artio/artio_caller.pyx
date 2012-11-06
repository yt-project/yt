#struct artio_context_struct {
#        MPI_Comm comm;
#};
cdef struct artio_context_struct: #in artio_mpi.h MPI_Comm
    int comm

#typedef struct artio_context_struct * artio_context;
ctypedef artio_context_struct * artio_context

cdef extern from "artio_headers/artio.c":

     ctypedef int int64_t 

     ctypedef struct param: 
         int key_length
         char key[64]
         int val_length
         int type 
         char *value
         
     ctypedef struct list: 
         param * head
         param * tail
         param * cursor
         int iterate_flag 
         int endian_swap 

     #    artio_internal.h
     ctypedef struct artio_particle_file:
        #	artio_fh *ffh
        int num_particle_files
        int64_t *file_sfc_index
        int64_t cache_sfc_begin
        int64_t cache_sfc_end
        int64_t *sfc_offset_table
        
        #/* maintained for consistency and user-error detection */
        int num_species
        int cur_file
        int cur_species
        int cur_particle
        int64_t cur_sfc
        int *num_primary_variables
        int *num_secondary_variables
        int *num_particles_per_species
        
     ctypedef struct artio_grid_file:
        #        artio_fh *ffh
        int num_grid_variables
        int num_grid_files
        int64_t *file_sfc_index
        int64_t cache_sfc_begin
        int64_t cache_sfc_end
        int64_t *sfc_offset_table
        
        int file_max_level
        #/* maintained for consistency and user-error detection */
        int cur_file
        int cur_num_levels
        int cur_level
        int cur_octs
        int64_t cur_sfc
        int *octs_per_level
        
        
     ctypedef struct artio_file:
        char file_prefix[256]
        int endian_swap
        int open_type
        int open_mode
        int rank
        int num_procs
        artio_context context #artio_mpi.h MPI_Comm
        
        int64_t *proc_sfc_index
        int64_t proc_sfc_begin
        int64_t proc_sfc_end
        int64_t num_root_cells
        
        list param_list
        artio_grid_file grid
        artio_particle_file particle
        
        
     artio_file artio_fileset_open(char * file_prefix, int type, artio_context context) 

     void wrap_artio_fileset_open(char *file_prefix, int type)
   

#python only passes numbers and strings... no pointers or structures
cpdef read_header(char * file_prefix, int type): 
     wrap_artio_fileset_open(file_prefix, type)
        

