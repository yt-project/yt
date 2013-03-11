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

cdef extern from "stdlib.h":
    void *alloca(int)

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
        print 'failure with status', status, 'in function',fname,'from caller', callername, nline 
        sys.exit(1)


cdef class artio_fileset :
    cdef public object parameters 
    cdef artio_fileset_handle *handle

    # grid attributes
    cdef public int num_grid
    cdef int64_t num_root_cells
    cdef int min_level, max_level
    cdef int64_t sfc_file_min, sfc_file_max
    cdef int num_grid_variables, num_species
    
    cdef public object fnART_primary, fnART_secondary, fnART_static
    cdef public object duplicate_species, labels_species 

    def __init__(self, char *file_prefix) :
        cdef int artio_type = ARTIO_OPEN_HEADER
        cdef int64_t num_root

        self.handle = artio_fileset_open( file_prefix, artio_type, artio_context_global ) 
        self.read_parameters()
        #print 'print parameters in caller.pyx',self.parameters
        print 'done reading header parameters'

        self.num_root_cells = self.parameters['num_root_cells'][0]
        self.num_grid = 1
        num_root = self.num_root_cells
        while num_root > 1 :
            self.num_grid <<= 1
            num_root >>= 3

        self.min_level = 0
        self.max_level = self.parameters['grid_max_level'][0]
        self.sfc_file_min = 0
        self.sfc_file_max = self.parameters['grid_file_sfc_index'][1]-1
        self.num_grid_variables = self.parameters['num_grid_variables'][0]

        # ART fieldnames 
        self.num_species = self.parameters['num_particle_species'][0]
        self.labels_species = self.parameters['particle_species_labels']
        self.fnART_primary={}
        self.fnART_secondary={}
        self.fnART_static={}
        for ispec in range(self.num_species) : 
            listdict = "species_%02d_primary_variable_labels" % ispec
            self.fnART_primary[ispec] = self.parameters[listdict]
            if self.parameters["num_secondary_variables"][ispec] > 0 :
                listdict = "species_%02d_secondary_variable_labels" % ispec
                self.fnART_secondary[ispec] = self.parameters[listdict]
            else : 
                self.fnART_secondary[ispec] = []
                
            # N-BODY species have a static label called N-Body for MASS    
            if self.labels_species[ispec] == 'N-BODY' :
                self.fnART_static[ispec] = ["MASS"]
            else : 
                self.fnART_static[ispec] = [] 
            self.fnART_static[ispec].append("particle_index") 
        prev_specie = '0'
        self.duplicate_species = {}
        print self.labels_species  
        for specie in self.labels_species : 
            print 'hi',specie
            if prev_specie == specie : 
                self.duplicate_species[specie].append(specie)
            else : 
                self.duplicate_species[specie] = []
            prev_specie = specie

        #kln - add particle detection code
        status = artio_fileset_open_particles( self.handle )
        check_artio_status(status)
 
        # dhr - add grid detection code 
        status = artio_fileset_open_grid( self.handle )
        check_artio_status(status)

            
    def all_indices(value, qlist):
        indices = []
        idx = -1
        while True:
            try:
                idx = qlist.index(value, idx+1)
                indices.append(idx)
            except ValueError:
                break
        return indices

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
                print "ERROR: invalid type!"

            self.parameters[key] = parameter

#    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def particle_var_fill_mask(self, SelectorObject selector, int64_t sfc_start, int64_t sfc_end, accessed_species, fields) :
        # mask should be a function of accessed species 
        # (and indexed by a count that is species dependent) 
        # ispec is index out of all specs and aspecs is index out of accessed
        cdef double **primary_variables
        cdef float **secondary_variables
        cdef int **field_to_index
        cdef int *iacc_to_ispec 
        cdef int status
        cdef np.ndarray[np.float32_t, ndim=1] arr
        cdef int **mask
        cdef int *num_particles_per_species 
        cdef int **pos_index

        cdef int *subspecies
        subspecies = <int*>malloc(sizeof(int))
        cdef int64_t *pid
        pid = <int64_t *>malloc(sizeof(int64_t))

        cdef int num_fields = len(fields)
        cdef int i, j, level
        cdef np.float64_t pos[3]
        cdef np.float64_t dds[3]
        cdef int eterm[3]
        for i in range(3) : dds[i] = 0

        selected_particles = {}


        # accessed species includes duplicate naming (e.g. all nbody specs)
        for aspec, specie in enumerate(accessed_species):
            if self.duplicate_species[specie] : #not empty list is true
                accessed_species.insert(aspec+1,self.duplicate_species[specie])
        # field names + counts as a function of aspec
        num_acc_species = len(accessed_species)
        num_fieldnames = np.zeros(num_acc_species,dtype="int32")
        fieldnames={}
        for aspec, specie in enumerate(accessed_species):
            fieldnames[aspec]=[]
            for fieldtype, fieldname in fields:
                if specie == fieldtype : 
                    num_fieldnames[aspec] += 1
                    fieldnames[aspec].append(fieldname)

        # translate from yt specie index to ART specie index
        pos_index = <int**>malloc(sizeof(int*)*self.num_species)
        if not pos_index: raise MemoryError
        for ispec in range(self.num_species) : 
            pos_index[ispec] = <int*>malloc(3*sizeof(int))
            pos_index[ispec][0] = self.fnART_primary[ispec].index('POSITION_X')
            pos_index[ispec][1] = self.fnART_primary[ispec].index('POSITION_Y')
            pos_index[ispec][2] = self.fnART_primary[ispec].index('POSITION_Z')
        iacc_to_ispec = <int*>malloc(sizeof(int)*num_acc_species)
        if not iacc_to_ispec: raise MemoryError
        for aspec, specie in enumerate(accessed_species):
            ispec = self.labels_species.index(specie) #find first instance of species
            iacc_to_ispec[aspec] = ispec
            # duplicated species are neighbors
            if aspec > 0 and iacc_to_ispec[aspec] == iacc_to_ispec[aspec-1] : 
                iacc_to_ispec[aspec] = ispec+1
        # double check that iacc_to_ispec points to uniq indices
        for i in range(num_acc_species): 
            for j in range(i+1,num_acc_species):  
                if iacc_to_ispec[i]==iacc_to_ispec[j]:
                    print iacc_to_ispec[i]
                    print 'some accessed species indices point to the same ispec'
                    sys.exit(1)

        # allocate io pointers 
        num_particles_per_species =  <int *>malloc(sizeof(int)*self.num_species) 
        if not num_particles_per_species : raise MemoryError
        primary_variables = <double **>malloc(sizeof(double**)*num_acc_species)  
        secondary_variables = <float **>malloc(sizeof(float**)*num_acc_species)  
        if (not primary_variables) or (not secondary_variables) : raise MemoryError
        for aspec in range(num_acc_species) : 
            primary_variables[aspec]   = <double *>malloc(self.parameters['num_primary_variables'][aspec]*sizeof(double))
            secondary_variables[aspec] = <float *>malloc(self.parameters['num_secondary_variables'][aspec]*sizeof(float))
            if (not primary_variables[aspec]) or \
                    (not secondary_variables[aspec]) : raise MemoryError

        # cache the range
#        status = artio_particle_cache_sfc_range( self.handle, sfc_start, sfc_end )
        status = artio_grid_cache_sfc_range( self.handle, self.sfc_file_min, self.sfc_file_max )
        check_artio_status(status)

        # determine max number of particles we could hit (optimize later)
        max_particles= np.zeros(num_acc_species,dtype="int32")
        for sfc in range( self.sfc_start, self.sfc_end+1 ) :
            status = artio_particle_read_root_cell_begin( self.handle, sfc,
                    num_particles_per_species )
            check_artio_status(status)	
            for aspec in num_acc_species : 
                ispec = iacc_to_ispec[aspec]
                max_particles[aspec] += num_particles_per_species[ispec] 
            status = artio_particle_read_root_cell_end( self.handle )
            check_artio_status(status)

        # mask begin  ###################
        mask = <int**>malloc(sizeof(int*)*num_acc_species)
        if not mask : raise MemoryError
        for aspec in range(num_acc_species) :
            mask[aspec] = <int*>malloc( sizeof(int)*max_particles[aspec])
            if not mask[aspec]: raise MemoryError
        count_mask = []
        ipspec = []
        for aspec in range(num_acc_species) :
            count_mask.append(0)
            ipspec.append(0)
            ispec=iacc_to_ispec[aspec]

        print "generating mask for particles"
        for sfc in range( self.sfc_start, self.sfc_end+1 ) :
            status = artio_particle_read_root_cell_begin( 
                self.handle, sfc,
                num_particles_per_species )
            check_artio_status(status)

            for aspec in range(num_acc_species ) :
                ispec = iacc_to_ispec[aspec]
                status = artio_particle_read_species_begin(
                    self.handle, ispec)
                check_artio_status(status)

                for particle in range( num_particles_per_species[ispec] ) :
                    status = artio_particle_read_particle(
                        self.handle,
                        pid, subspecies, primary_variables[aspec],
                        secondary_variables[aspec])
                    check_artio_status(status)
                    pos[0] = primary_variables[aspec][pos_index[aspec][0]]
                    pos[1] = primary_variables[aspec][pos_index[aspec][1]]
                    pos[2] = primary_variables[aspec][pos_index[aspec][2]]
                    mask[aspec][ipspec[aspec]] = selector.select_cell(pos,dds,eterm)
                    count_mask[aspec] += mask[aspec][count_mask[aspec]]
                    ipspec[aspec] += 1
                status = artio_particle_read_species_end( self.handle )
                check_artio_status(status)
            status = artio_particle_read_root_cell_end( self.handle )
            check_artio_status(status)
        free(pos_index)
        print 'done masking'


	##########################################################
        # attribute fields to primary/secondary/static/empty
        # field_to_index, how_to_read, 
        field_to_index = <int**>malloc(sizeof(int*)*self.num_species)
        if not field_to_index: raise MemoryError
        how_to_read = {}
        for aspec in range(self.num_species) : 
            field_to_index[aspec] = <int*>malloc(num_fieldnames[aspec]*sizeof(int))
            if not field_to_index[aspec] : raise MemoryError
        countnbody = 0 
        for aspec in range(num_acc_species) : 
            ispec = iacc_to_ispec[aspec]
            for i, f in enumerate(fieldnames[aspec]):
                if   f in self.fnART_primary[ispec]:
                    how_to_read[ispec,i]= 'primary'
                    field_to_index[ispec][i] = self.fnART_primary[ispec].index(f)
                elif f in self.fnART_secondary[ispec]:
                    how_to_read[ispec,i]= 'secondary'
                    field_to_index[ispec][i] = self.fnART_secondary[ispec].index(f)
                elif f in self.fnART_static[ispec]:
                    # each new N-BODY spec adds one to the static mass location
                    if self.labels_species[ispec] == 'N-BODY' and f == 'MASS' :
                        how_to_read[ispec,i]= 'staticNBODY'
                        field_to_index[ispec][i] = countnbody
                        countnbody += 1 # MASS happens once per N-BODY species
                        print 'count the nbody species',countnbody
                    else :
                        how_to_read[ispec,i]= 'staticINDEX'
                else : 
                    how_to_read[ispec,i]= 'empty'
                    field_to_index[ispec][i] = 9999999
                print 'ispec', ispec,'field',f, 'how_to_read', how_to_read[ispec,i] 

        cdef np.float32_t **fpoint
        fpoint = <np.float32_t**>malloc(sizeof(np.float32_t*)*num_fields)
        if not fpoint : raise MemoryError
        
        for i, f in enumerate(fields):
            aspec = self.all_indices(accessed_species, f[0]) #duplicate indices
            num_selected_particles = 0
            for a in aspec :
                num_selected_particles += count_mask[a] # deal with duplicates

            selected_particles[f] = np.empty(num_selected_particles,dtype="float32")    
            arr = selected_particles[f]
            fpoint[i] = <np.float32_t *>arr.data

	##########################################################
        # now use mask to read fields
        print "reading in particle fields"
        for aspec in range(num_acc_species) :
             ipspec[aspec] = 0
        ip_all = 0
        for sfc in range( self.sfc_start, self.sfc_end+1 ) :
                status = artio_particle_read_root_cell_begin( self.handle, sfc,
                    num_particles_per_species )
                check_artio_status(status)	

                aoffset = np.zeros(num_acc_species,dtype="int32")
                prev_accessed_species  = None
                for aspec in range(num_acc_species) :
                    ispec = iacc_to_ispec[aspec]
                    if prev_accessed_species != accessed_species[aspec]: #duplicates go to the same fieldtype
                        aoffset[aspec] = num_fieldnames[aspec]+aoffset[aspec]
                    status = artio_particle_read_species_begin(self.handle, ispec);
                    check_artio_status(status)
                    
                    for particle in range( num_particles_per_species[ispec] ) :
                        
                        status = artio_particle_read_particle(self.handle,
                                        pid, subspecies, primary_variables[aspec],
                                        secondary_variables[aspec])
                        check_artio_status(status)
                        
                        if mask[aspec][ipspec[aspec]] == 1 :
                             for i in range(num_fieldnames[aspec]):
                                 j = i + aoffset[aspec]
                                 if not (how_to_read[ispec,i] == 'empty') : 
                                     assert(field_to_index[ispec][i]<100)
                                 if   how_to_read[ispec,i] == 'primary' : 
                                     fpoint[j][ip_all] = primary_variables[aspec][field_to_index[ispec][i]]
                                 elif how_to_read[ispec,i] == 'secondary' :
                                     fpoint[j][ip_all] = secondary_variables[aspec][field_to_index[ispec][i]]
                                 elif how_to_read[ispec,i] == 'staticNBODY' : 
                                     fpoint[j][ip_all] = self.parameters["particle_species_mass"][field_to_index[ispec][i]]
                                 elif how_to_read[ispec,i] == 'staticINDEX' : 
                                     fpoint[j][ip_all] = ip_all
                                 elif how_to_read[ispec,i] == 'empty' : 
                                     fpoint[j][ip_all] = 0
                                 else : 
                                     print 'undefined how to read in caller', how_to_read[ispec,i]
                                     print 'this should be impossible.'
                                     sys.exit(1)
                                 # print 'reading into fpoint', ip_all,fpoint[i][ip_all], fields[i]
                             ip_all += 1
                        ipspec[aspec] += 1
                    prev_accessed_species = accessed_species[aspec]
                    status = artio_particle_read_species_end( self.handle )
                    check_artio_status(status)
                    
                status = artio_particle_read_root_cell_end( self.handle )
                check_artio_status(status)
 

        free(subspecies)
        free(pid)
        free(num_particles_per_species)
        free(iacc_to_ispec)
        free(mask)
        free(field_to_index)
        free(primary_variables)
        free(secondary_variables)
        free(fpoint)
        print 'done filling particle variables', ip_all
        return selected_particles


#    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def particle_var_fill(self, SelectorObject selector, int64_t sfc_start, int64_t sfc_end, accessed_species, fields) :
        # ispec is index out of all specs and aspecs is index out of accessed
        # fields is art naming for yt fields, BUT if species is duplicated in ART there 
        #   is only one occurrence of the particle type in fields

        cdef double **primary_variables
        cdef float **secondary_variables
        cdef int **field_to_index
        cdef int *iacc_to_ispec 
        cdef int status
        cdef int *num_particles_per_species 
        cdef int **pos_index

        cdef int *subspecies
        subspecies = <int*>malloc(sizeof(int))
        cdef int64_t *pid
        pid = <int64_t *>malloc(sizeof(int64_t))

        cdef int num_fields = len(fields)
        cdef int i, j, level
        cdef np.float64_t pos[3]
        cdef np.float64_t dds[3]
        cdef int eterm[3]
        for i in range(3) : dds[i] = 0

        # accessed species includes duplicate naming (e.g. all nbody specs)
        for aspec, specie in enumerate(accessed_species):
            if self.duplicate_species[specie] : #not empty list is true
                accessed_species.insert(aspec+1,self.duplicate_species[specie])
        # field names + counts as a function of aspec
        num_acc_species = len(accessed_species)
        num_fieldnames = np.zeros(num_acc_species,dtype="int32")
        fieldnames={}
        for aspec, specie in enumerate(accessed_species):
            fieldnames[aspec]=[]
            for fieldtype, fieldname in fields:
                if specie == fieldtype : 
                    num_fieldnames[aspec] += 1
                    fieldnames[aspec].append(fieldname)
        data = [ np.empty(0,dtype="float32") for i in range(num_fields)]

        # translate from yt index to ART index
        pos_index = <int**>malloc(sizeof(int*)*self.num_species)
        if not pos_index: raise MemoryError
        for ispec in range(self.num_species) : 
            pos_index[ispec] = <int*>malloc(3*sizeof(int))
            pos_index[ispec][0] = self.fnART_primary[ispec].index('POSITION_X')
            pos_index[ispec][1] = self.fnART_primary[ispec].index('POSITION_Y')
            pos_index[ispec][2] = self.fnART_primary[ispec].index('POSITION_Z')
        iacc_to_ispec = <int*>malloc(sizeof(int)*num_acc_species)
        if not iacc_to_ispec: raise MemoryError
        for aspec, specie in enumerate(accessed_species):
            ispec = self.labels_species.index(specie) #find first instance of species
            iacc_to_ispec[aspec] = ispec
            if aspec > 0 and iacc_to_ispec[aspec] == iacc_to_ispec[aspec-1] :  # duplicated species are neighbors
                iacc_to_ispec[aspec] = ispec+1
        # double check that iacc_to_ispec points to uniq indices
        for i in range(num_acc_species): 
            for j in range(i+1,num_acc_species):  
                if iacc_to_ispec[i]==iacc_to_ispec[j]:
                    print iacc_to_ispec[i]
                    print 'some accessed species indices point to the same ispec'
                    sys.exit(1)

        # allocate io pointers 
        num_particles_per_species =  <int *>malloc(sizeof(int)*self.num_species) 
        if not num_particles_per_species : raise MemoryError
        primary_variables = <double **>malloc(sizeof(double**)*num_acc_species)  
        secondary_variables = <float **>malloc(sizeof(float**)*num_acc_species)  
        if (not primary_variables) or (not secondary_variables) : raise MemoryError
        for aspec in range(num_acc_species) : 
            primary_variables[aspec]   = <double *>malloc(self.parameters['num_primary_variables'][aspec]*sizeof(double))
            secondary_variables[aspec] = <float *>malloc(self.parameters['num_secondary_variables'][aspec]*sizeof(float))
            if (not primary_variables[aspec]) or \
                    (not secondary_variables[aspec]) : raise MemoryError


        # cache the range
        status = artio_particle_cache_sfc_range( self.handle, sfc_start, sfc_end )
        check_artio_status(status)

	
	##########################################################
        # attach fieldnames to primary/secondary/static/empty
        # assign: field_to_index, how_to_read, 
        field_to_index = <int**>malloc(sizeof(int*)*self.num_species)
        if not field_to_index: raise MemoryError
        how_to_read = {}
        for aspec in range(self.num_species) : 
            field_to_index[aspec] = <int*>malloc(num_fieldnames[aspec]*sizeof(int))
            if not field_to_index[aspec] : raise MemoryError
        countnbody = 0 
        for aspec in range(num_acc_species) : 
            ispec = iacc_to_ispec[aspec]
            for i, f in enumerate(fieldnames[aspec]):
                if   f in self.fnART_primary[ispec]:
                    how_to_read[ispec,i]= 'primary'
                    field_to_index[ispec][i] = self.fnART_primary[ispec].index(f)
                elif f in self.fnART_secondary[ispec]:
                    how_to_read[ispec,i]= 'secondary'
                    field_to_index[ispec][i] = self.fnART_secondary[ispec].index(f)
                elif f in self.fnART_static[ispec]:
                    # each new N-BODY spec adds one to the static mass location
                    if self.labels_species[ispec] == 'N-BODY' and f == 'MASS' :
                        how_to_read[ispec,i]= 'staticNBODY'
                        field_to_index[ispec][i] = countnbody 
                        countnbody += 1 # MASS happens once per N-BODY species
                        print 'count the nbody species',countnbody
                    else :
                        how_to_read[ispec,i]= 'staticINDEX'
                else : 
                    how_to_read[ispec,i]= 'empty'
                    field_to_index[ispec][i] = 9999999
                print 'ispec', ispec,'field',f, 'how_to_read', how_to_read[ispec,i] 
            
        # now read fields
        print "reading in particle fields"
        ip_all = 0
        for sfc in range( self.sfc_start, self.sfc_end+1 ) :
                status = artio_particle_read_root_cell_begin( self.handle, sfc,
                    num_particles_per_species )
                check_artio_status(status)	

                aoffset = np.zeros(num_acc_species,dtype="int32")
                prev_accessed_species  = 0
                for aspec in range(num_acc_species) :
                    if prev_accessed_species != accessed_species[aspec]:  # duplicates go to the same fieldtype
                        aoffset[aspec] = num_fieldnames[aspec]+aoffset[aspec]

                    ispec = iacc_to_ispec[aspec]
                    status = artio_particle_read_species_begin(self.handle, ispec);
                    check_artio_status(status)
                    
                    for particle in range( num_particles_per_species[ispec] ) :
                        
                        status = artio_particle_read_particle(self.handle,
                                        pid, subspecies, primary_variables[aspec],
                                        secondary_variables[aspec])
                        check_artio_status(status)
                        
                        pos[0] = primary_variables[aspec][pos_index[aspec][0]]
                        pos[1] = primary_variables[aspec][pos_index[aspec][1]]
                        pos[2] = primary_variables[aspec][pos_index[aspec][2]]
                        if selector.select_cell(pos,dds,eterm) :
                             for i in range(num_fieldnames[aspec]):
                                 j = i + aoffset[aspec]
                                 data[j].resize(ip_all+1)
                                 if not (how_to_read[ispec,i] == 'empty') : 
                                     assert(field_to_index[ispec][i]<100)
                                 if   how_to_read[ispec,i] == 'primary' : 
                                     data[j][ip_all] = primary_variables[aspec][field_to_index[ispec][i]]
                                 elif how_to_read[ispec,i] == 'secondary' :
                                     data[j][ip_all] = secondary_variables[aspec][field_to_index[ispec][i]]
                                 elif how_to_read[ispec,i] == 'staticNBODY' : 
                                     data[j][ip_all] = self.parameters["particle_species_mass"][field_to_index[ispec][i]]
                                 elif how_to_read[ispec,i] == 'staticINDEX' : 
                                     data[j][ip_all] = ip_all
                                 elif how_to_read[ispec,i] == 'empty' : 
                                     data[j][ip_all] = 0
                                 else : 
                                     print 'undefined how to read in caller', how_to_read[ispec,i]
                                     print 'this should be impossible.'
                                     sys.exit(1)
                                 # print 'reading into data', ip_all,data[i][ip_all], fields[i]
                             ip_all += 1
                    prev_accessed_species = accessed_species[aspec]
                    status = artio_particle_read_species_end( self.handle )
                    check_artio_status(status)
                    
                status = artio_particle_read_root_cell_end( self.handle )
                check_artio_status(status)
 

        free(subspecies)
        free(pid)
        free(num_particles_per_species)
        free(iacc_to_ispec)
        free(field_to_index)
        free(primary_variables)
        free(secondary_variables)
        print 'done filling particle variables', ip_all
        return data



    #@cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def read_grid_chunk(self, SelectorObject selector, int64_t sfc_start, int64_t sfc_end, fields ):
        cdef int i
        cdef int level
        cdef int num_oct_levels
        cdef int *num_octs_per_level
        cdef float *variables
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

        #print "reading chunk ", sfc_start, sfc_end, sfc_end-sfc_start+1

        # translate fields from ARTIO names to indices
        var_labels = self.parameters['grid_variable_labels']
        for i, f in enumerate(fields):
            if f not in var_labels:
                print "This field is not known to ARTIO:", f
                raise RuntimeError
            field_order[i] = var_labels.index(f)

        # dhr - can push these mallocs to the object to save malloc/free
        num_octs_per_level = <int *>malloc(self.max_level*sizeof(int))
        variables = <float *>malloc(8*self.num_grid_variables*sizeof(float))

        # dhr - cache the entire domain (replace later)
        status = artio_grid_cache_sfc_range( self.handle, self.sfc_file_min, self.sfc_file_max )
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
                    dpos, variables, &num_oct_levels, num_octs_per_level )
            check_artio_status(status) 

            if num_oct_levels == 0 :
                for i in range(num_fields) :
                    data[i].resize(count+1)
                    data[i][count] = variables[field_order[i]]
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

                for oct in range(num_octs_per_level[level-1]) :
                    status = artio_grid_read_oct( self.handle, dpos, variables, refined )
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
                                    data[i][count] = variables[self.num_grid_variables*child+field_order[i]]
                                count += 1 
                status = artio_grid_read_level_end( self.handle )
                check_artio_status(status) 

            status = artio_grid_read_root_cell_end( self.handle )
            check_artio_status(status) 
        
        #status = artio_grid_clear_sfc_cache( self.handle )
        #check_artio_status(status)

        free(num_octs_per_level) 
        free(variables)
        free(field_order)

        #fcoords.resize((count,3))
        #ires.resize(count)
        #    
        #for i in range(num_fields) :
        #    data[i].resize(count)

        return (fcoords, ires, data)

    def root_sfc_ranges(self, SelectorObject selector) :
        cdef int max_range_size = 1024
        cdef int coords[3]
        cdef int64_t sfc_start, sfc_end
        cdef np.float64_t pos[3]
        cdef np.float64_t dds[3]
        cdef int eterm[3]
        cdef artio_selection *selection
        cdef int i, j, k

        dds[0] = dds[1] = dds[2] = 1.0

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

        print "Selected", artio_selection_size(selection), "root cells"
 
        while artio_selection_iterator(selection, max_range_size, 
                &sfc_start, &sfc_end) == ARTIO_SUCCESS :
            sfc_ranges.append([sfc_start, sfc_end])

        print "in", len(sfc_ranges), "ranges"
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
