"""

"""
import numpy as np
cimport numpy as np
import sys 

from libc.stdint cimport int32_t, int64_t
from libc.stdlib cimport malloc, free
import  data_structures  

ROOT_LEVEL = 1 #snlroot this has to be 1 currently 

cdef struct artio_context_struct: 
    int comm
ctypedef artio_context_struct * artio_context
cdef artio_context artio_context_global

cdef extern from "sfc.h":
    int sfc_index( int coords[3], int num_root_grid_refinements ) 
    int sfc_index0( int ix, int iy, int iz, int num_root_grid_refinements ) 

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
    cdef int ARTIO_READ_REFINED_NOT_ROOT "ARTIO_READ_REFINED_NOT_ROOT"
    

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

    ctypedef void (* GridCallBackYTPos)(float * variables, int level, int refined, int64_t isfc, double pos[3], void *)
    int artio_grid_read_sfc_range_ytpos(artio_file handle,\
                int64_t sfc_min, int64_t sfc_max,\
                int min_level_to_read, int max_level_to_read,\
                int options,\
                GridCallBackYTPos callback, void *pyobject)

    ctypedef void (* GridCallBackYT)(float * variables, int level, int refined, int64_t isfc, void *)
    int artio_grid_read_sfc_range_yt(artio_file handle,\
                int64_t sfc_min, int64_t sfc_max,\
                int min_level_to_read, int max_level_to_read,\
                int options,\
                GridCallBackYT callback, void *pyobject)

    ctypedef void (* GridCallBack)(float * variables, int level, int refined,int64_t isfc)
    int artio_grid_read_sfc_range(artio_file handle,\
                int64_t sfc_min, int64_t sfc_max,\
                int min_level_to_read, int max_level_to_read,\
                int options,\
                GridCallBack callback)
    
                                 
    ctypedef void (* ParticleCallBack)(int pid, double *primary_variables,\
                                           double *secondary_variables, int species,\
                                           int subspecies, int isfc, void *pyobject )
    int artio_particle_read_sfc_range_yt(artio_file handle,\
                                          int64_t sfc_min, int64_t sfc_max,\
                                          int64_t species_min, int64_t species_max,\
                                          ParticleCallBack callback, void *pyobject)
########## wrappers calling c
cdef extern from "artio.c":
    artio_file artio_fileset_open(char * file_prefix, int type, artio_context context) 


cdef class read_parameters_artio : 
    cdef public object parameters 
    cdef artio_file handle

    def __init__(self, char *file_prefix, int artio_type) :
        self.handle = artio_fileset_open( file_prefix, artio_type, artio_context_global ) 
        self.read_parameters()
        artio_fileset_close(self.handle)  

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
        print 'print parameters in caller.pyx',self.parameters
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
        self.art_level_count = None

        self.min_level_to_read = 0
        self.max_level_to_read = param_handle.parameters['grid_max_level'][0]
        #snl FIX: the sfc values used should come from "subset" and describe the domain for chunking
        # note the root level method may force chunking to be done on 0-level ytocts 
        self.sfc_min = 0   
        self.sfc_max = param_handle.parameters['grid_file_sfc_index'][1]-1
        self.num_grid_variables = param_handle.parameters['num_grid_variables'][0]
        self.num_root_grid_refinements = int( np.log10( 
                param_handle.parameters["num_root_cells"][0]+0.5) / 
                                              (3*np.log10(2)))
        self.param_handle=param_handle
        self.ng=1 
        self.cpu=0
        self.domain_id=0

        self.grid_variable_labels=param_handle.parameters['grid_variable_labels']
 
        # the dictionary from artio variable names to yt 
        #snl FIX: this should in fields.py ?
        self.arttoyt_label_var = {}
        for label in self.grid_variable_labels :
            if   label == 'HVAR_GAS_DENSITY' :
                self.arttoyt_label_var[label] = 'Density'
            elif label == 'HVAR_GAS_ENERGY' :
                self.arttoyt_label_var[label] = 'TotalGasEnergy'
            elif label == 'HVAR_INTERNAL_ENERGY' :
                self.arttoyt_label_var[label] = 'GasEnergy'
            elif label == 'HVAR_PRESSURE' :
                self.arttoyt_label_var[label] = 'Pressure'
            elif label == 'HVAR_MOMENTUM_X' :
                self.arttoyt_label_var[label] = 'XMomentumDensity'
            elif label == 'HVAR_MOMENTUM_Y' :
                self.arttoyt_label_var[label] = 'YMomentumDensity'
            elif label == 'HVAR_MOMENTUM_Z' :
                self.arttoyt_label_var[label] = 'ZMomentumDensity'
            elif label == 'HVAR_GAMMA' :
                self.arttoyt_label_var[label] = 'Gamma'
            elif label == 'HVAR_METAL_DENSITY_II' :
                self.arttoyt_label_var[label] = 'MetalDensitySNIa'
            elif label == 'HVAR_METAL_DENSITY_Ia' :
                self.arttoyt_label_var[label] = 'MetalDensitySNII'
            elif label == 'VAR_POTENTIAL' :
                self.arttoyt_label_var[label] = 'Potential'
            elif label == 'VAR_POTENTIAL_HYDRO' :
                self.arttoyt_label_var[label] = 'PotentialHydro'
            else :
                self.arttoyt_label_var[label] = label
        print 'arttoyt_label_var (in caller.pyx):', self.arttoyt_label_var

        self.label_index = {}
        self.matched_fieldnames = []
        #snl FIX not sure about file handling:
        # self.handle = <object> handle #<------seg faults
        # and you dont bother to close all of these handles
        
    def count_refined_octs(self) :
        cdef int min_level_to_read = self.min_level_to_read
        cdef int max_level_to_read = self.max_level_to_read
        cdef int64_t sfc_min = self.sfc_min
        cdef int64_t sfc_max = self.sfc_max
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
        art_level_count = np.zeros(n_levels, dtype='int64')

        status = artio_grid_cache_sfc_range(handle, sfc_min, sfc_max)
        check_artio_status(status, artio_grid_routines.__name__)
        for sfc in xrange(sfc_min,sfc_max):
            status = artio_grid_read_root_nocts(handle, sfc,
                                                variables, num_oct_levels,
                                                num_octs_per_level)
            check_artio_status(status, artio_grid_routines.__name__)
            noct_levels = num_oct_levels[0]
            count_level_octs = {}          
            count_level_octs = [ num_octs_per_level[i] for i in xrange(noct_levels) ]
            for level in xrange(noct_levels) : 
                art_level_count[level+1] += count_level_octs[level]*8 
            num_sfc_octs = sum(count_level_octs)
            num_total_octs += num_sfc_octs
            status = artio_grid_read_root_cell_end(handle)
            check_artio_status(status, artio_grid_routines.__name__)

        if ROOT_LEVEL == 1 :    
            #add root-level octs to the count
            art_level_count[0] = (self.sfc_max+1)
            num_total_octs += (self.sfc_max+1)/8 
        self.art_level_count = art_level_count 
        print '_artio_caller.pyx num_total_octs (too big compared to what pos counts?)', num_total_octs
        return num_total_octs

    def grid_pos_fill(self, oct_handler, domain_dimensions1D) :
        ''' adds all refined octs and a new array of ghost octs for  
        the "fully refined octs" at level=-1 in ART or 0 in yt convention 
        so that root level can consist of children
        '''
        print 'start filling oct positions'
        self.oct_handler = oct_handler
        self.oct_count=0
        # fill the root grid yt-only octs 
        if ROOT_LEVEL == 1 :    
            pos = np.empty((1,3), dtype='float64')
            self.domain_dimensions = domain_dimensions1D
            for iz in  range(self.domain_dimensions/2) :
                for iy in  range(self.domain_dimensions/2) :
                    for ix in range(self.domain_dimensions/2) :
                        pos[0,0]=ix*2+1
                        pos[0,1]=iy*2+1
                        pos[0,2]=iz*2+1
                        level=0
                        self.oct_count += 1
                        self.oct_handler.add(self.cpu+1, level, 
                                             self.ng, pos, self.domain_id)
            ###################################
        # Now do real ART octs
        print 'start filling oct positions children'
        cdef artio_file handle
        handle = artio_fileset_open( self.param_handle.file_prefix, 
                                     ARTIO_OPEN_GRID, artio_context_global ) 
        status = artio_grid_read_sfc_range_ytpos(handle,\
                    self.sfc_min, self.sfc_max,\
                    self.min_level_to_read, self.max_level_to_read,\
                    ARTIO_READ_REFINED,\
                    wrap_oct_pos_callback, <void*>self) 
        check_artio_status(status, artio_grid_routines.__name__)
        artio_fileset_close(handle) 
        print 'done filling oct positions; allocated octs:', self.oct_count
        # snl FIX assert oct_count matches num octs elsewhere
    def grid_var_fill(self, source, fields):
        print 'start filling grid vars the root grid fill takes too long...'
        self.source = source
        i=-1
        for artlabel in self.grid_variable_labels :
            label = self.arttoyt_label_var[artlabel]
            i=i+1
            for field in fields : 
                if field == label :
                    print 'match in fields?',artlabel, field,label, fields
                    self.label_index[field]=i
                    self.matched_fieldnames.append(field)
        print 'matched fields:',self.matched_fieldnames
        print 'art index of matched fields',self.label_index

        
#        i=-1
#        for field in fields : 
#            i+=1
#            print field
#            self.matched_fieldnames.append(field)
#            self.label_index[field] = i
#        print 'quitting in caller'
#        sys.exit(1)
        self.count=0
        #fill grid positions with values
        cdef artio_file handle
        
        if ROOT_LEVEL == 1 :    
            coords = np.empty(3, dtype='int64')
            child = np.empty(3, dtype='int64')
            if len(self.label_index) > 0 :
                # start with root run over root-level yt-only octs 
                # where yt-octs live on level 0 
                handle = artio_fileset_open( self.param_handle.file_prefix, 
                                             ARTIO_OPEN_GRID, artio_context_global ) 
                for iz in  range(self.domain_dimensions/2) :
                    for iy in  range(self.domain_dimensions/2) :
                        for ix in range(self.domain_dimensions/2) :
                            coords[0]=ix*2
                            coords[1]=iy*2
                            coords[2]=iz*2
                            for k in range(2):
                                for j in range(2):
                                    for i in range(2):
                                        child[0] = coords[0]+i
                                        child[1] = coords[1]+j
                                        child[2] = coords[2]+k
                                        isfc = sfc_index0(child[0], child[1], child[2], self.num_root_grid_refinements)
                                        status = artio_grid_read_sfc_range_yt(
                                            handle, isfc, isfc,\
                                                0, 0,\
                                                ARTIO_READ_ALL,\
                                                wrap_cell_var_callback,\
                                                <void*>self
                                            ) #only reads root
                                        check_artio_status(status, artio_grid_routines.__name__)
                artio_fileset_close(handle) 
                print 'snlroot done with root var fill'
            #now run through oct children
            handle = artio_fileset_open( self.param_handle.file_prefix, 
                                         ARTIO_OPEN_GRID, artio_context_global ) 
            status = artio_grid_read_sfc_range_yt(
                handle, self.sfc_min, self.sfc_max,\
                    self.min_level_to_read+1, self.max_level_to_read,\
                    ARTIO_READ_ALL,\
                    wrap_cell_var_callback,\
                    <void*>self
                ) #only octs!
            check_artio_status(status, artio_grid_routines.__name__)
            artio_fileset_close(handle) 
        print 'done buffering variables'
    def oct_pos_callback(self, level, refined, isfc, pos):
#        print 'callerpos ',self.oct_count*8,pos[0,0],pos[0,1],pos[0,2],vars, level
        self.oct_count += 1
        self.oct_handler.add(self.cpu+1, level-self.min_level_to_read, 
                             self.ng, pos, self.domain_id)
    def cell_var_callback(self, level, refined, ichild, cell_var):
        for field in self.matched_fieldnames : 
            self.source[field][self.count] = cell_var[self.label_index[field]] 
#        if (ichild == 0) and (level>0):
#            print 'callervar ',self.count,cell_var[0],level
        self.count += 1
 

###### callbacks (e.g. https://github.com/cython/cython/tree/master/Demos/callback) ######
cdef void wrap_oct_pos_callback(float *variables, int level, int refined, 
                                 int64_t isfc, double *pos, void *pyobject):
    position = np.empty((1, 3), dtype='float64')
    position[0,0] = pos[0]
    position[0,1] = pos[1]
    position[0,2] = pos[2]
    # add one to level because in yt, octs live on the same level as their children
    # 0-level ART octs do not exist in memory (the ART root cells are not children)
    ytlevel = level+1 
    artioroutines = <object>pyobject
    #    print '_artio_caller.pyx:octpositionsandvalues',ytlevel, pos[0],pos[1],pos[2],level,variables[0]
    artioroutines.oct_pos_callback(ytlevel, refined, isfc, position) #variables[0]

cdef void wrap_cell_var_callback(float *variables, int level, int refined, 
                                 int64_t ichild, void *pyobject):
    artioroutines = <object>pyobject
    cell_var={}
    cell_var = [ variables[i] for i in range(artioroutines.num_grid_variables) ]
    #    if level > 0:
    #        print '_artio_caller.pyx:sourcecallvalue', level, cell_var
    artioroutines.cell_var_callback(level, refined, ichild, cell_var)

class artio_particle_routines(object) : 
    def __init__(self, param_handle) :

        self.min_level_to_read = 0
        self.max_level_to_read = param_handle.parameters['grid_max_level'][0]
        # sfc values should come from a subset object and describe the chunked domains
        self.sfc_min = 0   
        self.sfc_max = param_handle.parameters['grid_file_sfc_index'][1]-1
        self.param_handle=param_handle
        self.ng=1 #RISM
        self.cpu=0 #RISM
        self.domain_id=0 #RISM

        # FIX: too much is hardcoded here, much of this should probably live in fields.py
        # choose particle types (species labels), species (self.num_particle_species) 
        # and fields (variable_labels) to be read.
        #     particle_species_labels: N-BODY STAR
        #     variable_labels: POSITION_X, VELOCITY_X, TIMESTEP
        self.num_particle_species = param_handle.parameters['num_particle_species'][0]
        self.num_primary_variables = param_handle.parameters['num_primary_variables']
        self.num_secondary_variables = param_handle.parameters['num_secondary_variables']
        # sfc values should come from subset and describe the domain for chunking
        self.species_min = 0   
        self.species_max = self.num_particle_species
 
        self.particle_variable_labels = {}
        # combine primary and secondary varibles into one particle_variable labels 
        # organized by species count
        self.num_particle_variable_labels=[]
        for ispec in xrange(self.num_particle_species) :
            speclabelprimary = "species_%02d_primary_variable_labels" % ispec           
            self.particle_variable_labels[ispec] = param_handle.parameters[speclabelprimary]
            if self.num_secondary_variables[ispec] > 0 :
                speclabelsecondary = "species_%02d_secondary_variable_labels" % ispec           
                self.particle_variable_labels[ispec] = param_handle.parameters[speclabelprimary]+\
                    param_handle.parameters[speclabelsecondary]
                #self.num_particle_variable_labels[ispec]=self.num_primary_variables[ispec]+self.num_secondary_variables[ispec]
                
        # ####### variable dictionary from artio to yt ############
        self.arttoyt_label_spec = {}
        self.particle_species_labels = param_handle.parameters['particle_species_labels']
        for label in self.particle_species_labels : 
            print 'particle labels in caller.pyx',label
            if label == 'N-BODY' : 
                self.arttoyt_label_spec[label] = 'nbody' 
            elif label == 'STARS' : 
                self.arttoyt_label_spec[label] = 'stars' 
            else :
                self.arttoyt_label_spec[label] = label
        self.arttoyt_label_var = {}
        # 
        for ispec in xrange(self.num_particle_species):
            for label in self.particle_variable_labels[ispec] :
                if label == 'BIRTH_TIME' :
                    self.arttoyt_label_var[(ispec,label)] = 'creation_time'
                else :
                    self.arttoyt_label_var[(ispec,label)] = label
                    print (ispec,label),label
        print '----------------- arttoyt_label_spec (caller.pyx)'
        print self.arttoyt_label_spec
        print '----------------- arttoyt_label_var (caller.pyx)'
        print self.arttoyt_label_var
        print '-----------------'
        print 'particle arttoyt translation'

    def particle_pos_fill(self,accessed_species, fieldnames) :
        #snl FIX the following should go in fields or be a subclass
        #what are the matched art species? 
        art_isubspec={}
        art_labelispec={}
        self.art_ispec={}
        ispec=-1 
        for artlabelspec in self.particle_species_labels : #N-BODY, N-BODY, STAR
            ispec+=1
            for ytlabelspec in accessed_species : 
                if self.arttoyt_label_spec[artlabelspec] == ytlabelspec : 
                    print 'match, in species labels?', ytlabelspec, accessed_species
                    if art_isubspec[ytlabelspec] == None :
                        isubspec = 0
                    else : 
                        isubspec += 1
                    self.art_ispec.append(ispec)
                    art_isubspec[ytlabelspec].append(isubspec)
                    art_labelispec[ytlabelspec].append(ispec)
#                    self.matched_labelspec.append(ytlabelspec)
                    
        # FIX Hardcoded subspecies selection
        # e.g., only read least massive N-BODY particles -- pop off the species corresponding 
        # to the subspecies you don't want
        min_mass = 1e30
        imin_mass = 0
        for isubspec in art_isubspec['nbody'] : 
            ispec = art_labelispec['nbody'][isubspec]
            subspecmass = self.param_handle.parameters['particle_species_mass'][ispec] 
            if subspecmass < min_mass :
                min_mass = subspecmass
                imin_mass = isubspec
        # pop elements in nbody that are not min_mass
        for isubspec in art_isubspec['nbody'] : 
            ispec = art_labelispec['nbody'][isubspec]
            if isubspec != imin_mass :
                art_isubspec['nbody'].pop(ispec)
                self.art_ispec.pop(ispec)
                art_labelispec['nbody'].pop(ispec)
        print 'art_isubspec[nbody] after pop...', art_isubspec['nbody']
        print self.art_ispec
        print art_labelispec['nbody'].pop(ispec)

        print 'exiting in particle_pos_fill' 
        sys.exit(1)
         
        #what are the matched art variable indices and (yt) labels for the matched species?
        self.art_ivar=[]
        self.matched_labelvar=[]
        for ispec in self.art_ispec :
            self.matched_labelvar[ispec]={}
            i=-1
            artivar=-1
            for artlabelvar in self.particle_variable_labels[ispec] :  
            #POSITION_X,... BIRTH_TIME, INITIAL_MASS, MASS, METALLICITY_SNII, METALLICITY_SNIa
                artivar+=1
                i+=1
                for ytlabelvar in fieldnames : #fieldnames could be a function of species
                    if self.arttoyt_label_var[artlabelvar] == ytlabelvar :
                        print 'match, in fieldnames?',ispec, ytlabelvar, fieldnames
                        self.art_ivar[(ispec,ytlabelvar)]=i
#                        self.art_ivar[(ispec,artivar)]=i
                        self.matched_labelvar[ispec].append(ytlabelvar)
            print 'ispec=',ispec,'matched variable labels', self.matched_labelvar[ispec]
            print 'art index of matched fieldnames',self.art_ivar
        print 'exiting in particle_pos_fill' 
        sys.exit(1)
################ The above should be part of a subclass that uses accessed_species #######       

        cdef artio_file handle
        self.accessed_species = accessed_species
        
        self.selection={}
        #ouch we are reading the particle file nspecies times?!
        for ispec in self.art_ispec : 
#              ispec->      self.species_min, self.species_max,\ 
            handle = artio_fileset_open( self.param_handle.file_prefix, 
                                         ARTIO_OPEN_PARTICLES, artio_context_global ) 
            status = artio_particle_read_sfc_range_yt(handle,\
                    self.sfc_min, self.sfc_max,\
                    ispec, ispec,\
                    wrap_particle_pos_callback, <void*>self) 
            check_artio_status(status, artio_grid_routines.__name__)
            artio_fileset_close(handle) 
            
        print 'done reading particle positions'
        return self.selection
    def particle_var_fill(self, source, fieldnames, accessed_species, particle_mask):
        cdef artio_file handle
        self.source = source
        if len(self.art_ivar) > 0 :
            self.count=0
            for ispec in self.art_ispec : 
                handle = artio_fileset_open( self.param_handle.file_prefix, 
                                             ARTIO_OPEN_PARTICLES, artio_context_global ) 
                #snl check this:
                status = artio_particle_read_sfc_range_yt(handle,\
                        self.sfc_min, self.sfc_max,\
                        self.species_min, self.species_max,\
                        wrap_particle_var_callback, <void*>self)
                check_artio_status(status, artio_grid_routines.__name__)
                artio_fileset_close(handle) 
        print 'done buffering variables'
    def particle_pos_callback(self, pos):
        self.selection['x'].append(pos[0])
        self.selection['y'].append(pos[1])
        self.selection['z'].append(pos[2])
    def particle_var_callback(self, particle_var, ispec):
        for labelvar in self.matched_labelvar[ispec] : 
            artivar = self.art_ivar[(ispec,labelvar)]
            self.source[labelvar][self.count] = particle_var[artivar] 
        self.count=self.count+1
 
#        callback(pid,
#                 primary_variables,
#                 secondary_variables,
#                 species, subspecies, sfc)
###### callbacks (e.g. https://github.com/cython/cython/tree/master/Demos/callback) ######
cdef void wrap_particle_pos_callback(int pid, double *primary_variables, 
                                     double *secondary_variables, int species, 
                                     int subspecies, int isfc, void *pyobject ) :
    artioroutines = <object>pyobject
    #FIX hardcoded assumption that position always comes first in primary variables
    #assert this elsewhere
    pos=[]
    pos[0:2]=[primary_variables[0],primary_variables[1],primary_variables[2]]
    artioroutines.particle_pos_callback(pos)

cdef void wrap_particle_var_callback(int pid, double *primary_variables, double *secondary_variables, int species, 
                                 int subspecies, int isfc, void *pyobject ) :
    artioroutines = <object>pyobject
    cell_var={}
    # subspecies is only relevant when there are multiple "types" of "star-particles"  (e.g. AGN) 
    ispec = species
    particle_var = [ primary_variables[i] for i in range(artioroutines.num_primary_variables[ispec]) ]\
        +        [ secondary_variables[i] for i in range(artioroutines.num_secondary_variables[ispec]) ]
    print 'snl count particle var',(artioroutines.num_primary_variables[ispec]+artioroutines.num_secondary_variables[ispec]), particle_var 
    artioroutines.cell_var_callback(particle_var)
        
###################################################
def artio_is_valid( char *file_prefix ) :
    cdef artio_file handle = artio_fileset_open( file_prefix, 
            ARTIO_OPEN_HEADER, artio_context_global )
    if handle == NULL :
        return False
    else :
        artio_fileset_close(handle) 
    return True

