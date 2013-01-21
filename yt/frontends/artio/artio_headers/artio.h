/*
 * artio.h
 *
 *  Created on: Feb 21, 2010
 *      Author: Yongen Yu
 *  Modified: Jun 6, 2010 - Doug Rudd
 *            Nov 18, 2010 - Doug Rudd
 *            Nov 14, 2012 - Doug Rudd
 */

#ifndef __ARTIO_H__
#define __ARTIO_H__

#include <stdint.h>
#ifndef int64_t
#ifdef _WIN32
typedef __int64 int64_t;
#endif
#else
#error "Undefined int64_t!"
#endif

#define ARTIO_OPEN_HEADER					0
#define ARTIO_OPEN_PARTICLES                1
#define ARTIO_OPEN_GRID                     2

#define ARTIO_READ_LEAFS                    1
#define ARTIO_READ_REFINED                  2
#define	ARTIO_READ_ALL                      3
#define	ARTIO_READ_REFINED_AND_ROOT         4

/* allocation strategy */
#define ARTIO_ALLOC_EQUAL_SFC               0
#define ARTIO_ALLOC_EQUAL_PROC              1
#define ARTIO_ALLOC_MAX_FILE_SIZE  	    2

#define ARTIO_TYPE_STRING                   0
#define ARTIO_TYPE_CHAR                     1
#define ARTIO_TYPE_INT                      2
#define ARTIO_TYPE_FLOAT                    3
#define ARTIO_TYPE_DOUBLE                   4
#define ARTIO_TYPE_LONG                     5

/* error codes */
#define ARTIO_SUCCESS                       0

#define ARTIO_ERR_PARAM_NOT_FOUND           1
#define ARTIO_ERR_PARAM_INVALID_LENGTH      2
#define ARTIO_ERR_PARAM_TYPE_MISMATCH       3
#define ARTIO_ERR_PARAM_LENGTH_MISMATCH     4
#define ARTIO_ERR_PARAM_LENGTH_INVALID      5
#define ARTIO_ERR_PARAM_DUPLICATE           6

#define ARTIO_ERR_INVALID_FILESET_MODE      100
#define	ARTIO_ERR_INVALID_FILE_NUMBER       101
#define ARTIO_ERR_INVALID_FILE_MODE         102
#define ARTIO_ERR_INVALID_SFC_RANGE         103
#define ARTIO_ERR_INVALID_SFC               104
#define ARTIO_ERR_INVALID_STATE             105
#define ARTIO_ERR_INVALID_SEEK              106
#define ARTIO_ERR_INVALID_OCT_LEVELS        107
#define ARTIO_ERR_INVALID_SPECIES           108
#define ARTIO_ERR_INVALID_ALLOC_STRATEGY    109
#define ARTIO_ERR_INVALID_LEVEL             110
#define ARTIO_ERR_INVALID_PARAM_LIST        111
#define ARTIO_ERR_INVALID_DATATYPE          112
#define ARTIO_ERR_INVALID_OCT_REFINED       113
#define ARTIO_ERR_INVALID_HANDLE            114

#define ARTIO_ERR_DATA_EXISTS               200
#define ARTIO_ERR_INSUFFICIENT_DATA         201
#define ARTIO_ERR_FILE_CREATE               202
#define ARTIO_ERR_PARTICLE_FILE_NOT_FOUND   203
#define ARTIO_ERR_GRID_FILE_NOT_FOUND       204

#define ARTIO_ERR_PARAM_CORRUPTED           207
#define ARTIO_ERR_PARAM_CORRUPTED_MAGIC     208

#define ARTIO_ERR_64_TO_32_BIT_TRUNCATION   209
#define ARTIO_ERR_MEMORY_ALLOCATION         210

#define ARTIO_PARAMETER_EXHAUSTED           300

#ifdef ARTIO_MPI
#include <mpi.h>

struct artio_context_struct {
    MPI_Comm comm;
};
#else
struct artio_context_struct {
    int comm;
};
#endif

typedef struct artio_file_struct * artio_file;
typedef struct artio_param_list * artio_parameters;
typedef struct artio_context_struct * artio_context;

extern artio_context artio_context_global;

/*
 * Description: Open the file
 *
 *  filename			The file prefix
 *  type			combination of ARTIO_OPEN_PARTICLES and ARTIO_OPEN_GRID flags
 */
artio_file artio_fileset_open(char * file_name, int type, artio_context context);

/**
 * Description: Create fileset and begin populating header information
 *
 *  file_name			file name of refined cells
 *  root_cells			the number of root level cells
 *  proc_sfc_begin-end		the range of local space-filling-curve indices
 *  handle			the artio file handle
 *
 */
artio_file artio_fileset_create(char * file_prefix, 
        int64_t root_cells, int64_t proc_sfc_begin, int64_t proc_sfc_end, artio_context context);

/*
 * Description	Close the file
 */
int artio_fileset_close(artio_file handle);

/* public parameter interface */
int artio_parameter_iterate( artio_file handle, char *key, int *type, int *length );
int artio_parameter_get_array_length(artio_file handle, char * key, int *length);

void artio_parameter_set_int(artio_file handle, char * key, int32_t value);
int artio_parameter_get_int(artio_file handle, char * key, int32_t * value);

void artio_parameter_set_int_array(artio_file handle, char * key, int length,
		int32_t *values);
int artio_parameter_get_int_array(artio_file handle, char * key, int length,
		int32_t *values);

void artio_parameter_set_string(artio_file handle, char * key, char * value);
int artio_parameter_get_string(artio_file handle, char * key, char * value, int max_length);

void artio_parameter_set_string_array(artio_file handle, char * key,
		int length, char ** values);
int artio_parameter_get_string_array(artio_file handle, char * key,
		int length, char ** values, int max_length);

void artio_parameter_set_float(artio_file handle, char * key, float value);
int artio_parameter_get_float(artio_file handle, char * key, float * value);

void artio_parameter_set_float_array(artio_file handle, char * key,
		int length, float *values);
int artio_parameter_get_float_array(artio_file handle, char * key,
		int length, float * values);

void artio_parameter_set_double(artio_file handle, char * key, double value);
int  artio_parameter_get_double(artio_file handle, char * key, double * value);

void artio_parameter_set_double_array(artio_file handle, char * key,
		int length, double * values);
int artio_parameter_get_double_array(artio_file handle, char * key,
        int length, double *values);

void artio_parameter_set_long(artio_file handle, char * key, int64_t value);
int artio_parameter_get_long(artio_file handle, char * key, int64_t *value);

void artio_parameter_set_long_array(artio_file handle, char * key,
        int length, int64_t *values);
int artio_parameter_get_long_array(artio_file handle, char * key,
        int length, int64_t *values);

/* public grid interface */
typedef void (* GridCallBack)(float * variables, int level, int refined,
		int64_t sfc_index);
typedef void (* GridCallBackYTPos)(double * variables, int level, int refined,
                                 int64_t sfc_index, double pos[3], void * pyobject);
typedef void (* GridCallBackYT)(double * variables, int level, int refined,
                                 int64_t sfc_index, void * pyobject);

/*
 * Description:	Add a grid component to a fileset open for writing
 *
 *  handle			The fileset handle
 *  num_grid_files		The number of grid files to create
 *  allocation_strategy		How to apportion sfc indices to each grid file
 *  num_grid_variables		The number of variables per cell
 *  grid_variable_labels	Identifying labels for each variable
 *  num_levels_per_root_tree	Maximum tree depth for each oct tree
 *  num_octs_per_root_tree	Total octs in each oct tree
 */
int artio_fileset_add_grid(artio_file handle,
        int num_grid_files, int allocation_strategy,
        int num_grid_variables,
        char ** grid_variable_labels,
        int * num_levels_per_root_tree,
        int * num_octs_per_root_tree );

int artio_fileset_open_grid(artio_file handle);
int artio_fileset_close_grid(artio_file handle);
int artio_fileset_open_particle(artio_file handle);
int artio_fileset_close_particle(artio_file handle);

/*
 * Description:	Output the variables of the root level cell and the hierarchy of the Oct tree correlated with this root level cell
 *
 *  handle			The File handle
 *  sfc				The sfc index of root cell
 *  variables			The variables of the root level cell
 *  level			The depth of the Oct tree correlated to the root level cell
 *  num_level_octs		The array store the number of Oct nodes each level
 */
int artio_grid_write_root_cell_begin(artio_file handle, int64_t sfc, 
		float * variables, int level, int * num_octs_per_level);

/*
 * Description:	Do something at the end of writing the root level cell
 */
int artio_grid_write_root_cell_end(artio_file handle);

/*
 * Description:	Do something at the beginning of each level
 */
int artio_grid_write_level_begin(artio_file handle, int level );

/*
 * Description:	Do something at the end of each level
 */
int artio_grid_write_level_end(artio_file handle);

/*
 * Description:	Output the data of a special oct tree node to the file
 *
 *  handle			The handle of the file
 *  variables 			The array recording the variables of the eight cells belonging to this Octree node.
 */
int artio_grid_write_oct(artio_file handle, float *variables, int *refined);

/*
 * Description:	Read the variables of the root level cell and the hierarchy of the Octtree
 *              correlated with this root level cell
 *
 *  handle			The File handle
 *  variables			The variables of the root level cell
 *  level 			The depth of the OCT tree
 *  num_octs_per_level		The number of node of each oct level
 *
 */
int artio_grid_read_root_nocts(artio_file handle, int64_t sfc, float *variables,
		int *num_tree_levels, int *num_octs_per_level);
int artio_grid_read_root_cell_begin(artio_file handle, int64_t sfc, float *variables,
		int *num_tree_levels, int *num_octs_per_level);

/*
 * Description:	Do something at the end of reading the root level cell
 */
int artio_grid_read_root_cell_end(artio_file handle);

/*
 * Description:	Do something at the beginning of each level
 */
int artio_grid_read_level_begin(artio_file handle, int level );

/*
 * Description:	Do something at the end of each level
 */
int artio_grid_read_level_end(artio_file handle);

/*
 * Description:	Read the data of a special oct tree node from the file
 */
int artio_grid_read_oct(artio_file handle, float *variables, int *refined);

int artio_grid_cache_sfc_range(artio_file handle, int64_t sfc_start, int64_t sfc_end);

/*
 * Description:	Read a segment of oct nodes
 *
 *  handle			file pointer
 *  sfc1			the start sfc index
 *  sfc2			the end sfc index
 *  max_level_to_read		max level to read for each oct tree
 *  option			1. refined nodes; 2 leaf nodes; 3 all nodes
 *  callback			callback function
 */
int artio_grid_read_sfc_range(       artio_file handle, int64_t sfc1, int64_t sfc2, int min_level_to_read, int max_level_to_read, int options, GridCallBack callback);
int artio_grid_read_sfc_range_ytpos(artio_file handle, int64_t sfc1, int64_t sfc2, int min_level_to_read,int max_level_to_read, int options, GridCallBackYTPos callback, void *pyobject);
int artio_grid_read_sfc_range_yt(artio_file handle, int64_t sfc1, int64_t sfc2, int min_level_to_read, int max_level_to_read, int options, GridCallBackYT callback, void *pyobject);
		

typedef void (* ParticleCallBack)(int64_t pid, 
		double *primary_variables, float *secondary_variables, 
		int species, int subspecies, int64_t sfc_index);
typedef void (* ParticleCallBackYT)(int64_t pid, 
		double *primary_variables, float *secondary_variables, 
                int species, int subspecies, int64_t sfc_index, void *pyobject);

/**
 *  header			head file name
 *  num_particle_files		the number of files to record refined cells
 *  allocation_strategy
 *  num_species			number of particle species
 *  species_labels		string identifier for each species
 *  handle			the artio file handle
 *
 */
int artio_fileset_add_particles(artio_file handle, 
        int num_particle_files, int allocation_strategy,
        int num_species, char **species_labels,
        int *num_primary_variables,
        int *num_secondary_variables,
        char ***primary_variable_labels_per_species,
        char ***secondary_variable_labels_per_species,
        int *num_particles_per_species_per_root_tree );

int artio_fileset_open_particles(artio_file handle);
int artio_fileset_close_particles(artio_file handle);

/*
 * Description:	Output the variables of the root level cell and the hierarchy of 
 *                  the oct-tree correlated with this root level cell
 *
 *  handle			The File handle
 *  sfc				The sfc index of root cell
 *  variables			The variables of the root level cell
 *  level			The depth of the Oct tree correlated to the root level cell
 *  num_level_octs		The array store the number of Oct nodes each level
 */
int artio_particle_write_root_cell_begin(artio_file handle, int64_t sfc,
		int *num_particles_per_species);

/*
 * Description:	Do something at the end of writing the root level cell
 */
int artio_particle_write_root_cell_end(artio_file handle);

/*
 * Description:	Do something at the beginning of each level
 */
int artio_particle_write_species_begin(artio_file handle, int species );

/*
 * Description:	Do something at the end of each level
 */
int artio_particle_write_species_end(artio_file handle);

/*
 * Description: Output the data of a special oct tree node to the file
 *
 *  handle			The handle of the file
 *  variables 			The array recording the variables of the eight cells belonging to this Octree node.
 */
int artio_particle_write_particle(artio_file handle, int64_t pid, int subspecies, 
			double* primary_variables, float *secondary_variables);

/*
 * Description:	Read the variables of the root level cell and the hierarchy of the Octtree
 *              correlated with this root level cell
 *
 *  handle			The File handle
 *  variables			The variables of the root level cell
 *  level 			The depth of the OCT tree
 *  num_octs_per_level		The number of node of each oct level
 *
 */
int artio_particle_read_root_cell_begin(artio_file handle, int64_t sfc, 
			int * num_particle_per_species);

/*
 * Description:	Do something at the end of reading the root level cell
 */
int artio_particle_read_root_cell_end(artio_file handle);

/*
 * Description:	Do something at the beginning of each level
 */
int artio_particle_read_species_begin(artio_file handle, int species );

/*
 * Description:  Do something at the end of each level
 */
int artio_particle_read_species_end(artio_file handle);

/*
 * Description:	Read the data of a single particle from the file
 */
int artio_particle_read_particle(artio_file handle, int64_t *pid, int *subspecies,
			double *primary_variables, float *secondary_variables);

int artio_particle_cache_sfc_range(artio_file handle, int64_t sfc_start, int64_t sfc_end);

/*
 * Description: Read a segment of particles
 *
 *  handle			file pointer
 *  sfc1			the start sfc index
 *  sfc2			the end sfc index
 *  start_species		the first particle species to read
 *  end_species			the last particle species to read
 *  callback			callback function
 */
int artio_particle_read_sfc_range(artio_file handle, 
		int64_t sfc1, int64_t sfc2, 
		int start_species, int end_species,
		ParticleCallBack callback);

int artio_particle_read_sfc_range_yt(artio_file handle, 
		int64_t sfc1, int64_t sfc2, 
		int start_species, int end_species,
                ParticleCallBackYT callback, void *pyobject);

#endif /* __ARTIO_H__ */
