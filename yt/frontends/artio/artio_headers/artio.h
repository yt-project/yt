/**********************************************************************
 * Copyright (c) 2012-2013, Douglas H. Rudd
 * All rights reserved.
 *
 * This file is part of the artio library.
 *
 * artio is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * artio is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * Copies of the GNU Lesser General Public License and the GNU General
 * Public License are available in the file LICENSE, included with this
 * distribution.  If you failed to receive a copy of this file, see
 * <http://www.gnu.org/licenses/>
 **********************************************************************/

#ifndef __ARTIO_H__
#define __ARTIO_H__

#define ARTIO_MAJOR_VERSION     1
#define ARTIO_MINOR_VERSION     2

#ifdef ARTIO_MPI
#include <mpi.h>
#endif

#ifdef _WIN32
typedef __int64 int64_t;
typedef __int32 int32_t;
#else
#include <stdint.h>
#endif

#define ARTIO_OPEN_HEADER					0
#define ARTIO_OPEN_PARTICLES                1
#define ARTIO_OPEN_GRID                     2

#define ARTIO_READ_LEAFS                    1
#define ARTIO_READ_REFINED                  2
#define	ARTIO_READ_ALL                      3

#define ARTIO_RETURN_OCTS					4
#define ARTIO_RETURN_CELLS					0

/* allocation strategy */
#define ARTIO_ALLOC_EQUAL_SFC               0
#define ARTIO_ALLOC_EQUAL_PROC              1
#define ARTIO_ALLOC_MAX_FILE_SIZE  	        2

#define ARTIO_TYPE_STRING                   0
#define ARTIO_TYPE_CHAR                     1
#define ARTIO_TYPE_INT                      2
#define ARTIO_TYPE_FLOAT                    3
#define ARTIO_TYPE_DOUBLE                   4
#define ARTIO_TYPE_LONG                     5

/* error codes */
#define ARTIO_SUCCESS                       0

#define ARTIO_ERR_PARAM_NOT_FOUND           1
#define ARTIO_PARAMETER_EXHAUSTED			2
#define ARTIO_ERR_PARAM_INVALID_LENGTH      3 
#define ARTIO_ERR_PARAM_TYPE_MISMATCH       4
#define ARTIO_ERR_PARAM_LENGTH_MISMATCH     5
#define ARTIO_ERR_PARAM_LENGTH_INVALID      6
#define ARTIO_ERR_PARAM_DUPLICATE           7
#define ARTIO_ERR_PARAM_CORRUPTED           8
#define ARTIO_ERR_PARAM_CORRUPTED_MAGIC     9
#define ARTIO_ERR_STRING_LENGTH             10

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
#define ARTIO_ERR_INVALID_PARAMETER_LIST    111
#define ARTIO_ERR_INVALID_DATATYPE          112
#define ARTIO_ERR_INVALID_OCT_REFINED       113
#define ARTIO_ERR_INVALID_HANDLE            114
#define ARTIO_ERR_INVALID_CELL_TYPES        115
#define ARTIO_ERR_INVALID_BUFFER_SIZE		116
#define ARTIO_ERR_INVALID_INDEX				117

#define ARTIO_ERR_DATA_EXISTS               200
#define ARTIO_ERR_INSUFFICIENT_DATA         201
#define ARTIO_ERR_FILE_CREATE               202
#define ARTIO_ERR_GRID_DATA_NOT_FOUND       203
#define ARTIO_ERR_GRID_FILE_NOT_FOUND       204
#define ARTIO_ERR_PARTICLE_DATA_NOT_FOUND   205
#define ARTIO_ERR_PARTICLE_FILE_NOT_FOUND   206
#define ARTIO_ERR_IO_OVERFLOW               207
#define ARTIO_ERR_IO_WRITE                  208
#define ARTIO_ERR_IO_READ                   209
#define ARTIO_ERR_BUFFER_EXISTS             210

#define ARTIO_SELECTION_EXHAUSTED           300
#define ARTIO_ERR_INVALID_SELECTION         301
#define ARTIO_ERR_INVALID_COORDINATES       302

#define ARTIO_ERR_MEMORY_ALLOCATION         400

#define ARTIO_ERR_VERSION_MISMATCH			500

#ifdef ARTIO_MPI
typedef struct {
    MPI_Comm comm;
} artio_context;
#else
typedef struct {
    int comm;
} artio_context;
#endif

#define ARTIO_MAX_STRING_LENGTH				256

typedef struct artio_fileset_struct artio_fileset;
typedef struct artio_selection_struct artio_selection;

extern const artio_context *artio_context_global;

/*
 * Description: Open the file
 *
 *  filename		The file prefix
 *  type			combination of ARTIO_OPEN_PARTICLES and ARTIO_OPEN_GRID flags
 */
artio_fileset *artio_fileset_open( char * file_name, int type, const artio_context *context);

/**
 * Description: Create fileset and begin populating header information
 *
 *  file_name			file name of refined cells
 *  root_cells			the number of root level cells
 *  proc_sfc_begin-end		the range of local space-filling-curve indices
 *  handle			the artio file handle
 *
 */
artio_fileset *artio_fileset_create(char * file_prefix, 
        int64_t root_cells, int64_t proc_sfc_begin, int64_t proc_sfc_end, const artio_context *context);

/*
 * Description	Close the file
 */
int artio_fileset_close(artio_fileset *handle);
int artio_fileset_set_buffer_size( int buffer_size );
int artio_fileset_has_grid( artio_fileset *handle );
int artio_fileset_has_particles( artio_fileset *handle );

/* public parameter interface */
int artio_parameter_iterate( artio_fileset *handle, char *key, int *type, int *length );
int artio_parameter_get_array_length(artio_fileset *handle, const char * key, int *length);

int artio_parameter_set_int(artio_fileset *handle, const char * key, int32_t value);
int artio_parameter_get_int(artio_fileset *handle, const char * key, int32_t * value);

int artio_parameter_set_int_array(artio_fileset *handle, const char * key, int length,
		int32_t *values);
int artio_parameter_get_int_array(artio_fileset *handle, const char * key, int length,
		int32_t *values);
int artio_parameter_get_int_array_index(artio_fileset *handle, const char * key, 
		int index, int32_t *values);

int artio_parameter_set_string(artio_fileset *handle, const char * key, char * value);
int artio_parameter_get_string(artio_fileset *handle, const char * key, char * value );

int artio_parameter_set_string_array(artio_fileset *handle, const char * key,
		int length, char ** values);
int artio_parameter_get_string_array(artio_fileset *handle, const char * key,
		int length, char ** values );
int artio_parameter_get_string_array_index(artio_fileset *handle, const char * key,
		int index, char * values );

int artio_parameter_set_float(artio_fileset *handle, const char * key, float value);
int artio_parameter_get_float(artio_fileset *handle, const char * key, float * value);

int artio_parameter_set_float_array(artio_fileset *handle, const char * key,
		int length, float *values);
int artio_parameter_get_float_array(artio_fileset *handle, const char * key,
		int length, float * values);
int artio_parameter_get_float_array_index(artio_fileset *handle, const char * key,
		int index, float * values);

int artio_parameter_set_double(artio_fileset *handle, const char * key, double value);
int  artio_parameter_get_double(artio_fileset *handle, const char * key, double * value);

int artio_parameter_set_double_array(artio_fileset *handle, const char * key,
		int length, double * values);
int artio_parameter_get_double_array(artio_fileset *handle, const char * key,
        int length, double *values);
int artio_parameter_get_double_array_index(artio_fileset *handle, const char * key,
		int index, double *values);

int artio_parameter_set_long(artio_fileset *handle, const char * key, int64_t value);
int artio_parameter_get_long(artio_fileset *handle, const char * key, int64_t *value);

int artio_parameter_set_long_array(artio_fileset *handle, const char * key,
        int length, int64_t *values);
int artio_parameter_get_long_array(artio_fileset *handle, const char * key,
        int length, int64_t *values);
int artio_parameter_get_long_array_index(artio_fileset *handle, const char * key,
		int index, int64_t *values);

/* public grid interface */
typedef void (* artio_grid_callback)( int64_t sfc_index, int level,
		double *pos, float * variables, int *refined, void *params );

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
int artio_fileset_add_grid(artio_fileset *handle,
        int num_grid_files, int allocation_strategy,
        int num_grid_variables,
        char ** grid_variable_labels,
        int * num_levels_per_root_tree,
        int * num_octs_per_root_tree );

int artio_fileset_open_grid(artio_fileset *handle);
int artio_fileset_close_grid(artio_fileset *handle);

/*
 * Description:	Output the variables of the root level cell and the index of the Oct tree correlated with this root level cell
 *
 *  handle			The File handle
 *  sfc				The sfc index of root cell
 *  variables			The variables of the root level cell
 *  level			The depth of the Oct tree correlated to the root level cell
 *  num_level_octs		The array store the number of Oct nodes each level
 */
int artio_grid_write_root_cell_begin(artio_fileset *handle, int64_t sfc, 
		float * variables, int level, int * num_octs_per_level);

/*
 * Description:	Do something at the end of writing the root level cell
 */
int artio_grid_write_root_cell_end(artio_fileset *handle);

/*
 * Description:	Do something at the beginning of each level
 */
int artio_grid_write_level_begin(artio_fileset *handle, int level );

/*
 * Description:	Do something at the end of each level
 */
int artio_grid_write_level_end(artio_fileset *handle);

/*
 * Description:	Output the data of a special oct tree node to the file
 *
 *  handle			The handle of the file
 *  variables 			The array recording the variables of the eight cells belonging to this Octree node.
 */
int artio_grid_write_oct(artio_fileset *handle, float *variables, int *refined);

/*
 * Description:	Read the variables of the root level cell and the index of the Octtree
 *              correlated with this root level cell
 *
 *  handle			The File handle
 *  variables			The variables of the root level cell
 *  level 			The depth of the OCT tree
 *  num_octs_per_level		The number of node of each oct level
 *
 */
int artio_grid_read_root_cell_begin(artio_fileset *handle, int64_t sfc, 
		double *pos, float *variables,
		int *num_tree_levels, int *num_octs_per_level);

/*
 * Description:	Do something at the end of reading the root level cell
 */
int artio_grid_read_root_cell_end(artio_fileset *handle);

/*
 * Description:	Do something at the beginning of each level
 */
int artio_grid_read_level_begin(artio_fileset *handle, int level );

/*
 * Description:	Do something at the end of each level
 */
int artio_grid_read_level_end(artio_fileset *handle);

/*
 * Description:	Read the data of a special oct tree node from the file
 */
int artio_grid_read_oct(artio_fileset *handle, double *pos, 
		float *variables, int *refined);

int artio_grid_cache_sfc_range(artio_fileset *handle, int64_t sfc_start, int64_t sfc_end);
int artio_grid_clear_sfc_cache(artio_fileset *handle );

int artio_grid_count_octs_in_sfc_range(artio_fileset *handle,
        int64_t start, int64_t end, int64_t *num_octs_in_range );

/*
 * Description:	Read a segment of oct nodes
 *
 *  handle			file pointer
 *  sfc1			the start sfc index
 *  sfc2			the end sfc index
 *  max_level_to_read		max level to read for each oct tree
 *  option			1. refined nodes; 2 leaf nodes; 3 all nodes
 *  callback        callback function
 *  params          a pointer to user-defined data passed to the callback
 */
int artio_grid_read_sfc_range_levels(artio_fileset *handle, 
		int64_t sfc1, int64_t sfc2, 
		int min_level_to_read, int max_level_to_read, 
		int options, artio_grid_callback callback,
		void *params );

int artio_grid_read_sfc_range(artio_fileset *handle,
        int64_t sfc1, int64_t sfc2, int options,
        artio_grid_callback callback,
		void *params );

int artio_grid_read_selection(artio_fileset *handle,
		artio_selection *selection, int options,
		artio_grid_callback callback,
		void *params );

int artio_grid_read_selection_levels( artio_fileset *handle,
		artio_selection *selection, 
		int min_level_to_read, int max_level_to_read,
		int options,
		artio_grid_callback callback,
		void *params );

/**
 *  header			head file name
 *  num_particle_files		the number of files to record refined cells
 *  allocation_strategy
 *  num_species			number of particle species
 *  species_labels		string identifier for each species
 *  handle			the artio file handle
 *
 */
int artio_fileset_add_particles(artio_fileset *handle, 
        int num_particle_files, int allocation_strategy,
        int num_species, char **species_labels,
        int *num_primary_variables,
        int *num_secondary_variables,
        char ***primary_variable_labels_per_species,
        char ***secondary_variable_labels_per_species,
        int *num_particles_per_species_per_root_tree );

int artio_fileset_open_particles(artio_fileset *handle);
int artio_fileset_close_particles(artio_fileset *handle);

/*
 * Description:	Output the variables of the root level cell and the index of 
 *                  the oct-tree correlated with this root level cell
 *
 *  handle			The File handle
 *  sfc				The sfc index of root cell
 *  variables			The variables of the root level cell
 *  level			The depth of the Oct tree correlated to the root level cell
 *  num_level_octs		The array store the number of Oct nodes each level
 */
int artio_particle_write_root_cell_begin(artio_fileset *handle, int64_t sfc,
		int *num_particles_per_species);

/*
 * Description:	Do something at the end of writing the root level cell
 */
int artio_particle_write_root_cell_end(artio_fileset *handle);

/*
 * Description:	Do something at the beginning of each level
 */
int artio_particle_write_species_begin(artio_fileset *handle, int species );

/*
 * Description:	Do something at the end of each level
 */
int artio_particle_write_species_end(artio_fileset *handle);

/*
 * Description: Output the data of a special oct tree node to the file
 *
 *  handle			The handle of the file
 *  variables 			The array recording the variables of the eight cells belonging to this Octree node.
 */
int artio_particle_write_particle(artio_fileset *handle, int64_t pid, int subspecies, 
			double* primary_variables, float *secondary_variables);

/*
 * Description:	Read the variables of the root level cell and the index of the Octtree
 *              correlated with this root level cell
 *
 *  handle			The File handle
 *  variables			The variables of the root level cell
 *  level 			The depth of the OCT tree
 *  num_octs_per_level		The number of node of each oct level
 *
 */
int artio_particle_read_root_cell_begin(artio_fileset *handle, int64_t sfc, 
			int * num_particle_per_species);

/*
 * Description:	Do something at the end of reading the root level cell
 */
int artio_particle_read_root_cell_end(artio_fileset *handle);

/*
 * Description:	Do something at the beginning of each level
 */
int artio_particle_read_species_begin(artio_fileset *handle, int species );

/*
 * Description:  Do something at the end of each level
 */
int artio_particle_read_species_end(artio_fileset *handle);

/*
 * Description:	Read the data of a single particle from the file
 */
int artio_particle_read_particle(artio_fileset *handle, int64_t *pid, int *subspecies,
			double *primary_variables, float *secondary_variables);

int artio_particle_cache_sfc_range(artio_fileset *handle, int64_t sfc_start, int64_t sfc_end);
int artio_particle_clear_sfc_cache(artio_fileset *handle );                                                          

typedef void (* artio_particle_callback)(int64_t sfc_index,
		int species, int subspecies, int64_t pid, 
		double *primary_variables, float *secondary_variables, void *params );

/*
 * Description: Read a segment of particles
 *
 *  handle			file pointer
 *  sfc1			the start sfc index
 *  sfc2			the end sfc index
 *  start_species   the first particle species to read
 *  end_species     the last particle species to read
 *  callback        callback function
 *  params          user defined data passed to the callback function
 */
int artio_particle_read_sfc_range(artio_fileset *handle, 
		int64_t sfc1, int64_t sfc2, 
		artio_particle_callback callback,
		void *params);

int artio_particle_read_sfc_range_species( artio_fileset *handle, 
        int64_t sfc1, int64_t sfc2, 
        int start_species, int end_species,
        artio_particle_callback callback,
		void *params);

int artio_particle_read_selection(artio_fileset *handle,
        artio_selection *selection, artio_particle_callback callback,
		void *params );

int artio_particle_read_selection_species( artio_fileset *handle,
        artio_selection *selection, int start_species, int end_species,
        artio_particle_callback callback,
		void *params );

artio_selection *artio_selection_allocate( artio_fileset *handle );
artio_selection *artio_select_all( artio_fileset *handle );
artio_selection *artio_select_volume( artio_fileset *handle, double lpos[3], double rpos[3] );
artio_selection *artio_select_cube( artio_fileset *handle, double center[3], double size );
int artio_selection_add_root_cell( artio_selection *selection, int coords[3] );                   
int artio_selection_destroy( artio_selection *selection );
void artio_selection_print( artio_selection *selection );
int artio_selection_iterator( artio_selection *selection,
         int64_t max_range_size, int64_t *start, int64_t *end );
int artio_selection_iterator_reset( artio_selection *selection );
int64_t artio_selection_size( artio_selection *selection );

#endif /* __ARTIO_H__ */
