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

#ifndef __ARTIO_INTERNAL_H__
#define __ARTIO_INTERNAL_H__

#ifdef ARTIO_MPI
#include <mpi.h>
#endif

#include <stdlib.h>
#include <limits.h>

#ifdef _WIN32
typedef __int64 int64_t;
typedef __int32 int32_t;
#else
#include <stdint.h>
#endif

#include "artio_endian.h"

#ifndef ARTIO_DEFAULT_BUFFER_SIZE
#define ARTIO_DEFAULT_BUFFER_SIZE	65536
#endif

extern int artio_fh_buffer_size;

#define nDim			3

#ifndef MIN
#define MIN(x,y)        (((x) < (y)) ? (x): (y))
#endif
#ifndef MAX
#define MAX(x,y)        (((x) > (y)) ? (x): (y))
#endif

/* limit individual writes to 32-bit safe quantity */
#define ARTIO_IO_MAX	(1<<30)

#ifdef INT64_MAX
#define ARTIO_INT64_MAX	INT64_MAX
#else
#define ARTIO_INT64_MAX 0x7fffffffffffffffLL
#endif

typedef struct ARTIO_FH artio_fh;

typedef struct artio_particle_file_struct {
	artio_fh **ffh;

	void *buffer;
	int buffer_size;

	int num_particle_files;
	int64_t *file_sfc_index;
	int64_t cache_sfc_begin;
	int64_t cache_sfc_end;
	int64_t *sfc_offset_table;

	/* maintained for consistency and user-error detection */
	int num_species;
	int cur_file;
	int cur_species;
	int cur_particle;
	int64_t cur_sfc;
	int *num_primary_variables;
	int *num_secondary_variables;
	int *num_particles_per_species;
} artio_particle_file;

typedef struct artio_grid_file_struct {
	artio_fh **ffh;

	void *buffer;
	int buffer_size;

	int num_grid_variables;
	int num_grid_files;
	int64_t *file_sfc_index;
	int64_t cache_sfc_begin;
	int64_t cache_sfc_end;
	int64_t *sfc_offset_table;

	int file_max_level;
	/* maintained for consistency and user-error detection */
	int cur_file;
	int cur_num_levels;
	int cur_level;
	int cur_octs;
	int64_t cur_sfc;
	int *octs_per_level;

	int pos_flag;
	int pos_cur_level;
	int next_level_size;
	int cur_level_size;
	double cell_size_level;
	double *next_level_pos;
	double *cur_level_pos;
	int next_level_oct;
	
} artio_grid_file;

typedef struct parameter_struct {
	int key_length;
	char key[64];
	int val_length;
	int type;
	char *value;
	struct parameter_struct * next;
} parameter;

typedef struct parameter_list_struct {
	parameter * head;
	parameter * tail;
	parameter * cursor;
	int iterate_flag;
} parameter_list;

struct artio_fileset_struct {
	char file_prefix[256];
	int endian_swap;
	int open_type;
	int open_mode;
	int rank;
	int num_procs;
	artio_context *context;

	int64_t *proc_sfc_index;
	int64_t proc_sfc_begin;
	int64_t proc_sfc_end;
	int64_t num_root_cells;
	int sfc_type;
	int nBitsPerDim;
	int num_grid;
	
	parameter_list *parameters;
	artio_grid_file *grid;
	artio_particle_file *particle;
};

struct artio_selection_struct {
    int64_t *list;
    int size;
    int num_ranges; 
	int cursor;
	int64_t subcycle;
	artio_fileset *fileset;
};

#define ARTIO_FILESET_READ      0
#define ARTIO_FILESET_WRITE     1

#define ARTIO_MODE_READ         1
#define ARTIO_MODE_WRITE        2
#define ARTIO_MODE_ACCESS       4
#define ARTIO_MODE_ENDIAN_SWAP  8

#define ARTIO_SEEK_SET          0
#define ARTIO_SEEK_CUR          1
#define ARTIO_SEEK_END			2

/* wrapper functions for profiling and debugging */
artio_fh *artio_file_fopen( char * filename, int amode, const artio_context *context );
int artio_file_attach_buffer( artio_fh *handle, void *buf, int buf_size );
int artio_file_detach_buffer( artio_fh *handle );
int artio_file_fwrite(artio_fh *handle, const void *buf, int64_t count, int type );
int artio_file_ftell( artio_fh *handle, int64_t *offset );
int artio_file_fflush(artio_fh *handle);
int artio_file_fseek(artio_fh *ffh, int64_t offset, int whence);
int artio_file_fread(artio_fh *handle, void *buf, int64_t count, int type );
int artio_file_fclose(artio_fh *handle);
void artio_file_set_endian_swap_tag(artio_fh *handle);

/* internal versions */
artio_fh *artio_file_fopen_i( char * filename, int amode, const artio_context *context );
int artio_file_attach_buffer_i( artio_fh *handle, void *buf, int buf_size );
int artio_file_detach_buffer_i( artio_fh *handle );
int artio_file_fwrite_i(artio_fh *handle, const void *buf, int64_t count, int type );
int artio_file_ftell_i( artio_fh *handle, int64_t *offset );
int artio_file_fflush_i(artio_fh *handle);
int artio_file_fseek_i(artio_fh *ffh, int64_t offset, int whence);
int artio_file_fread_i(artio_fh *handle, void *buf, int64_t count, int type );
int artio_file_fclose_i(artio_fh *handle);
void artio_file_set_endian_swap_tag_i(artio_fh *handle);

#define ARTIO_ENDIAN_MAGIC	0x1234

parameter_list *artio_parameter_list_init(void);

parameter *artio_parameter_list_search(parameter_list *parameters, const char *key);

int artio_parameter_array_length( parameter *item );

int artio_parameter_list_insert(parameter_list *parameters, const char *key, int length,
		void * value, int type);

int artio_parameter_read(artio_fh *handle, parameter_list *parameters);

int artio_parameter_write(artio_fh *handle, parameter_list *parameters);

int artio_parameter_list_print(parameter_list *parameters);
int artio_parameter_list_free(parameter_list *parameters);
int artio_parameter_list_print(parameter_list *parameters);

size_t artio_type_size( int type );
int artio_parameter_list_unpack(parameter_list *parameters, const char *key, int length, void *value, int type );

#define ARTIO_SFC_SLAB_X	0
#define ARTIO_SFC_MORTION	1
#define ARTIO_SFC_HILBERT	2
#define ARTIO_SFC_SLAB_Y	3
#define ARTIO_SFC_SLAB_Z	4

int64_t artio_sfc_index_position( artio_fileset *handle, double position[nDim] );
int64_t artio_sfc_index( artio_fileset *handle, int coords[nDim] );
void artio_sfc_coords( artio_fileset *handle, int64_t index, int coords[nDim] );

#endif /* __ARTIO_INTERNAL_H__ */
