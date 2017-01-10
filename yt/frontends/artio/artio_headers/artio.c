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
#include "artio.h"
#include "artio_internal.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _WIN32
typedef __int64 int64_t;
typedef __int32 int32_t;
#else
#include <stdint.h>
#endif

artio_fileset *artio_fileset_allocate( char *file_prefix, int mode,
		const artio_context *context );
void artio_fileset_destroy( artio_fileset *handle );

int artio_fh_buffer_size = ARTIO_DEFAULT_BUFFER_SIZE;

int artio_fileset_set_buffer_size( int buffer_size ) {
	if ( buffer_size < 0 ) {
		return ARTIO_ERR_INVALID_BUFFER_SIZE;
	}

	artio_fh_buffer_size = buffer_size;
	return ARTIO_SUCCESS;
}

artio_fileset *artio_fileset_open(char * file_prefix, int type, const artio_context *context) {
	artio_fh *head_fh;
	char filename[256];
	int ret;
	int64_t tmp;
	int artio_major, artio_minor;

	artio_fileset *handle = 
		artio_fileset_allocate( file_prefix, ARTIO_FILESET_READ, context );
	if ( handle == NULL ) {
		return NULL;
	}

	/* open header file */
	sprintf(filename, "%s.art", handle->file_prefix);
	head_fh = artio_file_fopen(filename, 
			ARTIO_MODE_READ | ARTIO_MODE_ACCESS, context);

	if ( head_fh == NULL ) {
		artio_fileset_destroy(handle);
		return NULL;
	}

	ret = artio_parameter_read(head_fh, handle->parameters );
	if ( ret != ARTIO_SUCCESS ) {
		artio_fileset_destroy(handle);
		return NULL;
	}

	artio_file_fclose(head_fh);

	/* check versions */
	if ( artio_parameter_get_int(handle, "ARTIO_MAJOR_VERSION", &artio_major ) == 
			ARTIO_ERR_PARAM_NOT_FOUND ) {
		/* version pre 1.0 */
		artio_major = 0;
		artio_minor = 9;
	} else {
		artio_parameter_get_int(handle, "ARTIO_MINOR_VERSION", &artio_minor );
	}

	if ( artio_major > ARTIO_MAJOR_VERSION ) {
		fprintf(stderr,"ERROR: artio file version newer than library (%u.%u vs %u.%u).\n",
			artio_major, artio_minor, ARTIO_MAJOR_VERSION, ARTIO_MINOR_VERSION );
		artio_fileset_destroy(handle);
		return NULL;
	}
	
	artio_parameter_get_long(handle, "num_root_cells", &handle->num_root_cells);
	
	if ( artio_parameter_get_int(handle, "sfc_type", &handle->sfc_type ) != ARTIO_SUCCESS ) {
		handle->sfc_type = ARTIO_SFC_HILBERT;
	}

	handle->nBitsPerDim = 0;
	tmp = handle->num_root_cells >> 3;
	while ( tmp ) {
		handle->nBitsPerDim++;
		tmp >>= 3;
	}
	handle->num_grid = 1<<handle->nBitsPerDim;

	/* default to accessing all sfc indices */
	handle->proc_sfc_begin = 0;
	handle->proc_sfc_end = handle->num_root_cells-1;

	/* open data files */
	if (type & ARTIO_OPEN_PARTICLES) {
		ret = artio_fileset_open_particles(handle);
		if ( ret != ARTIO_SUCCESS ) {
			artio_fileset_destroy(handle);
			return NULL;
		}
	}

	if (type & ARTIO_OPEN_GRID) {
		ret = artio_fileset_open_grid(handle);
		if ( ret != ARTIO_SUCCESS ) {
			artio_fileset_destroy(handle);
			return NULL;
		}
	}

	return handle;
}

artio_fileset *artio_fileset_create(char * file_prefix, int64_t root_cells, 
		int64_t proc_sfc_begin, int64_t proc_sfc_end, const artio_context *context) {
    artio_fileset *handle = 
		artio_fileset_allocate( file_prefix, ARTIO_FILESET_WRITE, context );
    if ( handle == NULL ) {
        return NULL;
    }

	handle->proc_sfc_index = 
		(int64_t*)malloc((handle->num_procs+1)*sizeof(int64_t));
	if ( handle->proc_sfc_index == NULL ) {
		artio_fileset_destroy(handle);
		return NULL;
	}

#ifdef ARTIO_MPI
	MPI_Allgather( &proc_sfc_begin, 1, MPI_LONG_LONG, 
			handle->proc_sfc_index, 1, MPI_LONG_LONG, handle->context->comm );
#else
	handle->proc_sfc_index[0] = 0;
#endif /* ARTIO_MPI */
	handle->proc_sfc_index[handle->num_procs] = root_cells;

	handle->proc_sfc_begin = proc_sfc_begin;
	handle->proc_sfc_end = proc_sfc_end;
	handle->num_root_cells = root_cells;

	artio_parameter_set_long(handle, "num_root_cells", root_cells);

	artio_parameter_set_int(handle, "ARTIO_MAJOR_VERSION", ARTIO_MAJOR_VERSION );
	artio_parameter_set_int(handle, "ARTIO_MINOR_VERSION", ARTIO_MINOR_VERSION );

	return handle;
}

int artio_fileset_close(artio_fileset *handle) {
	char header_filename[256];
	artio_fh *head_fh;
	
	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if (handle->open_mode == ARTIO_FILESET_WRITE) {
		/* ensure we've flushed open particle and 
		 * grid files before writing header */
		if ( handle->grid != NULL ) {
			artio_fileset_close_grid(handle);
		}

		if ( handle->particle != NULL ) {
			artio_fileset_close_particles(handle);
		}

		sprintf(header_filename, "%s.art", handle->file_prefix);
		head_fh = artio_file_fopen(header_filename, 
				ARTIO_MODE_WRITE | ((handle->rank == 0) ? ARTIO_MODE_ACCESS : 0), 
				handle->context);

		if (head_fh == NULL) {
			return ARTIO_ERR_FILE_CREATE;
		}

		if (0 == handle->rank) {
			artio_parameter_write(head_fh, handle->parameters );
		}

		artio_file_fclose(head_fh);
	}

	artio_fileset_destroy(handle);

	return ARTIO_SUCCESS;
}

artio_fileset *artio_fileset_allocate( char *file_prefix, int mode, 
		const artio_context *context ) {
	int my_rank;
	int num_procs;

    artio_fileset *handle = (artio_fileset *)malloc(sizeof(artio_fileset));
	if ( handle != NULL ) {
		handle->parameters = artio_parameter_list_init();

#ifdef ARTIO_MPI
		handle->context = (artio_context *)malloc(sizeof(artio_context));
		if ( handle->context == NULL ) {
			return NULL;
		}
		memcpy( handle->context, context, sizeof(artio_context) );

		MPI_Comm_size(handle->context->comm, &num_procs);
		MPI_Comm_rank(handle->context->comm, &my_rank);
#else
		handle->context = NULL;

		num_procs = 1;
		my_rank = 0;
#endif /* MPI */

		strncpy(handle->file_prefix, file_prefix, 250);

		handle->open_mode = mode;
		handle->open_type = ARTIO_OPEN_HEADER;

		handle->rank = my_rank;
		handle->num_procs = num_procs;
		handle->endian_swap = 0;

		handle->proc_sfc_index = NULL;
		handle->proc_sfc_begin = -1;
		handle->proc_sfc_end = -1;
		handle->num_root_cells = -1;
		
		handle->grid = NULL;
		handle->particle = NULL;
	}	
	return handle;
}

void artio_fileset_destroy( artio_fileset *handle ) {
	if ( handle == NULL ) return;

	if ( handle->proc_sfc_index != NULL ) free( handle->proc_sfc_index );
	
	if ( handle->grid != NULL ) {
        artio_fileset_close_grid(handle);
    }

    if ( handle->particle != NULL ) {
        artio_fileset_close_particles(handle);
    }

	if ( handle->context != NULL ) free( handle->context );

	artio_parameter_list_free(handle->parameters);

	free(handle);
}

int artio_fileset_has_grid( artio_fileset *handle ) {
	int num_grid_files = 0;
	return ( handle->grid != NULL ||
		( artio_parameter_get_int( handle, "num_grid_files", &num_grid_files ) == ARTIO_SUCCESS &&
		  num_grid_files > 0 ) );
}

int artio_fileset_has_particles( artio_fileset *handle ) {
	int num_particle_files = 0;
	return ( handle->particle != NULL ||
			( artio_parameter_get_int( handle, "num_particle_files", &num_particle_files ) == ARTIO_SUCCESS &&
			  num_particle_files > 0 ) );
}
