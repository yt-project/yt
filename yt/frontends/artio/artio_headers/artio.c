/*
 * artio.c
 *
 *  Created on: Feb 21, 2010
 *  Author: Yongen Yu
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "artio.h"
#include "artio_internal.h"

artio_fileset *artio_fileset_allocate( char *file_prefix, int mode,
		const artio_context *context );
void artio_fileset_destroy( artio_fileset *handle );

artio_fileset *artio_fileset_open(char * file_prefix, int type, const artio_context *context) {
	artio_fh *head_fh;
	char filename[256];
	int ret;
	int64_t tmp;

	artio_fileset *handle = 
		artio_fileset_allocate( file_prefix, ARTIO_FILESET_READ, context );
	if ( handle == NULL ) {
		return NULL;
	}

	/* open header file */
	sprintf(filename, "%s.art", handle->file_prefix);
	head_fh = artio_file_fopen(filename, 
			ARTIO_MODE_READ | ARTIO_MODE_ACCESS | ARTIO_MODE_DIRECT, context);

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
				ARTIO_MODE_WRITE | ARTIO_MODE_DIRECT |
					   ((handle->rank == 0) ? ARTIO_MODE_ACCESS : 0), 
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

		handle->context = (artio_context *)malloc(sizeof(artio_context));
		if ( handle->context == NULL ) {
			return NULL;
		}
		memcpy( handle->context, context, sizeof(artio_context) );

#ifdef ARTIO_MPI
		MPI_Comm_size(handle->context->comm, &num_procs);
		MPI_Comm_rank(handle->context->comm, &my_rank);
#else
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
