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

#ifdef ARTIO_MPI
#include "artio_mpi.h"
#endif

artio_file artio_fileset_open(char * file_prefix, int type, artio_context context) ;

/* simpler function for cpdef */
void wrap_artio_fileset_open(char * file_prefix, int type) {
    artio_file junk;
    artio_context context=NULL;
    printf("file_prefix %s type %d\n",file_prefix, type);
    junk = artio_fileset_open(file_prefix, type, context);
}

artio_file artio_fileset_open(char * file_prefix, int type, artio_context context) {
	artio_fh head_fh;
	char filename[256];
	int my_rank, num_procs;
	int re;

	artio_file handle = (artio_file)malloc(sizeof(struct artio_file_struct));
	artio_parameter_list_init(&handle->param_list);

#ifdef ARTIO_MPI
	if(context == NULL)
	  {
	    fprintf(stderr,"Context cannot be NULL in the MPI mode.\n");
	    return NULL;
	  }
	MPI_Comm_size(context->comm, &num_procs);
	MPI_Comm_rank(context->comm, &my_rank);
#else
	num_procs = 1;
	my_rank = 0;
#endif /* MPI */

	handle->open_mode = ARTIO_FILESET_READ;
	handle->open_type = type;

	handle->rank = my_rank;
	handle->num_procs = num_procs;
	handle->context = context;
	strncpy(handle->file_prefix, file_prefix, 255);

	/* open header file */
	snprintf(filename, 255, "%s.art", file_prefix);
	head_fh = artio_file_fopen(filename, 
			ARTIO_MODE_READ | ARTIO_MODE_ACCESS | ARTIO_MODE_DIRECT, context);

	if ( head_fh == NULL ) {
            fprintf(stderr,"snl quitting fileset open: failed to open header\n");
		return NULL;
	}

	re = artio_parameter_read(head_fh, &handle->param_list);
	if ( re != ARTIO_SUCCESS ) {
            fprintf(stderr,"snl quitting fileset open: failed to read header\n");
		return NULL;
	}

	artio_file_fclose(head_fh);

	artio_parameter_get_long(handle, "num_root_cells", &handle->num_root_cells);
        printf("num_root_cells %d\n",handle->num_root_cells); //snl temporary

	/* default to accessing all sfc indices */
	handle->proc_sfc_begin = 0;
	handle->proc_sfc_end = handle->num_root_cells-1;
	handle->proc_sfc_index = NULL;

	/* open data files */
	if (type & ARTIO_OPEN_PARTICLES) {
		re = artio_particle_open(handle);
		if ( re != ARTIO_SUCCESS ) {
			return NULL;
		}
	}

	if (type & ARTIO_OPEN_GRID) {
		re = artio_grid_open(handle);
		if ( re != ARTIO_SUCCESS ) {
			return NULL;
		}
	}

	return handle;
}

artio_file artio_fileset_create(char * file_prefix, int64_t root_cells, 
		int64_t proc_sfc_begin, int64_t proc_sfc_end, artio_context context) {
	int my_rank, num_procs;
	int64_t *proc_sfc_index;
	artio_file handle;

#ifdef ARTIO_MPI
	if(context == NULL)
	  {
	    fprintf(stderr,"Context cannot be NULL in the MPI mode.\n");
	    return NULL;
	  }

	MPI_Comm_size(context->comm, &num_procs);
	MPI_Comm_rank(context->comm, &my_rank);

	if ( proc_sfc_begin < 0 || proc_sfc_end > root_cells ) {
		return NULL;
	}
#else
	num_procs = 1;
	my_rank = 0;

	if ( proc_sfc_begin != 0 || proc_sfc_end != root_cells-1 ) {
		return NULL;
	}
#endif /* MPI */

	handle = (artio_file)malloc(sizeof(struct artio_file_struct));
	artio_parameter_list_init(&handle->param_list);

	handle->open_mode = ARTIO_FILESET_WRITE;
	handle->open_type = 0;

	handle->rank = my_rank;
	handle->num_procs = num_procs;
	handle->context = context;
	handle->num_root_cells = root_cells;

	strncpy(handle->file_prefix, file_prefix, 255);

	proc_sfc_index = (int64_t*)malloc((num_procs+1)*sizeof(int64_t));
	if ( proc_sfc_index == NULL ) {
		fprintf(stderr, "ERROR ALLOCATING MEMORY!\n");
		exit(1);
	}
#ifdef ARTIO_MPI
	MPI_Allgather( &proc_sfc_begin, 1, MPI_LONG_LONG, proc_sfc_index, 1, MPI_LONG_LONG, context->comm );
#else
	proc_sfc_index[0] = 0;
#endif /* ARTIO_MPI */
	proc_sfc_index[handle->num_procs] = root_cells;
	handle->proc_sfc_index = proc_sfc_index;
	handle->proc_sfc_begin = proc_sfc_begin;
	handle->proc_sfc_end = proc_sfc_end;

	artio_parameter_set_long(handle, "num_root_cells", root_cells);

	return handle;
}

int artio_fileset_close(artio_file handle) {
	char header_filename[256];
	artio_fh head_fh;

	if (handle->open_mode == ARTIO_FILESET_WRITE) {
		snprintf(header_filename, 255, "%s.art", handle->file_prefix);
		head_fh = artio_file_fopen(header_filename, 
				ARTIO_MODE_WRITE | ARTIO_MODE_DIRECT |
					   ((handle->rank == 0) ? ARTIO_MODE_ACCESS : 0), handle->context);

		if (0 == handle->rank) {
			artio_parameter_write(head_fh, &handle->param_list);
		}

		artio_file_fclose(head_fh);
		free(handle->proc_sfc_index);
	}

	if (handle->open_type & ARTIO_OPEN_GRID) {
		artio_grid_close(handle);
		free(handle->grid);
	}

	if (handle->open_type & ARTIO_OPEN_PARTICLES) {
		artio_particle_close(handle);
		free(handle->particle);
	}

	artio_parameter_free_list(&handle->param_list);
	free(handle);

	return ARTIO_SUCCESS;
}
