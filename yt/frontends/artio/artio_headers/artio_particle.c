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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _WIN32
typedef __int64 int64_t;
typedef __int32 int32_t;
#else
#include <stdint.h>
#endif

int artio_particle_find_file(artio_particle_file *phandle, int start, int end, int64_t sfc);
artio_particle_file *artio_particle_file_allocate(void);
void artio_particle_file_destroy( artio_particle_file *phandle );

/*
 * Open existing particle files and add to fileset
 */
int artio_fileset_open_particles(artio_fileset *handle) {
	int i;
	char filename[256];
	int first_file, last_file;
	int mode;
	artio_particle_file *phandle;

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if ( handle->open_type & ARTIO_OPEN_PARTICLES ||
			handle->open_mode != ARTIO_FILESET_READ ||
			handle->particle != NULL ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}
	handle->open_type |= ARTIO_OPEN_PARTICLES;

	phandle = artio_particle_file_allocate();
	if ( phandle == NULL ) {
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}
	
	if ( artio_parameter_get_int( handle, "num_particle_files", &phandle->num_particle_files ) != ARTIO_SUCCESS ||
			artio_parameter_get_int(handle, "num_particle_species", &phandle->num_species) != ARTIO_SUCCESS ) {
		return ARTIO_ERR_PARTICLE_DATA_NOT_FOUND;
	}

	phandle->num_primary_variables = (int *)malloc(sizeof(int) * phandle->num_species);
	if ( phandle->num_primary_variables == NULL ) {
		artio_particle_file_destroy(phandle);
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}

	phandle->num_secondary_variables = (int *)malloc(sizeof(int) * phandle->num_species);
	if ( phandle->num_primary_variables == NULL ) {
		artio_particle_file_destroy(phandle);
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}

	phandle->num_particles_per_species =
		(int *)malloc(phandle->num_species * sizeof(int));
	if ( phandle->num_particles_per_species == NULL ) {
		artio_particle_file_destroy(phandle);
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}

	artio_parameter_get_int_array(handle, "num_primary_variables",
			phandle->num_species, phandle->num_primary_variables);
	artio_parameter_get_int_array(handle, "num_secondary_variables",
			phandle->num_species, phandle->num_secondary_variables);

	phandle->file_sfc_index = (int64_t *)malloc(sizeof(int64_t) * (phandle->num_particle_files + 1));
	if ( phandle->file_sfc_index == NULL ) {
		artio_particle_file_destroy(phandle);
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}

	artio_parameter_get_long_array(handle, "particle_file_sfc_index",
			phandle->num_particle_files + 1, phandle->file_sfc_index);

	first_file = artio_particle_find_file(phandle, 0, 
			phandle->num_particle_files, handle->proc_sfc_begin);
	last_file = artio_particle_find_file(phandle, first_file,
			phandle->num_particle_files, handle->proc_sfc_end);

	/* allocate file handles */
	phandle->ffh = (artio_fh **)malloc(phandle->num_particle_files * sizeof(artio_fh *));
	if ( phandle->ffh == NULL ) {
		artio_particle_file_destroy(phandle);
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}

	for ( i = 0; i < phandle->num_particle_files; i++ ) {
		phandle->ffh[i] = NULL;
	}

	/* open files on all processes */
	for (i = 0; i < phandle->num_particle_files; i++) {
		sprintf(filename, "%s.p%03d", handle->file_prefix, i);

		mode = ARTIO_MODE_READ;
		if (i >= first_file && i <= last_file) {
			mode |= ARTIO_MODE_ACCESS;
		}
		if (handle->endian_swap) {
			mode |= ARTIO_MODE_ENDIAN_SWAP;
		}

		phandle->ffh[i] = artio_file_fopen(filename, mode, handle->context);
		if ( phandle->ffh[i] == NULL ) {
			artio_particle_file_destroy(phandle);
			return ARTIO_ERR_PARTICLE_FILE_NOT_FOUND;
		}
	}

	handle->particle = phandle;
	return ARTIO_SUCCESS;
}

int artio_fileset_add_particles( artio_fileset *handle,
        int num_particle_files, int allocation_strategy,
        int num_species, char ** species_labels,
        int * num_primary_variables,
        int * num_secondary_variables,
        char *** primary_variable_labels_per_species,
        char *** secondary_variable_labels_per_species,
		int * num_particles_per_species_per_root_tree ) {
	int i, k;
	int ret;
	int64_t l, cur;
	int64_t first_file_sfc, last_file_sfc;

	int first_file, last_file;
	char filename[256];
	char species_label[64];
	int mode;
	int64_t *local_particles_per_species, *total_particles_per_species;
	artio_particle_file *phandle;

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if ( handle->open_mode != ARTIO_FILESET_WRITE ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	if ( handle->open_type & ARTIO_OPEN_PARTICLES ||
			handle->particle != NULL ) {
		return ARTIO_ERR_DATA_EXISTS;
	}
	handle->open_type |= ARTIO_OPEN_PARTICLES;

	/* compute total number of particles per species */
	local_particles_per_species = (int64_t *)malloc( num_species * sizeof(int64_t));
	if ( local_particles_per_species == NULL ) {
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}

	total_particles_per_species = (int64_t *)malloc( num_species * sizeof(int64_t));
	if ( total_particles_per_species == NULL ) {
		free( local_particles_per_species );
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}

	for ( i = 0; i < num_species; i++ ) {
		local_particles_per_species[i] = 0;
	}

	for ( l = 0; l < handle->proc_sfc_end-handle->proc_sfc_begin+1; l++ ) {
		for ( i = 0; i < num_species; i++ ) {
			local_particles_per_species[i] += 
					num_particles_per_species_per_root_tree[num_species*l+i];
		}
	}

#ifdef ARTIO_MPI
	MPI_Allreduce( local_particles_per_species, total_particles_per_species, 
			num_species, MPI_LONG_LONG_INT, MPI_SUM, handle->context->comm );
#else
	for ( i = 0; i < num_species; i++ ) {
		total_particles_per_species[i] = local_particles_per_species[i];
	}
#endif

	artio_parameter_set_long_array(handle, 
			"num_particles_per_species", num_species,
			total_particles_per_species );

	free( local_particles_per_species );
	free( total_particles_per_species );

	artio_parameter_set_int(handle, "num_particle_files", num_particle_files);
	artio_parameter_set_int(handle, "num_particle_species", num_species);
	artio_parameter_set_string_array(handle, "particle_species_labels",
			num_species, species_labels);
	artio_parameter_set_int_array(handle, "num_primary_variables", num_species,
			num_primary_variables);
	artio_parameter_set_int_array(handle, "num_secondary_variables",
			num_species, num_secondary_variables);

	for(i=0;i<num_species;i++) {
		sprintf( species_label, "species_%02u_primary_variable_labels", i );
		artio_parameter_set_string_array( handle, species_label,
				num_primary_variables[i], primary_variable_labels_per_species[i] );

		sprintf( species_label, "species_%02u_secondary_variable_labels", i );
		artio_parameter_set_string_array( handle, species_label,
				num_secondary_variables[i], secondary_variable_labels_per_species[i] );
	}

	phandle = artio_particle_file_allocate();
	if ( phandle == NULL ) {
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}

	phandle->file_sfc_index = (int64_t *)malloc(sizeof(int64_t) * (num_particle_files + 1));
	if ( phandle->file_sfc_index == NULL ) {
		artio_particle_file_destroy(phandle);
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}

	switch (allocation_strategy) {
		case ARTIO_ALLOC_EQUAL_PROC:
			if (num_particle_files > handle->num_procs) {
				return ARTIO_ERR_INVALID_FILE_NUMBER;
			}

			for (i = 0; i < num_particle_files; i++) {
				phandle->file_sfc_index[i] = handle->proc_sfc_index[ 
					(handle->num_procs*i+num_particle_files-1) / num_particle_files ];
			}
			phandle->file_sfc_index[num_particle_files] = handle->proc_sfc_index[handle->num_procs];
			break;
		case ARTIO_ALLOC_EQUAL_SFC:
			for (i = 0; i < num_particle_files; i++) {
				phandle->file_sfc_index[i] = (handle->num_root_cells*i+num_particle_files-1) / 
													num_particle_files;
			}
			phandle->file_sfc_index[num_particle_files] = handle->num_root_cells;
			break;
		default:
			artio_particle_file_destroy(phandle);
			return ARTIO_ERR_INVALID_ALLOC_STRATEGY;
	}

	phandle->num_particle_files = num_particle_files;
	phandle->num_species = num_species;

	phandle->num_primary_variables = (int *)malloc(sizeof(int) * num_species);
	if ( phandle->num_primary_variables == NULL ) {
		artio_particle_file_destroy(phandle);
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}

	phandle->num_secondary_variables = (int *)malloc(sizeof(int) * num_species);
	if ( phandle->num_secondary_variables == NULL ) {
		artio_particle_file_destroy(phandle);
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}

	phandle->num_particles_per_species = 
		(int *)malloc(phandle->num_species * sizeof(int));
	if ( phandle->num_particles_per_species == NULL ) {
		artio_particle_file_destroy(phandle);
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}

	for (i = 0; i < num_species; i++) {
		phandle->num_primary_variables[i] = num_primary_variables[i];
		phandle->num_secondary_variables[i] = num_secondary_variables[i];
	}

	/* allocate space for sfc offset cache */
	phandle->cache_sfc_begin = handle->proc_sfc_begin;
	phandle->cache_sfc_end = handle->proc_sfc_end;
	phandle->sfc_offset_table = 
		(int64_t *)malloc( (size_t)(handle->proc_sfc_end - handle->proc_sfc_begin + 1) * sizeof(int64_t));
	if ( phandle->sfc_offset_table == NULL ) {
		artio_particle_file_destroy(phandle);
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}

	/* allocate file handles */
	phandle->ffh = (artio_fh **)malloc(num_particle_files * sizeof(artio_fh *));
	if ( phandle->ffh == NULL ) {
		artio_particle_file_destroy(phandle);
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}

	for ( i = 0; i < num_particle_files; i++ ) {
		phandle->ffh[i] = NULL;
	}

	/* open file handles */
	first_file = artio_particle_find_file(phandle, 0, 
					num_particle_files, handle->proc_sfc_begin);
	last_file = artio_particle_find_file(phandle, first_file, 
					num_particle_files, handle->proc_sfc_end);

	for (i = 0; i < num_particle_files; i++) {
		sprintf(filename, "%s.p%03d", handle->file_prefix, i);
		mode = ARTIO_MODE_WRITE;
		if (i >= first_file && i <= last_file) {
			mode |= ARTIO_MODE_ACCESS;
		}

		phandle->ffh[i] = artio_file_fopen(filename, mode, handle->context);
		if ( phandle->ffh[i] == NULL ) {
			artio_particle_file_destroy(phandle);
			return ARTIO_ERR_FILE_CREATE;
		}

		/* write sfc offset header if we contribute to this file */
		if (i >= first_file && i <= last_file) {
#ifdef ARTIO_MPI
			if ( phandle->file_sfc_index[i] >= handle->proc_sfc_index[ handle->rank ] &&
					phandle->file_sfc_index[i] < handle->proc_sfc_index[ handle->rank + 1 ] ) {
				cur = (phandle->file_sfc_index[i+1] - phandle->file_sfc_index[i]) * sizeof(int64_t);
			} else {
				/* obtain offset from previous process */
				MPI_Recv( &cur, 1, MPI_LONG_LONG_INT, handle->rank - 1, i,
						handle->context->comm, MPI_STATUS_IGNORE );
			}
#else
			cur = (phandle->file_sfc_index[i+1] - phandle->file_sfc_index[i]) * sizeof(int64_t);
#endif /* ARTIO_MPI */

			first_file_sfc = MAX( phandle->cache_sfc_begin, phandle->file_sfc_index[i] );
			last_file_sfc = MIN( phandle->cache_sfc_end, phandle->file_sfc_index[i+1]-1 );

			for (l = first_file_sfc - phandle->cache_sfc_begin; 
					l < last_file_sfc - phandle->cache_sfc_begin + 1; l++) {

				phandle->sfc_offset_table[l] = cur;
				cur += sizeof(int) * num_species;

				for (k = 0; k < num_species; k++) {
					cur += num_particles_per_species_per_root_tree[l*num_species+k]
							* (sizeof(int64_t) + sizeof(int) +
									num_primary_variables[k] * sizeof(double) + 
									num_secondary_variables[k] * sizeof(float));
				}
			}

#ifdef ARTIO_MPI
			if ( phandle->file_sfc_index[i+1] > handle->proc_sfc_end+1 ) {
				MPI_Send( &cur, 1, MPI_LONG_LONG_INT, handle->rank+1, i, handle->context->comm );
			}
#endif /* ARTIO_MPI */

			/* seek and write our portion of sfc table */
			ret = artio_file_fseek(phandle->ffh[i], 
					(first_file_sfc - phandle->file_sfc_index[i]) * sizeof(int64_t), ARTIO_SEEK_SET);
			if ( ret != ARTIO_SUCCESS ) {
				artio_particle_file_destroy(phandle);
				return ret;
			}

			ret = artio_file_fwrite(phandle->ffh[i],
					&phandle->sfc_offset_table[first_file_sfc - phandle->cache_sfc_begin],
					last_file_sfc - first_file_sfc + 1, ARTIO_TYPE_LONG);
			if ( ret != ARTIO_SUCCESS ) {
				artio_particle_file_destroy(phandle);
				return ret;
			}
		}
	}

	artio_parameter_set_long_array(handle, "particle_file_sfc_index",
			phandle->num_particle_files + 1, 
			phandle->file_sfc_index);

	handle->particle = phandle;
	return ARTIO_SUCCESS;
}

artio_particle_file *artio_particle_file_allocate(void) {
	artio_particle_file *phandle = 
		(artio_particle_file *)malloc(sizeof(artio_particle_file));
	if ( phandle != NULL ) {
		phandle->ffh = NULL;
		phandle->num_particle_files = -1;
		phandle->file_sfc_index = NULL;
		phandle->cache_sfc_begin = -1;
		phandle->cache_sfc_end = -1;
		phandle->sfc_offset_table = NULL;
		phandle->num_species = -1;
		phandle->cur_particle = -1;
		phandle->cur_sfc = -1;
		phandle->num_primary_variables = NULL;
		phandle->num_secondary_variables = NULL;
		phandle->num_particles_per_species = NULL;
		phandle->cur_file = -1;
		phandle->buffer_size = artio_fh_buffer_size;
		phandle->buffer = malloc(phandle->buffer_size);
		if ( phandle->buffer == NULL ) {
			free(phandle);
			return NULL;
		}	
	}
	return phandle;
}

void artio_particle_file_destroy( artio_particle_file *phandle ) {
	int i;

	if ( phandle == NULL ) return;
	
	if ( phandle->ffh != NULL ) {
	    for (i = 0; i < phandle->num_particle_files; i++) {
			if ( phandle->ffh[i] != NULL ) {
		        artio_file_fclose(phandle->ffh[i]);
			}
    	}
	    free(phandle->ffh);
	}

    if (phandle->sfc_offset_table != NULL) free(phandle->sfc_offset_table);
    if (phandle->num_particles_per_species != NULL) free(phandle->num_particles_per_species);
	if (phandle->num_primary_variables != NULL) free(phandle->num_primary_variables);
	if (phandle->num_secondary_variables != NULL) free(phandle->num_secondary_variables);
	if (phandle->file_sfc_index != NULL) free(phandle->file_sfc_index);
	if (phandle->buffer != NULL) free(phandle->buffer);

	free(phandle);
}

int artio_fileset_close_particles(artio_fileset *handle) {
	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if ( !(handle->open_type & ARTIO_OPEN_PARTICLES) ||
            handle->particle == NULL ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}
	artio_particle_file_destroy(handle->particle);
	handle->particle = NULL;
	return ARTIO_SUCCESS;
}

int artio_particle_cache_sfc_range(artio_fileset *handle, 
		int64_t start, int64_t end) {
	int i;
	int ret;
	int first_file, last_file;
	int64_t min, count, cur;
	artio_particle_file *phandle;

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if (handle->open_mode != ARTIO_FILESET_READ || 
			!(handle->open_type & ARTIO_OPEN_PARTICLES) ||
			handle->particle == NULL ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	phandle = handle->particle;

	if ( start > end || start < handle->proc_sfc_begin || 
			end > handle->proc_sfc_end) {
		return ARTIO_ERR_INVALID_SFC_RANGE;
	}

	/* check if we've already cached the range */
	if ( start >= phandle->cache_sfc_begin &&
			end <= phandle->cache_sfc_end ) {
		return ARTIO_SUCCESS;
	}

	artio_grid_clear_sfc_cache(handle);

	first_file = artio_particle_find_file(phandle, 0, 
			phandle->num_particle_files, start);
	last_file = artio_particle_find_file(phandle, first_file,
			phandle->num_particle_files, end);

	phandle->cache_sfc_begin = start;
	phandle->cache_sfc_end = end;
	phandle->sfc_offset_table = (int64_t *)malloc(sizeof(int64_t) * (size_t)(end - start + 1));
	if ( phandle->sfc_offset_table == NULL ) {
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}

	if ( phandle->cur_file != -1 ) {
		artio_file_detach_buffer( phandle->ffh[phandle->cur_file]);
		phandle->cur_file = -1;
	}

	cur = 0;
	for (i = first_file; i <= last_file; i++) {
		min = MAX( 0, start - phandle->file_sfc_index[i] );
		count = MIN( phandle->file_sfc_index[i+1], end+1 )
				- MAX( start, phandle->file_sfc_index[i] );

		artio_file_attach_buffer( phandle->ffh[i], 
				phandle->buffer, phandle->buffer_size );

		ret = artio_file_fseek(phandle->ffh[i], 
				sizeof(int64_t) * min,
				ARTIO_SEEK_SET);
		if ( ret != ARTIO_SUCCESS ) return ret;

		ret = artio_file_fread(phandle->ffh[i], 
				&phandle->sfc_offset_table[cur],
				count, ARTIO_TYPE_LONG);
		if ( ret != ARTIO_SUCCESS ) return ret;

		artio_file_detach_buffer( phandle->ffh[i] );
		cur += count;
	}

	return ARTIO_SUCCESS;
}

int artio_particle_clear_sfc_cache( artio_fileset *handle ) {
	artio_particle_file *phandle;

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if (handle->open_mode != ARTIO_FILESET_READ ||
			!(handle->open_type & ARTIO_OPEN_GRID) ||
			handle->grid == NULL ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	phandle = handle->particle;

	if ( phandle->sfc_offset_table != NULL ) {
		free(phandle->sfc_offset_table);
		phandle->sfc_offset_table = NULL;
	}

	phandle->cache_sfc_begin = -1;
	phandle->cache_sfc_end = -1;

	return ARTIO_SUCCESS;
}


int artio_particle_find_file(artio_particle_file *phandle,
		int start, int end, int64_t sfc) {
	int j;

	if ( start < 0 || start > phandle->num_particle_files ||
			end < 0 || end > phandle->num_particle_files || 
			sfc < phandle->file_sfc_index[start] ||
            sfc >= phandle->file_sfc_index[end] ) {
		return -1;
	}

	if (start == end || 
			sfc == phandle->file_sfc_index[start] ) {
		return start;
	}

	if (1 == end - start) {
		if (sfc < phandle->file_sfc_index[end]) {
			return start;
		} else {
			return end;
		}
	}

	j = start + (end - start) / 2;
	if (sfc > phandle->file_sfc_index[j]) {
		return artio_particle_find_file(phandle, j, end, sfc);
	} else if (sfc < phandle->file_sfc_index[j]) {
		return artio_particle_find_file(phandle, start, j, sfc);
	} else {
		return j;
	}
}

int artio_particle_seek_to_sfc(artio_fileset *handle, int64_t sfc) {
	int64_t offset;
	artio_particle_file *phandle;
	int file;

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if ( !(handle->open_type & ARTIO_OPEN_PARTICLES) ||
			handle->particle == NULL ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	phandle = handle->particle;

	if (phandle->cache_sfc_begin == -1 || 
			sfc < phandle->cache_sfc_begin || 
			sfc > phandle->cache_sfc_end) {
		return ARTIO_ERR_INVALID_SFC;
	}

	file = artio_particle_find_file(phandle, 0, phandle->num_particle_files, sfc);
	if ( file != phandle->cur_file ) {
		if ( phandle->cur_file != -1 ) {
			artio_file_detach_buffer( phandle->ffh[phandle->cur_file] );
		}
		if ( phandle->buffer_size > 0 ) {
			artio_file_attach_buffer( phandle->ffh[file],
					phandle->buffer, phandle->buffer_size );
		}
		phandle->cur_file = file;
	}
	offset = phandle->sfc_offset_table[sfc - phandle->cache_sfc_begin];
	return artio_file_fseek(phandle->ffh[phandle->cur_file], offset, ARTIO_SEEK_SET);
}

int artio_particle_write_root_cell_begin(artio_fileset *handle, int64_t sfc,
		int * num_particles_per_species) {
	int i;
	int ret;
	artio_particle_file *phandle;

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if (handle->open_mode != ARTIO_FILESET_WRITE || 
			!(handle->open_type & ARTIO_OPEN_PARTICLES) ||
			handle->particle == NULL ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	phandle = handle->particle;

	if ( phandle->cur_sfc != -1 ) {
		return ARTIO_ERR_INVALID_STATE;
	}

	ret = artio_particle_seek_to_sfc(handle, sfc);
	if ( ret != ARTIO_SUCCESS ) return ret;

	ret = artio_file_fwrite(phandle->ffh[phandle->cur_file], num_particles_per_species,
			phandle->num_species, ARTIO_TYPE_INT);
	if ( ret != ARTIO_SUCCESS ) return ret;

	for (i = 0; i < phandle->num_species; i++) {
		phandle->num_particles_per_species[i] = num_particles_per_species[i];
	}

	phandle->cur_sfc = sfc;
	phandle->cur_species = -1;
	phandle->cur_particle = -1;

	return ARTIO_SUCCESS;
}

int artio_particle_write_root_cell_end(artio_fileset *handle) {
	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if (handle->open_mode != ARTIO_FILESET_WRITE || 
			!(handle->open_type & ARTIO_OPEN_PARTICLES) ||
			handle->particle == NULL ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	if ( handle->particle->cur_sfc == -1 || 
			handle->particle->cur_species != -1 ) {
		return ARTIO_ERR_INVALID_STATE;
	}

	handle->particle->cur_sfc = -1;
	return ARTIO_SUCCESS;
}

int artio_particle_write_species_begin(artio_fileset *handle,
		int species) {
	artio_particle_file *phandle;

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if (handle->open_mode != ARTIO_FILESET_WRITE || 
			!(handle->open_type & ARTIO_OPEN_PARTICLES) ||
			handle->particle == NULL ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	phandle = handle->particle;

	if (phandle->cur_sfc == -1 || phandle->cur_species != -1 ) {
		return ARTIO_ERR_INVALID_STATE;
	}

	if ( species < 0 || species >= phandle->num_species) {
		return ARTIO_ERR_INVALID_SPECIES;
	}

	phandle->cur_species = species;
	phandle->cur_particle = 0;

	return ARTIO_SUCCESS;
}

int artio_particle_write_species_end(artio_fileset *handle) {
	artio_particle_file *phandle;

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if (handle->open_mode != ARTIO_FILESET_WRITE || 
			!(handle->open_type & ARTIO_OPEN_PARTICLES) ||
			handle->particle == NULL ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	phandle = handle->particle;

	if (phandle->cur_species == -1 || 
			phandle->cur_particle != 
				phandle->num_particles_per_species[phandle->cur_species]) {
		return ARTIO_ERR_INVALID_STATE;
	}

	phandle->cur_species = -1;
	phandle->cur_particle = -1;

	return ARTIO_SUCCESS;
}

int artio_particle_write_particle(artio_fileset *handle, int64_t pid, int subspecies,
		double * primary_variables, float *secondary_variables) {

	int ret;
	artio_particle_file *phandle;

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if (handle->open_mode != ARTIO_FILESET_WRITE || 
			!(handle->open_type & ARTIO_OPEN_PARTICLES) ||
			handle->particle == NULL ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	phandle = handle->particle;

	if (phandle->cur_species == -1 ||
			phandle->cur_particle >= phandle->num_particles_per_species[phandle->cur_species]) {
		return ARTIO_ERR_INVALID_STATE;
	}

	ret = artio_file_fwrite(phandle->ffh[phandle->cur_file], &pid, 1, ARTIO_TYPE_LONG);
	if ( ret != ARTIO_SUCCESS ) return ret;

	ret = artio_file_fwrite(phandle->ffh[phandle->cur_file], &subspecies, 1, ARTIO_TYPE_INT);
	if ( ret != ARTIO_SUCCESS ) return ret;

	ret = artio_file_fwrite(phandle->ffh[phandle->cur_file], primary_variables,
			phandle->num_primary_variables[phandle->cur_species],
			ARTIO_TYPE_DOUBLE);
	if ( ret != ARTIO_SUCCESS ) return ret;

	ret = artio_file_fwrite(phandle->ffh[phandle->cur_file], secondary_variables,
			phandle->num_secondary_variables[phandle->cur_species],
			ARTIO_TYPE_FLOAT);
	if ( ret != ARTIO_SUCCESS ) return ret;

	phandle->cur_particle++;
	return ARTIO_SUCCESS;
}

/*
 *
 */
int artio_particle_read_root_cell_begin(artio_fileset *handle, int64_t sfc,
		int * num_particles_per_species) {
	int i;
	int ret;
	artio_particle_file *phandle;

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if (handle->open_mode != ARTIO_FILESET_READ || 
			!(handle->open_type & ARTIO_OPEN_PARTICLES) ||
			handle->particle == NULL ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	phandle = handle->particle;

	ret = artio_particle_seek_to_sfc(handle, sfc);
	if ( ret != ARTIO_SUCCESS ) return ret;

	ret = artio_file_fread(phandle->ffh[phandle->cur_file], num_particles_per_species,
			phandle->num_species, ARTIO_TYPE_INT);
	if ( ret != ARTIO_SUCCESS ) return ret;

	for (i = 0; i < phandle->num_species; i++) {
		phandle->num_particles_per_species[i] = num_particles_per_species[i];
	}

	phandle->cur_sfc = sfc;
	phandle->cur_species = -1;
	phandle->cur_particle = 0;
	return ARTIO_SUCCESS;
}

/* Description  */
int artio_particle_read_particle(artio_fileset *handle, int64_t * pid, int *subspecies,
		double * primary_variables, float * secondary_variables) {
	int ret;
	artio_particle_file *phandle;

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if (handle->open_mode != ARTIO_FILESET_READ || 
			!(handle->open_type & ARTIO_OPEN_PARTICLES) ||
			handle->particle == NULL ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	phandle = handle->particle; 

	if (phandle->cur_species == -1 || 
			phandle->cur_particle >= phandle->num_particles_per_species[phandle->cur_species]) {
		return ARTIO_ERR_INVALID_STATE;
	}

	ret = artio_file_fread(phandle->ffh[phandle->cur_file], pid, 1, ARTIO_TYPE_LONG);
	if ( ret != ARTIO_SUCCESS ) return ret;

	ret = artio_file_fread(phandle->ffh[phandle->cur_file], subspecies, 1, ARTIO_TYPE_INT);
	if ( ret != ARTIO_SUCCESS ) return ret;

	ret = artio_file_fread(phandle->ffh[phandle->cur_file], primary_variables,
			phandle->num_primary_variables[phandle->cur_species],
			ARTIO_TYPE_DOUBLE);
	if ( ret != ARTIO_SUCCESS ) return ret;

	ret = artio_file_fread(phandle->ffh[phandle->cur_file], secondary_variables,
			phandle->num_secondary_variables[phandle->cur_species],
			ARTIO_TYPE_FLOAT);
	if ( ret != ARTIO_SUCCESS ) return ret;

	phandle->cur_particle++;
	return ARTIO_SUCCESS;
}

/*
 * Description        Start reading particle species
 */
int artio_particle_read_species_begin(artio_fileset *handle, int species) {
	int i;
	int ret;
	int64_t offset = 0;
	artio_particle_file *phandle;

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if (handle->open_mode != ARTIO_FILESET_READ || 
			!(handle->open_type & ARTIO_OPEN_PARTICLES) ||
			handle->particle == NULL ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	phandle = handle->particle;

	if (phandle->cur_sfc == -1) {
		return ARTIO_ERR_INVALID_STATE;
	}

	if (species < 0 || species >= phandle->num_species) {
		return ARTIO_ERR_INVALID_SPECIES;
	}

	offset = phandle->sfc_offset_table[phandle->cur_sfc - phandle->cache_sfc_begin];
	offset += sizeof(int32_t) * (phandle->num_species);

	for (i = 0; i < species; i++) {
		offset += ( sizeof(int64_t) + sizeof(int) +
					phandle->num_primary_variables[i] * sizeof(double) + 
					phandle->num_secondary_variables[i] * sizeof(float) ) *
						phandle->num_particles_per_species[i];
	}

	ret = artio_file_fseek(phandle->ffh[phandle->cur_file], offset, ARTIO_SEEK_SET);
	if ( ret != ARTIO_SUCCESS ) return ret;

	phandle->cur_species = species;
    phandle->cur_particle = 0;

	return ARTIO_SUCCESS; 
}

/*
 * Description        Do something at the end of each kind of read operation
 */
int artio_particle_read_species_end(artio_fileset *handle) {
	artio_particle_file *phandle;

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if (handle->open_mode != ARTIO_FILESET_READ || 
			!(handle->open_type & ARTIO_OPEN_PARTICLES) ||
			handle->particle == NULL ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	phandle = handle->particle;

	if (phandle->cur_species == -1) {
		return ARTIO_ERR_INVALID_STATE;
	}

	phandle->cur_species = -1;
	phandle->cur_particle = 0;

	return ARTIO_SUCCESS;
}

int artio_particle_read_root_cell_end(artio_fileset *handle) {
	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if (handle->open_mode != ARTIO_FILESET_READ || 
			!(handle->open_type & ARTIO_OPEN_PARTICLES) ||
			handle->particle == NULL ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	if ( handle->particle->cur_sfc == -1 ) {
		return ARTIO_ERR_INVALID_STATE;
	}

	handle->particle->cur_sfc = -1;
	return ARTIO_SUCCESS;
}

int artio_particle_read_selection(artio_fileset *handle,
		artio_selection *selection, artio_particle_callback callback, 
		void *params ) {

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if (handle->open_mode != ARTIO_FILESET_READ || 
			!(handle->open_type & ARTIO_OPEN_PARTICLES) ||
			handle->particle == NULL ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	return artio_particle_read_selection_species( handle,
		selection, 0, handle->particle->num_species-1,
		callback, params );		
}

int artio_particle_read_selection_species( artio_fileset *handle,
        artio_selection *selection, int start_species, int end_species,
		artio_particle_callback callback, void *params ) {
	int ret;
	int64_t start, end;

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
    }
    
	if (handle->open_mode != ARTIO_FILESET_READ || 
			!(handle->open_type & ARTIO_OPEN_PARTICLES) ||
			handle->particle == NULL ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	artio_selection_iterator_reset( selection );
	while ( artio_selection_iterator( selection,
				handle->num_root_cells, &start, &end ) == ARTIO_SUCCESS ) {
		ret = artio_particle_read_sfc_range_species( handle,
			start, end, start_species, end_species, callback, params );
		if ( ret != ARTIO_SUCCESS ) return ret;
	}

	return ARTIO_SUCCESS;
}

int artio_particle_read_sfc_range(artio_fileset *handle, 
		int64_t sfc1, int64_t sfc2, 
		artio_particle_callback callback,
		void *params ) {
	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}
	
	if (handle->open_mode != ARTIO_FILESET_READ || 
            !(handle->open_type & ARTIO_OPEN_PARTICLES) ||
			handle->particle == NULL ) {
        return ARTIO_ERR_INVALID_FILESET_MODE;
    }

	return artio_particle_read_sfc_range_species( handle,
		sfc1, sfc2, 0, handle->particle->num_species-1,
		callback, params );
}

int artio_particle_read_sfc_range_species(artio_fileset *handle, 
        int64_t sfc1, int64_t sfc2, 
        int start_species, int end_species, 
        artio_particle_callback callback,
		void *params ) {
	int64_t sfc;
	int particle, species;
	int *num_particles_per_species;
	artio_particle_file *phandle;
	int64_t pid = 0l;
	int subspecies;
	double * primary_variables = NULL;
	float * secondary_variables = NULL;
	int num_primary, num_secondary;
	int ret;

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if (handle->open_mode != ARTIO_FILESET_READ || 
			!(handle->open_type & ARTIO_OPEN_PARTICLES) ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	phandle = handle->particle;

	if ( start_species < 0 || start_species > end_species || 
			end_species > phandle->num_species-1 ) {
		return ARTIO_ERR_INVALID_SPECIES;
	}

	num_particles_per_species = (int *)malloc(phandle->num_species * sizeof(int));
	if ( num_particles_per_species == NULL ) {
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}

	ret = artio_particle_cache_sfc_range(handle, sfc1, sfc2);
	if ( ret != ARTIO_SUCCESS ) {
		free( num_particles_per_species );
		return ret;
	}

	num_primary = num_secondary = 0;
	for ( species = start_species; species <= end_species; species++ ) {
		num_primary = MAX( phandle->num_primary_variables[species], num_primary );
		num_secondary = MAX( phandle->num_secondary_variables[species], num_secondary );
	}

	primary_variables = (double *)malloc(num_primary * sizeof(double));
	if ( primary_variables == NULL ) {
		free( num_particles_per_species );
		return ARTIO_ERR_MEMORY_ALLOCATION;
    }

	secondary_variables = (float *)malloc(num_secondary * sizeof(float));
	if ( secondary_variables == NULL ) {
		free( num_particles_per_species );
		free( primary_variables );
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}

	for ( sfc = sfc1; sfc <= sfc2; sfc++ ) {
		ret = artio_particle_read_root_cell_begin(handle, sfc,
				num_particles_per_species);
		if ( ret != ARTIO_SUCCESS ) {
			free( num_particles_per_species );
			free( primary_variables );
			free( secondary_variables );
			return ret;
		}

		for ( species = start_species; species <= end_species; species++) {
			ret = artio_particle_read_species_begin(handle, species);
			if ( ret != ARTIO_SUCCESS ) {
				free( num_particles_per_species );
				free( primary_variables );
				free( secondary_variables );
				return ret;
			}

			for (particle = 0; particle < num_particles_per_species[species]; particle++) {
				ret = artio_particle_read_particle(handle, 
						&pid, 
						&subspecies,
						primary_variables,
						secondary_variables);
				if ( ret != ARTIO_SUCCESS ) {
					free( num_particles_per_species );
					free( primary_variables );
					free( secondary_variables );
					return ret;
				}

				callback(sfc, species, subspecies,
						pid, 
						primary_variables,
						secondary_variables,
						params  );
			}
			artio_particle_read_species_end(handle);
		}
		artio_particle_read_root_cell_end(handle);
	}

	free(primary_variables);
	free(secondary_variables);
	free(num_particles_per_species);

	return ARTIO_SUCCESS;
}
