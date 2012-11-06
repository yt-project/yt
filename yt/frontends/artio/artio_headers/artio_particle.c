/*
 * artio_particle.c
 *
 *  Created on: May 10, 2011
 *      Author: eric
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#ifdef ARTIO_MPI
#include "artio_mpi.h"
#endif /* ARTIO_MPI */

#include "artio.h"
#include "artio_internal.h"

int artio_particle_find_file(artio_file handle, int start, int end, int64_t sfc);

/*
 * Open existing particle files and add to fileset
 */
int artio_particle_open(artio_file handle) {
	int i;
	char filename[256];
	int first_file, last_file;
	int mode;

	if ( !(handle->open_type & ARTIO_OPEN_PARTICLES) ||
			handle->open_mode != ARTIO_FILESET_READ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	artio_particle_file phandle = (artio_particle_file)malloc(sizeof(struct artio_particle_file_struct));
	handle->particle = phandle;

	artio_parameter_get_int(handle, "num_particle_species", &phandle->num_species);

	phandle->num_primary_variables = (int *)malloc(sizeof(int) * phandle->num_species);
	phandle->num_secondary_variables = (int *)malloc(sizeof(int) * phandle->num_species);
	phandle->num_particles_per_species
	  = (int *)malloc(phandle->num_species * sizeof(int));

	artio_parameter_get_int_array(handle, "num_primary_variables",
			phandle->num_species, phandle->num_primary_variables);
	artio_parameter_get_int_array(handle, "num_secondary_variables",
			phandle->num_species, phandle->num_secondary_variables);

	artio_parameter_get_int(handle, "num_particle_files", &phandle->num_particle_files);

	phandle->file_sfc_index = (int64_t *)malloc(sizeof(int64_t) * (phandle->num_particle_files + 1));
	artio_parameter_get_long_array(handle, "particle_file_sfc_index",
			phandle->num_particle_files + 1, phandle->file_sfc_index);

	first_file = artio_particle_find_file(handle, 0, 
			phandle->num_particle_files, handle->proc_sfc_begin);
	last_file = artio_particle_find_file(handle, first_file,
			phandle->num_particle_files, handle->proc_sfc_end);

	/* allocate file handles */
	phandle->ffh = (artio_fh *)malloc(phandle->num_particle_files * sizeof(artio_fh));

	/* open files on all processes */
	for (i = 0; i < phandle->num_particle_files; i++) {
		snprintf(filename, 255, "%s.p%03d", handle->file_prefix, i);

		mode = ARTIO_MODE_READ;
		if (i >= first_file && i <= last_file) {
			mode |= ARTIO_MODE_ACCESS;
		}
		if (handle->param_list.endian_swap) {
			mode |= ARTIO_MODE_ENDIAN_SWAP;
		}
		phandle->ffh[i] = artio_file_fopen(filename, mode, handle->context);

		if ( phandle->ffh[i] == NULL ) {
			return ARTIO_ERR_PARTICLE_FILE_NOT_FOUND;
		}
	}

	phandle->cache_sfc_begin = -1;
	phandle->cache_sfc_end = -1;
	phandle->sfc_offset_table = NULL;

	phandle->cur_file = -1;
	phandle->cur_sfc = -1;
	phandle->cur_species = -1;
	phandle->cur_particle = -1;

	return ARTIO_SUCCESS;
}

int artio_fileset_add_particles( artio_file handle,                                                                                  
        int num_particle_files, int allocation_strategy,
        int num_species, char ** species_labels,
        int * num_primary_variables,
        int * num_secondary_variables,
        char *** primary_variable_labels_per_species,
        char *** secondary_variable_labels_per_species,
	int * num_particles_per_species_per_root_tree ) {

	int i, j, k;

	int *file_parent;
	int64_t *file_sfc_index;
	int64_t cur;
	int64_t first_file_sfc, last_file_sfc;

	int first_file, last_file;
	char filename[256];
	char species_label[64];
	int mode;
	int64_t *local_particles_per_species, *total_particles_per_species;

	if ( handle->open_mode != ARTIO_FILESET_WRITE ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	if ( handle->open_type & ARTIO_OPEN_PARTICLES ) {                                                         
		return ARTIO_ERR_DATA_EXISTS;
	}
	handle->open_type |= ARTIO_OPEN_PARTICLES;

	/* compute total number of particles per species */
	local_particles_per_species = (int64_t *)malloc( num_species * sizeof(int64_t));
	total_particles_per_species = (int64_t *)malloc( num_species * sizeof(int64_t));

	for ( i = 0; i < num_species; i++ ) {
		local_particles_per_species[i] = 0;
	}

	for ( i = 0; i < handle->proc_sfc_end-handle->proc_sfc_begin+1; i++ ) {
		for ( j = 0; j < num_species; j++ ) {
			local_particles_per_species[j] += num_particles_per_species_per_root_tree[num_species*i+j];
		}
	}

#ifdef ARTIO_MPI
	MPI_Allreduce( local_particles_per_species, total_particles_per_species, num_species,
			MPI_LONG_LONG_INT, MPI_SUM, handle->context->comm );
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

	artio_particle_file phandle = (artio_particle_file)malloc(sizeof(struct artio_particle_file_struct));
	handle->particle = phandle;

	file_parent = (int *)malloc(sizeof(int) * num_particle_files);
	file_sfc_index = (int64_t *)malloc(sizeof(int64_t) * (num_particle_files + 1));
	if ( file_parent == NULL || file_sfc_index == NULL ) {
		fprintf(stderr, "ERROR allocating memory!");
		exit(1);
	}

	switch (allocation_strategy) {
		case ARTIO_ALLOC_EQUAL_PROC:
			if (num_particle_files > handle->num_procs) {
				return ARTIO_ERR_INVALID_FILE_NUMBER;
			}

			file_parent[0] = 0;
			file_sfc_index[0] = 0;

			for (i = 1; i < num_particle_files; i++) {
				file_parent[i] = file_parent[i - 1] + 
							ceil((handle->num_procs - file_parent[i - 1]) / 
								(num_particle_files - i + 1));
				file_sfc_index[i] = handle->proc_sfc_index[file_parent[i]];
			}
			file_sfc_index[num_particle_files] = handle->proc_sfc_index[handle->num_procs];
			break;
		case ARTIO_ALLOC_EQUAL_SFC:
			file_sfc_index[0] = 0;
			file_parent[0] = 0;

			for (i = 1; i < num_particle_files; i++) {
				file_sfc_index[i] = file_sfc_index[i - 1] + 
							ceil((handle->num_root_cells - file_sfc_index[i - 1]) / 
								(num_particle_files - i + 1));
				file_parent[i] = file_parent[i - 1];
				while (file_parent[i] < handle->num_procs && 
						handle->proc_sfc_index[file_parent[i] + 1] <= file_sfc_index[i]) {
					file_parent[i]++;
				}
			}
			file_sfc_index[num_particle_files] = handle->num_root_cells;
			break;
		default:
			return ARTIO_ERR_INVALID_ALLOC_STRATEGY;		
	}

	phandle->num_particle_files = num_particle_files;
	phandle->num_species = num_species;
	phandle->file_sfc_index = file_sfc_index;
	phandle->num_primary_variables = (int *)malloc(sizeof(int) * num_species);
	phandle->num_secondary_variables = (int *)malloc(sizeof(int) * num_species);
	phandle->num_particles_per_species = 
	  (int *)malloc(phandle->num_species * sizeof(int));

	for (i = 0; i < num_species; i++) {
		phandle->num_primary_variables[i] = num_primary_variables[i];
		phandle->num_secondary_variables[i] = num_secondary_variables[i];
	}

	/* allocate space for sfc offset cache */
	phandle->cache_sfc_begin = handle->proc_sfc_begin;
	phandle->cache_sfc_end = handle->proc_sfc_end;
	phandle->sfc_offset_table = 
	  (int64_t *)malloc( (handle->proc_sfc_end - handle->proc_sfc_begin + 1) * sizeof(int64_t));

	/* allocate file handles */
	phandle->ffh = (artio_fh *)malloc(num_particle_files * sizeof(artio_fh));

	/* open file handles */
	first_file = artio_particle_find_file(handle, 0, 
					num_particle_files, handle->proc_sfc_begin);
	last_file = artio_particle_find_file(handle, first_file, 
					num_particle_files, handle->proc_sfc_end);

	for (i = 0; i < num_particle_files; i++) {
		snprintf(filename, 255, "%s.p%03d", handle->file_prefix, i);
		mode = ARTIO_MODE_WRITE;
		if (i >= first_file && i <= last_file) {
			mode |= ARTIO_MODE_ACCESS;
		}

		phandle->ffh[i] = artio_file_fopen(filename, mode, handle->context);

		if ( phandle->ffh[i] == NULL ) {
			return ARTIO_ERR_FILE_CREATE;
		}

		/* write sfc offset header if we contribute to this file */
		if (i >= first_file && i <= last_file) {
#ifdef ARTIO_MPI
			if (file_parent[i] == handle->rank) {
				cur = (file_sfc_index[i+1] - file_sfc_index[i]) * sizeof(int64_t);
			} else {
				/* obtain offset from previous process */
				MPI_Recv( &cur, 1, MPI_LONG_LONG_INT, handle->rank - 1, i,
						handle->context->comm, MPI_STATUS_IGNORE );
			}
#else
			cur = (file_sfc_index[i+1] - file_sfc_index[i]) * sizeof(int64_t);
#endif /* ARTIO_MPI */

			first_file_sfc = MAX( phandle->cache_sfc_begin, file_sfc_index[i] );
			last_file_sfc = MIN( phandle->cache_sfc_end, file_sfc_index[i+1]-1 );

			for (j = first_file_sfc - phandle->cache_sfc_begin; 
					j < last_file_sfc - phandle->cache_sfc_begin + 1; j++) {

				phandle->sfc_offset_table[j] = cur;
				cur += sizeof(int) * num_species;

				for (k = 0; k < num_species; k++) {
					cur += num_particles_per_species_per_root_tree[j*num_species+k]
							* (sizeof(int64_t) + sizeof(int) +
									num_primary_variables[k] * sizeof(double) + 
									num_secondary_variables[k] * sizeof(float));
				}
			}

#ifdef ARTIO_MPI
			if ( file_sfc_index[i+1] > handle->proc_sfc_end+1 ) {
				MPI_Send( &cur, 1, MPI_LONG_LONG_INT, handle->rank+1, i, handle->context->comm );
			}
#endif /* ARTIO_MPI */

			/* seek and write our portion of sfc table */
			artio_file_fseek(phandle->ffh[i], 
					(first_file_sfc - file_sfc_index[i]) * sizeof(int64_t), ARTIO_SEEK_SET);
			artio_file_fwrite(phandle->ffh[i],
					&phandle->sfc_offset_table[first_file_sfc - phandle->cache_sfc_begin],
					last_file_sfc - first_file_sfc + 1, ARTIO_TYPE_LONG);
		}
	}

	free(file_parent);

	phandle->cur_file = -1;
	phandle->cur_sfc = -1;
	phandle->cur_species = -1;
	phandle->cur_particle = -1;

	artio_parameter_set_long_array(handle, "particle_file_sfc_index",
			handle->particle->num_particle_files + 1, 
			handle->particle->file_sfc_index);

	return ARTIO_SUCCESS;
}

int artio_particle_close(artio_file handle) {
	int i;
	artio_particle_file phandle = handle->particle;
	for (i = 0; i < phandle->num_particle_files; i++) {
		artio_file_fclose(phandle->ffh[i]);
	}
	free(phandle->ffh);

	if (phandle->sfc_offset_table != NULL) {
		free(phandle->sfc_offset_table);
	}

	free(phandle->num_particles_per_species);
	free(phandle->num_primary_variables);
	free(phandle->num_secondary_variables);
	free(phandle->file_sfc_index);
	return ARTIO_SUCCESS;
}

int artio_particle_cache_sfc_range(artio_file handle, 
		int64_t start, int64_t end) {
	int i;
	int first_file, last_file, min, count, cur;
	artio_particle_file phandle = handle->particle;

	if ( start > end || start < handle->proc_sfc_begin || 
			end > handle->proc_sfc_end) {
		return ARTIO_ERR_INVALID_SFC_RANGE;
	}

	if (handle->open_mode != ARTIO_FILESET_READ || 
			!(handle->open_type & ARTIO_OPEN_PARTICLES) ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	if ( phandle->sfc_offset_table != NULL ) {
		free(phandle->sfc_offset_table);
	}

	first_file = artio_particle_find_file(handle, 0, 
			phandle->num_particle_files, start);
	last_file = artio_particle_find_file(handle, first_file,
			phandle->num_particle_files, end);

	phandle->cache_sfc_begin = start;
	phandle->cache_sfc_end = end;
	phandle->sfc_offset_table = (int64_t *)malloc(sizeof(int64_t) * (end - start + 1));

	cur = 0;
	for (i = first_file; i <= last_file; i++) {
		min = MAX( 0, start - phandle->file_sfc_index[i] );
		count = MIN( phandle->file_sfc_index[i+1], end+1 )
				- MAX( start, phandle->file_sfc_index[i] );

		artio_file_fseek(phandle->ffh[i], 
				sizeof(int64_t) * min,
				ARTIO_SEEK_SET);
		artio_file_fread(phandle->ffh[i], 
				&phandle->sfc_offset_table[cur],
				count, ARTIO_TYPE_LONG);
		cur += count;
	}

	return ARTIO_SUCCESS;
}

int artio_particle_find_file(artio_file handle,
		int start, int end, int64_t sfc) {
	int j;
	artio_particle_file phandle = handle->particle;

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
		return artio_particle_find_file(handle, j, end, sfc);
	} else if (sfc < phandle->file_sfc_index[j]) {
		return artio_particle_find_file(handle, start, j, sfc);
	} else {
		return j;
	}
}

int artio_particle_seek_to_sfc(artio_file handle, int64_t sfc) {
	int64_t offset;
	artio_particle_file phandle = handle->particle;

	if (phandle->cache_sfc_begin == -1 || 
			sfc < phandle->cache_sfc_begin || 
			sfc > phandle->cache_sfc_end) {
		return ARTIO_ERR_INVALID_SFC;
	}

	phandle->cur_file = artio_particle_find_file(handle, 0, phandle->num_particle_files, sfc);
	offset = phandle->sfc_offset_table[sfc - phandle->cache_sfc_begin];
	artio_file_fseek(phandle->ffh[phandle->cur_file], offset, ARTIO_SEEK_SET);

	return ARTIO_SUCCESS;
}

int artio_particle_write_root_cell_begin(artio_file handle, int64_t sfc,
		int * num_particles_per_species) {
	int i;
	artio_particle_file phandle=handle->particle;

	if (handle->open_mode != ARTIO_FILESET_WRITE || 
			!(handle->open_type & ARTIO_OPEN_PARTICLES) ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	if ( phandle->cur_sfc != -1 ) {
		return ARTIO_ERR_INVALID_STATE;
	}

	artio_particle_seek_to_sfc(handle, sfc);
	artio_file_fwrite(phandle->ffh[phandle->cur_file], num_particles_per_species,
			phandle->num_species, ARTIO_TYPE_INT);
	for (i = 0; i < phandle->num_species; i++) {
		phandle->num_particles_per_species[i] = num_particles_per_species[i];
	}

	phandle->cur_sfc = sfc;
	phandle->cur_species = -1;
	phandle->cur_particle = -1;

	return ARTIO_SUCCESS;
}

int artio_particle_write_root_cell_end(artio_file handle) {
	if (handle->open_mode != ARTIO_FILESET_WRITE || 
			!(handle->open_type & ARTIO_OPEN_PARTICLES) ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	if ( handle->particle->cur_sfc == -1 || 
			handle->particle->cur_species != -1 ) {
		return ARTIO_ERR_INVALID_STATE;
	}
;
	handle->particle->cur_sfc = -1;
	return ARTIO_SUCCESS;
}

int artio_particle_write_species_begin(artio_file handle,
		int species) {
	artio_particle_file phandle = handle->particle;

	if (handle->open_mode != ARTIO_FILESET_WRITE || 
			!(handle->open_type & ARTIO_OPEN_PARTICLES) ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

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

int artio_particle_write_species_end(artio_file handle) {
	artio_particle_file phandle = handle->particle;

	if (handle->open_mode != ARTIO_FILESET_WRITE || 
			!(handle->open_type & ARTIO_OPEN_PARTICLES) ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	if (phandle->cur_species == -1 || 
			phandle->cur_particle != 
				phandle->num_particles_per_species[phandle->cur_species]) {
		return ARTIO_ERR_INVALID_STATE;
	}

	phandle->cur_species = -1;
	phandle->cur_particle = -1;

	return ARTIO_SUCCESS;
}

int artio_particle_write_particle(artio_file handle, int64_t pid, int subspecies,
		double * primary_variables, float *secondary_variables) {

	artio_particle_file phandle = handle->particle;

	if (handle->open_mode != ARTIO_FILESET_WRITE || 
			!(handle->open_type & ARTIO_OPEN_PARTICLES) ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	if (phandle->cur_species == -1 ||
			phandle->cur_particle >= phandle->num_particles_per_species[phandle->cur_species]) {
		return ARTIO_ERR_INVALID_STATE;
	}

	artio_file_fwrite(phandle->ffh[phandle->cur_file], &pid, 1, ARTIO_TYPE_LONG);
	artio_file_fwrite(phandle->ffh[phandle->cur_file], &subspecies, 1, ARTIO_TYPE_INT);

	artio_file_fwrite(phandle->ffh[phandle->cur_file], primary_variables,
			phandle->num_primary_variables[phandle->cur_species],
			ARTIO_TYPE_DOUBLE);

	artio_file_fwrite(phandle->ffh[phandle->cur_file], secondary_variables,
			phandle->num_secondary_variables[phandle->cur_species],
			ARTIO_TYPE_FLOAT);

	phandle->cur_particle++;
	return ARTIO_SUCCESS;
}

/*
 *
 */
int artio_particle_read_root_cell_begin(artio_file handle, int64_t sfc,
		int * num_particles_per_species) {
	int i;
	artio_particle_file phandle = handle->particle;

	if (handle->open_mode != ARTIO_FILESET_READ || 
			!(handle->open_type & ARTIO_OPEN_PARTICLES) ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	artio_particle_seek_to_sfc(handle, sfc);
	artio_file_fread(phandle->ffh[phandle->cur_file], num_particles_per_species,
			phandle->num_species, ARTIO_TYPE_INT);

	for (i = 0; i < phandle->num_species; i++) {
		phandle->num_particles_per_species[i] = num_particles_per_species[i];
	}

	phandle->cur_sfc = sfc;
	phandle->cur_species = -1;
	phandle->cur_particle = 0;
	return ARTIO_SUCCESS;
}

/* Description  */
int artio_particle_read_particle(artio_file handle, int64_t * pid, int *subspecies,
		double * primary_variables, float * secondary_variables) {
	artio_particle_file phandle = handle->particle;

	if (handle->open_mode != ARTIO_FILESET_READ || 
			!(handle->open_type & ARTIO_OPEN_PARTICLES) ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	if (phandle->cur_species == -1 || 
			phandle->cur_particle >= phandle->num_particles_per_species[phandle->cur_species]) {
		return ARTIO_ERR_INVALID_STATE;
	}

	artio_file_fread(phandle->ffh[phandle->cur_file], pid, 1, ARTIO_TYPE_LONG);
	artio_file_fread(phandle->ffh[phandle->cur_file], subspecies, 1, ARTIO_TYPE_INT);
	artio_file_fread(phandle->ffh[phandle->cur_file], primary_variables,
			phandle->num_primary_variables[phandle->cur_species],
			ARTIO_TYPE_DOUBLE);
	artio_file_fread(phandle->ffh[phandle->cur_file], secondary_variables,
			phandle->num_secondary_variables[phandle->cur_species],
			ARTIO_TYPE_FLOAT);

	phandle->cur_particle++;
	return ARTIO_SUCCESS;
}

/*
 * Description        Start reading particle species
 */
int artio_particle_read_species_begin(artio_file handle, int species) {
	int i;
	int64_t offset = 0;
	artio_particle_file phandle = handle->particle;

	if (handle->open_mode != ARTIO_FILESET_READ || 
			!(handle->open_type & ARTIO_OPEN_PARTICLES) ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

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

	artio_file_fseek(phandle->ffh[phandle->cur_file], offset, ARTIO_SEEK_SET);
	phandle->cur_species = species;
	phandle->cur_particle = 0;

	return ARTIO_SUCCESS;
}

/*
 * Description        Do something at the end of each kind of read operation
 */
int artio_particle_read_species_end(artio_file handle) {
	artio_particle_file phandle = handle->particle;

	if (handle->open_mode != ARTIO_FILESET_READ || 
			!(handle->open_type & ARTIO_OPEN_PARTICLES) ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	if (phandle->cur_species == -1) {
		return ARTIO_ERR_INVALID_STATE;
	}

	phandle->cur_species = -1;
	phandle->cur_particle = 0;

	return ARTIO_SUCCESS;
}

int artio_particle_read_root_cell_end(artio_file handle) {
	if (handle->open_mode != ARTIO_FILESET_READ || 
			!(handle->open_type & ARTIO_OPEN_PARTICLES) ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	if ( handle->particle->cur_sfc == -1 ) {
		return ARTIO_ERR_INVALID_STATE;
	}

	handle->particle->cur_sfc = -1;
	return ARTIO_SUCCESS;
}

int artio_particle_read_sfc_range(artio_file handle, 
		int64_t sfc1, int64_t sfc2, 
		int start_species, int end_species, 
		ParticleCallBack callback) {

	int64_t sfc;
	int particle, species;
	int *num_particles_per_species;
	artio_particle_file phandle = handle->particle;
	int64_t pid = 0l;
	int subspecies;
	double * primary_variables = NULL;
	float * secondary_variables = NULL;
	int num_primary, num_secondary;
	int ret;

	if (handle->open_mode != ARTIO_FILESET_READ || 
			!(handle->open_type & ARTIO_OPEN_PARTICLES) ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	if ( start_species < 0 || start_species > end_species || 
			end_species > phandle->num_species-1 ) {
		return ARTIO_ERR_INVALID_SPECIES;
	}

	num_particles_per_species = (int *)malloc(phandle->num_species * sizeof(int));
	artio_particle_cache_sfc_range(handle, sfc1, sfc2);

	num_primary = num_secondary = 0;
	for ( species = start_species; species <= end_species; species++ ) {
		num_primary = MAX( phandle->num_primary_variables[species], num_primary );
		num_secondary = MAX( phandle->num_secondary_variables[species], num_secondary );
	}

	primary_variables = (double *)malloc(num_primary * sizeof(double));
	secondary_variables = (float *)malloc(num_secondary * sizeof(float));

	for ( sfc = sfc1; sfc <= sfc2; sfc++ ) {
		artio_particle_read_root_cell_begin(handle, sfc,
				num_particles_per_species);

		for ( species = start_species; species <= end_species; species++) {
			artio_particle_read_species_begin(handle, species);

			for (particle = 0; particle < num_particles_per_species[species]; particle++) {
				ret = artio_particle_read_particle(handle, 
						&pid, 
						&subspecies,
						primary_variables,
						secondary_variables);

				if ( ret != ARTIO_SUCCESS ) {
					return ARTIO_ERR_PARTICLE_CORRUPTED;
				} 

				callback(pid, 
						primary_variables,
						secondary_variables, 
						species, subspecies, sfc);
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
