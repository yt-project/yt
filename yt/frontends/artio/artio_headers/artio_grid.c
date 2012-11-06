/*
 * artio_grid.c
 *
 *  Created on: May 10, 2011
 *      Author: Yongen Yu
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

int artio_grid_find_file(artio_file handle, int start, int end, int64_t sfc);

/*
 * Open grid component of the fileset
 */
int artio_grid_open(artio_file handle) {
	int i;
	char filename[256];
	int first_file, last_file;
	int mode;

	/* check that the fileset contains a grid component */
	if ( !(handle->open_type & ARTIO_OPEN_GRID) ||
			handle->open_mode != ARTIO_FILESET_READ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	artio_grid_file ghandle = (artio_grid_file)malloc(sizeof(struct artio_grid_file_struct));
	handle->grid = ghandle;

	/* load grid parameters from header file (should be doing error handling...) */
    artio_parameter_get_int(handle, "num_grid_variables", &ghandle->num_grid_variables);
    artio_parameter_get_int(handle, "num_grid_files", &ghandle->num_grid_files);
    ghandle->file_sfc_index = (int64_t *)malloc(sizeof(int64_t) * (ghandle->num_grid_files + 1));
    artio_parameter_get_long_array(handle, "grid_file_sfc_index",
            ghandle->num_grid_files + 1, ghandle->file_sfc_index);
    artio_parameter_get_int(handle, "grid_max_level",
            &ghandle->file_max_level);

    ghandle->octs_per_level = (int *)malloc(ghandle->file_max_level * sizeof(int));

	/* allocate file handles */
    ghandle->ffh = (artio_fh *)malloc(ghandle->num_grid_files * sizeof(artio_fh));

	first_file = artio_grid_find_file(handle, 0, 
			ghandle->num_grid_files, handle->proc_sfc_begin);
	last_file = artio_grid_find_file(handle, first_file,
			ghandle->num_grid_files, handle->proc_sfc_end);

	/* open files on all processes */
	for (i = 0; i < ghandle->num_grid_files; i++) {
		snprintf(filename, 255, "%s.g%03d", handle->file_prefix, i);

		mode = ARTIO_MODE_READ;
		if (i >= first_file && i <= last_file) {
			mode |= ARTIO_MODE_ACCESS;
		}

		if(handle->param_list.endian_swap) {
			mode |= ARTIO_MODE_ENDIAN_SWAP;
		}

		ghandle->ffh[i] = artio_file_fopen(filename, mode, handle->context);

		if ( ghandle->ffh[i] == NULL ) {
			return ARTIO_ERR_GRID_FILE_NOT_FOUND;
		}
	}

	ghandle->cache_sfc_begin = -1;
	ghandle->cache_sfc_end = -1;
	ghandle->sfc_offset_table = NULL;

	ghandle->cur_file = -1;
	ghandle->cur_sfc = -1;
	ghandle->cur_level = -1;
	ghandle->cur_num_levels = -1;
	ghandle->cur_octs = -1;

	handle->grid = ghandle;
	return ARTIO_SUCCESS;
}

int artio_fileset_add_grid(artio_file handle, 
        int num_grid_files, int allocation_strategy, 
        int num_grid_variables, 
        char ** grid_variable_labels,
        int * num_levels_per_root_tree,
        int * num_octs_per_root_tree ) {

	int i, j;
	int *file_parent;
	int64_t *file_sfc_index;
	int file_max_level, local_max_level;
	int64_t cur, sfc;
	int64_t first_file_sfc, last_file_sfc;
	int first_file, last_file;
	char filename[256];
	int mode;

	if ( handle->open_mode != ARTIO_FILESET_WRITE ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	if ( handle->open_type & ARTIO_OPEN_GRID) {
		return ARTIO_ERR_DATA_EXISTS;
	}
	handle->open_type |= ARTIO_OPEN_GRID;

	artio_parameter_set_int(handle, "num_grid_files", num_grid_files);
	artio_parameter_set_int(handle, "num_grid_variables", num_grid_variables);
	artio_parameter_set_string_array(handle, "grid_variable_labels",
			num_grid_variables, grid_variable_labels);

	artio_grid_file ghandle = (artio_grid_file)malloc(sizeof(struct artio_grid_file_struct));
	handle->grid = ghandle;

	file_parent = (int *)malloc(sizeof(int) * num_grid_files);
	file_sfc_index = (int64_t *)malloc(sizeof(int64_t) * (num_grid_files + 1));

	/* compute global maximum level */
	local_max_level = 0;
	for (sfc = 0; sfc < handle->proc_sfc_end - handle->proc_sfc_begin + 1; sfc++) {
		if (num_levels_per_root_tree[sfc] > local_max_level) {
			local_max_level = num_levels_per_root_tree[sfc];
		}
	}

#ifdef ARTIO_MPI
	MPI_Allreduce( &local_max_level, &file_max_level, 1, MPI_INT, MPI_MAX, handle->context->comm );
#else
	file_max_level = local_max_level;
#endif /* ARTIO_MPI */

	switch (allocation_strategy) {
		case ARTIO_ALLOC_EQUAL_PROC:
			if (num_grid_files > handle->num_procs) {
				return ARTIO_ERR_INVALID_FILE_NUMBER;
			}

			file_parent[0] = 0;
			file_sfc_index[0] = 0;
			for (i = 1; i < num_grid_files; i++) {
				file_parent[i] = file_parent[i - 1] + 
						ceil((handle->num_procs - file_parent[i - 1]) / (num_grid_files - i + 1));
				file_sfc_index[i] = handle->proc_sfc_index[file_parent[i]];
			}
			file_sfc_index[num_grid_files] = handle->proc_sfc_index[handle->num_procs];
			break;
		case ARTIO_ALLOC_EQUAL_SFC:
			file_sfc_index[0] = 0;
			file_parent[0] = 0;
			for (i = 1; i < num_grid_files; i++) {
				file_sfc_index[i] = file_sfc_index[i - 1] + 
					ceil((handle->num_root_cells - file_sfc_index[i - 1]) / (num_grid_files - i + 1));
				file_parent[i] = file_parent[i - 1];
				while (file_parent[i] < handle->num_procs && 
						handle->proc_sfc_index[file_parent[i]+ 1] <= file_sfc_index[i]) {
					file_parent[i]++;
				}
			}
			file_sfc_index[num_grid_files] = handle->num_root_cells;
			break;
		default:
			return ARTIO_ERR_INVALID_ALLOC_STRATEGY;
	}

	ghandle->num_grid_files = num_grid_files;
	ghandle->file_sfc_index = file_sfc_index;
	ghandle->num_grid_variables = num_grid_variables;
	ghandle->file_max_level = file_max_level;

	/* allocate space for sfc offset cache */
	ghandle->cache_sfc_begin = handle->proc_sfc_begin;
	ghandle->cache_sfc_end = handle->proc_sfc_end;
	ghandle->sfc_offset_table = 
	  (int64_t *)malloc((ghandle->cache_sfc_end - ghandle->cache_sfc_begin + 1) * sizeof(int64_t));

	ghandle->octs_per_level = (int *)malloc(ghandle->file_max_level * sizeof(int));

	/* allocate file handles */
	ghandle->ffh = (artio_fh *)malloc(num_grid_files * sizeof(artio_fh));

	/* open file handles */
	first_file = artio_grid_find_file(handle, 0, num_grid_files, 
					handle->proc_sfc_begin);
	last_file = artio_grid_find_file(handle, first_file, num_grid_files, 
					handle->proc_sfc_end);
	
	for (i = 0; i < num_grid_files; i++) {
		snprintf(filename, 255, "%s.g%03d", handle->file_prefix, i);

		mode = ARTIO_MODE_WRITE;
		if (i >= first_file && i <= last_file) {
			mode |= ARTIO_MODE_ACCESS;
		}

		ghandle->ffh[i] = artio_file_fopen(filename, mode, handle->context);

		/* write sfc offset header if we contribute to this file */
		if (i >= first_file && i <= last_file) {
#ifdef ARTIO_MPI
			if (file_parent[i] == handle->rank) {
				cur = (file_sfc_index[i + 1] - file_sfc_index[i]) * sizeof(int64_t);
			} else {
				/* obtain offset from previous process */
				MPI_Recv( &cur, 1, MPI_LONG_LONG_INT, handle->rank - 1, i,
						handle->context->comm, MPI_STATUS_IGNORE );
			}
#else
			cur = (file_sfc_index[i + 1] - file_sfc_index[i]) * sizeof(int64_t);
#endif /* ARTIO_MPI */

			first_file_sfc = MAX( handle->proc_sfc_begin, file_sfc_index[i] );
			last_file_sfc = MIN( handle->proc_sfc_end, file_sfc_index[i+1]-1 );

			for (j = first_file_sfc - ghandle->cache_sfc_begin; 
					j < last_file_sfc - ghandle->cache_sfc_begin + 1; j++) {
				ghandle->sfc_offset_table[j] = cur;
				cur += sizeof(float) * ghandle->num_grid_variables + sizeof(int) * (1
						+ num_levels_per_root_tree[j])
						+ num_octs_per_root_tree[j] * 8 * (sizeof(float)
								* ghandle->num_grid_variables + sizeof(int));
			}

#ifdef ARTIO_MPI
			if ( file_sfc_index[i+1] > handle->proc_sfc_end+1 ) {
				MPI_Send( &cur, 1, MPI_LONG_LONG_INT, handle->rank + 1, i, handle->context->comm );
			}
#endif /* ARTIO_MPI */

			/* seek and write our portion of sfc table */
			artio_file_fseek(ghandle->ffh[i], 
					(first_file_sfc - file_sfc_index[i]) * sizeof(int64_t), ARTIO_SEEK_SET);
			artio_file_fwrite(ghandle->ffh[i],
					&ghandle->sfc_offset_table[first_file_sfc - ghandle->cache_sfc_begin],
					last_file_sfc - first_file_sfc + 1, ARTIO_TYPE_LONG);
		}
	}

	free(file_parent);

	ghandle->cur_file = -1;
	ghandle->cur_sfc = -1;
	ghandle->cur_level = -1;
	ghandle->cur_num_levels = -1;
	ghandle->cur_octs = -1;

	handle->grid = ghandle;

	artio_parameter_set_long_array(handle, "grid_file_sfc_index",
			ghandle->num_grid_files + 1, ghandle->file_sfc_index);
	artio_parameter_set_int(handle, "grid_max_level", ghandle->file_max_level);
	
	return ARTIO_SUCCESS;
}

int artio_grid_close(artio_file handle) {
	int i;
	artio_grid_file ghandle = handle->grid;
	for (i = 0; i < ghandle->num_grid_files; i++) {
		artio_file_fclose(ghandle->ffh[i]);
	}
	free(ghandle->ffh);

	if ( ghandle->sfc_offset_table != NULL ) {
		free(ghandle->sfc_offset_table);
	}

	free(ghandle->octs_per_level);
	free(ghandle->file_sfc_index);
	return ARTIO_SUCCESS;
}

int artio_grid_cache_sfc_range(artio_file handle, int64_t start, int64_t end) {
	int i;
	int first_file, last_file, first, count, cur;
	artio_grid_file ghandle = handle->grid;

	if ( start > end || start < handle->proc_sfc_begin ||
			end > handle->proc_sfc_end ) {
		return ARTIO_ERR_INVALID_SFC_RANGE;
	}

	if (handle->open_mode != ARTIO_FILESET_READ || 
			!(handle->open_type & ARTIO_OPEN_GRID) ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	if (ghandle->sfc_offset_table != NULL) {
		free(ghandle->sfc_offset_table);
	}

	first_file = artio_grid_find_file(handle, 0, ghandle->num_grid_files, start);
	last_file = artio_grid_find_file(handle, first_file, ghandle->num_grid_files, end);

	ghandle->cache_sfc_begin = start;
	ghandle->cache_sfc_end = end;
	ghandle->sfc_offset_table = (int64_t *)malloc(sizeof(int64_t) * (end - start + 1));

	cur = 0;
	for (i = first_file; i <= last_file; i++) {
		first = MAX( 0, start - ghandle->file_sfc_index[i] );
		count = MIN( ghandle->file_sfc_index[i+1], end+1 )
				- MAX( start, ghandle->file_sfc_index[i]);
		artio_file_fseek(ghandle->ffh[i], 
				sizeof(int64_t) * first,
				ARTIO_SEEK_SET);
		artio_file_fread(ghandle->ffh[i], 
				&ghandle->sfc_offset_table[cur],
				count, ARTIO_TYPE_LONG);
		cur += count;
	}

	return ARTIO_SUCCESS;
}

int artio_grid_find_file(artio_file handle, int start, int end, int64_t sfc) {
	int j;
	artio_grid_file ghandle = handle->grid;

	if ( start < 0 || start > ghandle->num_grid_files ||
			end < 0 || end > ghandle->num_grid_files || 
			sfc < ghandle->file_sfc_index[start] ||
			sfc >= ghandle->file_sfc_index[end] ) {
		return -1;
	}

	if (start == end || sfc == ghandle->file_sfc_index[start]) {
		return start;
	}

	if (1 == end - start) {
		if (sfc < ghandle->file_sfc_index[end]) {
			return start;
		} else {
			return end;
		}
	}

	j = start + (end - start) / 2;
	if (sfc > ghandle->file_sfc_index[j]) {
		return artio_grid_find_file(handle, j, end, sfc);
	} else if (sfc < ghandle->file_sfc_index[j]) {
		return artio_grid_find_file(handle, start, j, sfc);
	} else {
		return j;
	}
}

int artio_grid_seek_to_sfc(artio_file handle, int64_t sfc) {
	int re;
	int64_t offset;
	artio_grid_file ghandle = handle->grid;

	if (ghandle->cache_sfc_begin == -1 || 
			sfc < ghandle->cache_sfc_begin || 
			sfc > ghandle->cache_sfc_end) {
		return ARTIO_ERR_INVALID_SFC;
	}

	ghandle->cur_file = artio_grid_find_file(handle, 0, ghandle->num_grid_files, sfc);
	offset = ghandle->sfc_offset_table[sfc - ghandle->cache_sfc_begin];
	re = artio_file_fseek(ghandle->ffh[ghandle->cur_file], offset, ARTIO_SEEK_SET);

	if ( re != ARTIO_SUCCESS ) {
		fprintf(stderr,"unable to seek to %ld in file %d for sfc %ld\n", 
			offset, ghandle->cur_file, sfc ); fflush(stderr);
		exit(1);
	}

	return ARTIO_SUCCESS;
}

int artio_grid_write_root_cell_begin(artio_file handle, int64_t sfc,
		float *variables, int num_oct_levels, int *num_octs_per_level) {
	int i;

	artio_grid_file ghandle = handle->grid;

	if (handle->open_mode != ARTIO_FILESET_WRITE || 
			!(handle->open_type & ARTIO_OPEN_GRID) ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	if (num_oct_levels < 0 || num_oct_levels > ghandle->file_max_level) {
		return ARTIO_ERR_INVALID_OCT_LEVELS;
	}	

	artio_grid_seek_to_sfc(handle, sfc);
	artio_file_fwrite(ghandle->ffh[ghandle->cur_file], variables, ghandle->num_grid_variables,
			ARTIO_TYPE_FLOAT);
	artio_file_fwrite(ghandle->ffh[ghandle->cur_file], &num_oct_levels, 1, ARTIO_TYPE_INT);
	artio_file_fwrite(ghandle->ffh[ghandle->cur_file], num_octs_per_level, num_oct_levels,
			ARTIO_TYPE_INT);

	for (i = 0; i < num_oct_levels; i++) {
		ghandle->octs_per_level[i] = num_octs_per_level[i];
	}

	ghandle->cur_sfc = sfc;
	ghandle->cur_num_levels = num_oct_levels;
	ghandle->cur_level = -1;
	ghandle->cur_octs = 0;

	return ARTIO_SUCCESS;
}

int artio_grid_write_root_cell_end(artio_file handle) {
	artio_grid_file ghandle = handle->grid;

	if (handle->open_mode != ARTIO_FILESET_WRITE || 
			!(handle->open_type & ARTIO_OPEN_GRID) ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	ghandle->cur_sfc = -1;
	return ARTIO_SUCCESS;
}

int artio_grid_write_level_begin(artio_file handle, int level) {
	artio_grid_file ghandle = handle->grid;

	if (handle->open_mode != ARTIO_FILESET_WRITE || 
			!(handle->open_type & ARTIO_OPEN_GRID) ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	if (ghandle->cur_sfc == -1 ||
			level <= 0 || level > ghandle->cur_num_levels) {
		return ARTIO_ERR_INVALID_STATE;
	}

	ghandle->cur_level = level;
	return ARTIO_SUCCESS;
}

int artio_grid_write_level_end(artio_file handle) {
	artio_grid_file ghandle = handle->grid;

	if (handle->open_mode != ARTIO_FILESET_WRITE || 
			!(handle->open_type & ARTIO_OPEN_GRID) ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	if (ghandle->cur_level == -1 || 
			ghandle->cur_octs != ghandle->octs_per_level[ghandle->cur_level - 1] ) {
		return ARTIO_ERR_INVALID_STATE;
	}

	ghandle->cur_level = -1;
	ghandle->cur_octs = 0;

	return ARTIO_SUCCESS;
}

int artio_grid_write_oct(artio_file handle, float *variables,
		int *cellrefined) {
	artio_grid_file ghandle = handle->grid;

	if (handle->open_mode != ARTIO_FILESET_WRITE || 
			!(handle->open_type & ARTIO_OPEN_GRID) ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	if (ghandle->cur_level == -1 || 
			ghandle->cur_octs >= ghandle->octs_per_level[ghandle->cur_level - 1]) {
		return ARTIO_ERR_INVALID_STATE;
	}

	artio_file_fwrite(ghandle->ffh[ghandle->cur_file], 
			variables, 8 * ghandle->num_grid_variables,
			ARTIO_TYPE_FLOAT);
	artio_file_fwrite(ghandle->ffh[ghandle->cur_file], 
			cellrefined, 8, ARTIO_TYPE_INT);
	ghandle->cur_octs++;
	return ARTIO_SUCCESS;
}

/*
 *
 */
int artio_grid_read_root_cell_begin(artio_file handle, int64_t sfc,
		float *variables, int *num_oct_levels, int *num_octs_per_level) {
	int i;
	artio_grid_file ghandle = handle->grid;

	if (handle->open_mode != ARTIO_FILESET_READ || 
			!(handle->open_type & ARTIO_OPEN_GRID) ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	artio_grid_seek_to_sfc(handle, sfc);
	artio_file_fread(ghandle->ffh[ghandle->cur_file], 
			variables, ghandle->num_grid_variables,
			ARTIO_TYPE_FLOAT);
	artio_file_fread(ghandle->ffh[ghandle->cur_file], 
			num_oct_levels, 1, ARTIO_TYPE_INT);
	if ( *num_oct_levels > ghandle->file_max_level ) {
		fprintf(stderr,"%d: error reading root num_oct_levels\n", handle->rank );
		exit(1);
	}

	if (*num_oct_levels > 0) {
		artio_file_fread(ghandle->ffh[ghandle->cur_file], 
				num_octs_per_level, *num_oct_levels,
				ARTIO_TYPE_INT);
		for (i = 0; i < *num_oct_levels; i++) {
			ghandle->octs_per_level[i] = num_octs_per_level[i];
		}
	}

	ghandle->cur_sfc = sfc;
	ghandle->cur_num_levels = *num_oct_levels;
	ghandle->cur_level = -1;

	return ARTIO_SUCCESS;
}

/* Description  */
int artio_grid_read_oct(artio_file handle, float *variables,
		int *refined) {
	artio_grid_file ghandle = handle->grid;

	if (handle->open_mode != ARTIO_FILESET_READ ||                                                                                    
			!(handle->open_type & ARTIO_OPEN_GRID) ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	if (ghandle->cur_level == -1 || 
			ghandle->cur_octs > ghandle->octs_per_level[ghandle->cur_level - 1]) {
		return ARTIO_ERR_INVALID_STATE;
	}

	artio_file_fread(ghandle->ffh[ghandle->cur_file], 
			variables, 8 * ghandle->num_grid_variables,
			ARTIO_TYPE_FLOAT);
	artio_file_fread(ghandle->ffh[ghandle->cur_file], 
			refined, 8, ARTIO_TYPE_INT);
	ghandle->cur_octs++;

	return ARTIO_SUCCESS;
}

/*
 * Description        Obtain data from an appointed level tree node
 */
int artio_grid_read_level_begin(artio_file handle, int level) {
	int i;
	int64_t offset = 0;
	artio_grid_file ghandle = handle->grid;

	if (handle->open_mode != ARTIO_FILESET_READ ||                                                                                    
			!(handle->open_type & ARTIO_OPEN_GRID) ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	if ( ghandle->cur_sfc == -1 || level <= 0 || level > ghandle->cur_num_levels ) {
		return ARTIO_ERR_INVALID_STATE;
	}

	offset = ghandle->sfc_offset_table[ghandle->cur_sfc - ghandle->cache_sfc_begin];
	offset += sizeof(float) * ghandle->num_grid_variables + sizeof(int)
			* (ghandle->cur_num_levels + 1);
	for (i = 0; i < level - 1; i++) {
		offset += 8 * (sizeof(float) * ghandle->num_grid_variables + sizeof(int))
				* ghandle->octs_per_level[i];
	}

	artio_file_fseek(ghandle->ffh[ghandle->cur_file], offset, ARTIO_SEEK_SET);
	ghandle->cur_level = level;
	ghandle->cur_octs = 0;

	return ARTIO_SUCCESS;
}

/*
 * Description        Do something at the end of each kind of read operation
 */
int artio_grid_read_level_end(artio_file handle) {
	artio_grid_file ghandle = handle->grid;

	if (handle->open_mode != ARTIO_FILESET_READ ||                                                                                    
			!(handle->open_type & ARTIO_OPEN_GRID) ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	if (ghandle->cur_level == -1) {
		return ARTIO_ERR_INVALID_STATE;
	}

	ghandle->cur_level = -1;
	ghandle->cur_octs = -1;

	return ARTIO_SUCCESS;
}

int artio_grid_read_root_cell_end(artio_file handle) {
	if (handle->open_mode != ARTIO_FILESET_READ ||                                                                                    
			!(handle->open_type & ARTIO_OPEN_GRID) ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}
	handle->grid->cur_sfc = -1;
	return ARTIO_SUCCESS;
}

int artio_grid_read_sfc_range(artio_file handle, 
		int64_t sfc1, int64_t sfc2,
		int min_level_to_read, int max_level_to_read, 
		int options,
		GridCallBack callback) {
	int64_t sfc;
	int oct, level, j;
	int *octs_per_level;
	int refined;
	int oct_refined[8];
	int root_tree_levels;
	float * variable_of_root;
	float * variable_of_oct = 0;

	artio_grid_file ghandle = handle->grid;

	if (handle->open_mode != ARTIO_FILESET_READ ||                                                                                    
			!(handle->open_type & ARTIO_OPEN_GRID) ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	if ((min_level_to_read < 0) | (min_level_to_read > max_level_to_read)) {
		return ARTIO_ERR_INVALID_LEVEL;
	}

	octs_per_level = (int *)malloc(ghandle->file_max_level * sizeof(int));
	variable_of_root = (float *)malloc(ghandle->num_grid_variables * sizeof(float));

	if (max_level_to_read > 0) {
		variable_of_oct = (float *)malloc(8 * ghandle->num_grid_variables * sizeof(float));
	}

	artio_grid_cache_sfc_range(handle, sfc1, sfc2);

	for (sfc = sfc1; sfc <= sfc2; sfc++) {
		artio_grid_read_root_cell_begin(handle, sfc, variable_of_root,
				&root_tree_levels, octs_per_level);

		if (min_level_to_read == 0 && (options == ARTIO_READ_ALL || 
				options == ARTIO_READ_REFINED && root_tree_levels > 0 || 
				options == ARTIO_READ_LEAFS && root_tree_levels == 0)) {
			refined = (root_tree_levels > 0) ? 1 : 0;
			callback(variable_of_root, 0, refined, sfc);
		}

		for (level = MAX(min_level_to_read,1); 
				level <= MIN(root_tree_levels,max_level_to_read); level++) {
			artio_grid_read_level_begin(handle, level);
			for (oct = 0; oct < octs_per_level[level - 1]; oct++) {
				artio_grid_read_oct(handle, variable_of_oct, oct_refined);

				for (j = 0; j < 8; j++) {
					if (options == ARTIO_READ_ALL || 
							options == ARTIO_READ_REFINED && oct_refined[j] ||
							options == ARTIO_READ_LEAFS && !oct_refined[j] ) {
						callback(&variable_of_oct[j * ghandle->num_grid_variables],
								level, oct_refined[j], sfc);
					}
				}
			}
			artio_grid_read_level_end(handle);
		}
		artio_grid_read_root_cell_end(handle);
	}

	if (max_level_to_read > 0) {
		free(variable_of_oct);
	}

	free(variable_of_root);
	free(octs_per_level);

	return ARTIO_SUCCESS;
}
