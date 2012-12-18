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

#include "sfc.h"
#include "artio.h"
#include "artio_internal.h"

int artio_grid_find_file(artio_grid_file ghandle, int start, int end, int64_t sfc);
artio_grid_file artio_grid_file_allocate(void);
void artio_grid_file_destroy(artio_grid_file ghandle);

/*
 * Open grid component of the fileset
 */
int artio_fileset_open_grid(artio_file handle) {
	int i;
	char filename[256];
	int first_file, last_file;
	int mode;
	artio_grid_file ghandle;

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	/* check that the fileset doesn't already contain a grid component */
	if ( handle->open_type & ARTIO_OPEN_GRID ||
			handle->open_mode != ARTIO_FILESET_READ ||
			handle->grid != NULL ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}
	handle->open_type |= ARTIO_OPEN_GRID;

	ghandle = artio_grid_file_allocate();
	if ( ghandle == NULL ) {
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}

	/* load grid parameters from header file (should be doing error handling...) */
	artio_parameter_get_int(handle, "num_grid_variables", &ghandle->num_grid_variables);
	artio_parameter_get_int(handle, "num_grid_files", &ghandle->num_grid_files);

	ghandle->file_sfc_index = (int64_t *)malloc(sizeof(int64_t) * (ghandle->num_grid_files + 1));
	if ( ghandle->file_sfc_index == NULL ) {
		artio_grid_file_destroy(ghandle);
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}

	artio_parameter_get_long_array(handle, "grid_file_sfc_index",
				       ghandle->num_grid_files + 1, ghandle->file_sfc_index);
	artio_parameter_get_int(handle, "grid_max_level",
				&ghandle->file_max_level);

	ghandle->octs_per_level = (int *)malloc(ghandle->file_max_level * sizeof(int));
	if ( ghandle->octs_per_level == NULL ) {
		artio_grid_file_destroy(ghandle);
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}

	ghandle->ffh = (artio_fh *)malloc(ghandle->num_grid_files * sizeof(artio_fh));
	if ( ghandle->ffh == NULL ) {
		artio_grid_file_destroy(ghandle);
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}
 
	for ( i = 0; i < ghandle->num_grid_files; i++ ) {
		ghandle->ffh[i] = NULL;
	}

	first_file = artio_grid_find_file(ghandle, 0, 
			ghandle->num_grid_files, handle->proc_sfc_begin);
	last_file = artio_grid_find_file(ghandle, first_file,
			ghandle->num_grid_files, handle->proc_sfc_end);

	/* open files on all processes */
	for (i = 0; i < ghandle->num_grid_files; i++) {
		sprintf(filename, "%s.g%03d", handle->file_prefix, i);

		mode = ARTIO_MODE_READ;
		if (i >= first_file && i <= last_file) {
			mode |= ARTIO_MODE_ACCESS;
		}

		if (handle->endian_swap) {
			mode |= ARTIO_MODE_ENDIAN_SWAP;
		}

		ghandle->ffh[i] = artio_file_fopen(filename, mode, handle->context);
		if ( ghandle->ffh[i] == NULL ) {
			artio_grid_file_destroy(ghandle);
			return ARTIO_ERR_GRID_FILE_NOT_FOUND;
		}
	}

	handle->grid = ghandle;
	return ARTIO_SUCCESS;
}

int artio_fileset_add_grid(artio_file handle, 
		int num_grid_files, int allocation_strategy, 
		int num_grid_variables, 
		char ** grid_variable_labels,
		int * num_levels_per_root_tree,
		int * num_octs_per_root_tree ) {

	int i;
	int file_max_level, local_max_level;
	int64_t cur, sfc, l;
	int64_t first_file_sfc, last_file_sfc;
	int first_file, last_file;
	char filename[256];
	int mode;
	int ret;
	artio_grid_file ghandle;

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

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

	ghandle = artio_grid_file_allocate();
	if ( ghandle == NULL ) {
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}

	ghandle->file_sfc_index = (int64_t *)malloc(sizeof(int64_t) * (num_grid_files + 1));
	if ( ghandle->file_sfc_index == NULL ) {
		artio_grid_file_destroy(ghandle);
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}

	/* compute global maximum level */
	local_max_level = 0;
	for (sfc = 0; sfc < handle->proc_sfc_end - handle->proc_sfc_begin + 1; sfc++) {
		if (num_levels_per_root_tree[sfc] > local_max_level) {
			local_max_level = num_levels_per_root_tree[sfc];
		}
	}

#ifdef ARTIO_MPI
	MPI_Allreduce( &local_max_level, &file_max_level, 
			1, MPI_INT, MPI_MAX, handle->context->comm );
#else
	file_max_level = local_max_level;
#endif /* ARTIO_MPI */

	switch (allocation_strategy) {
		case ARTIO_ALLOC_EQUAL_PROC:
			if (num_grid_files > handle->num_procs) {
				return ARTIO_ERR_INVALID_FILE_NUMBER;
			}

			for (i = 0; i < num_grid_files; i++) {
				ghandle->file_sfc_index[i] = 
					handle->proc_sfc_index[(handle->num_procs*i+num_grid_files-1) / num_grid_files];
			}
			ghandle->file_sfc_index[num_grid_files] = 
				handle->proc_sfc_index[handle->num_procs];
			break;
		case ARTIO_ALLOC_EQUAL_SFC:
			if ( num_grid_files > handle->num_root_cells ) {
				return ARTIO_ERR_INVALID_FILE_NUMBER;
			}

			for (i = 0; i < num_grid_files; i++) {
				ghandle->file_sfc_index[i] = 
					(handle->num_root_cells*i+num_grid_files-1) / num_grid_files;
			}
			ghandle->file_sfc_index[num_grid_files] = handle->num_root_cells;
			break;
		default:
			artio_grid_file_destroy(ghandle);
			return ARTIO_ERR_INVALID_ALLOC_STRATEGY;
	}

	ghandle->num_grid_files = num_grid_files;
	ghandle->num_grid_variables = num_grid_variables;
	ghandle->file_max_level = file_max_level;

	/* allocate space for sfc offset cache */
	ghandle->cache_sfc_begin = handle->proc_sfc_begin;
	ghandle->cache_sfc_end = handle->proc_sfc_end;
	ghandle->sfc_offset_table = 
		(int64_t *)malloc((size_t)(ghandle->cache_sfc_end - 
					ghandle->cache_sfc_begin + 1) * sizeof(int64_t));
	if ( ghandle->sfc_offset_table == NULL ) {
		artio_grid_file_destroy(ghandle);
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}

	ghandle->octs_per_level = (int *)malloc(ghandle->file_max_level * sizeof(int));
	if ( ghandle->octs_per_level == NULL ) {
		artio_grid_file_destroy(ghandle);
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}

	/* allocate file handles */
	ghandle->ffh = (artio_fh *)malloc(num_grid_files * sizeof(artio_fh));
	if ( ghandle->ffh == NULL ) {
		artio_grid_file_destroy(ghandle);
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}
	for ( i = 0; i < num_grid_files; i++ ) {
		ghandle->ffh[i] = NULL;
	}

	/* open file handles */
	first_file = artio_grid_find_file(ghandle, 0, num_grid_files, 
					handle->proc_sfc_begin);
	last_file = artio_grid_find_file(ghandle, first_file, num_grid_files, 
					handle->proc_sfc_end);

	if ( first_file < 0 || first_file >= num_grid_files ||
			last_file < first_file || last_file >= num_grid_files ) {
		return ARTIO_ERR_INVALID_FILE_NUMBER;
	}

	for (i = 0; i < num_grid_files; i++) {
		sprintf(filename, "%s.g%03d", handle->file_prefix, i);

		mode = ARTIO_MODE_WRITE;
		if (i >= first_file && i <= last_file) {
			mode |= ARTIO_MODE_ACCESS;
		}

		ghandle->ffh[i] = artio_file_fopen(filename, mode, handle->context);
		if ( ghandle->ffh[i] == NULL ) {
			artio_grid_file_destroy(ghandle);
			return ARTIO_ERR_FILE_CREATE;
		}

		/* write sfc offset header if we contribute to this file */
		if (i >= first_file && i <= last_file) {
#ifdef ARTIO_MPI
			if (ghandle->file_sfc_index[i] >= handle->proc_sfc_index[ handle->rank ] &&
					ghandle->file_sfc_index[i] < handle->proc_sfc_index[ handle->rank + 1] ) {
				cur = (ghandle->file_sfc_index[i + 1] - ghandle->file_sfc_index[i]) * sizeof(int64_t);
			} else {
				/* obtain offset from previous process */
				MPI_Recv( &cur, 1, MPI_LONG_LONG_INT, handle->rank - 1, i,
						handle->context->comm, MPI_STATUS_IGNORE );
			}
#else
			cur = (ghandle->file_sfc_index[i + 1] - ghandle->file_sfc_index[i]) * sizeof(int64_t);
#endif /* ARTIO_MPI */

			first_file_sfc = MAX( handle->proc_sfc_begin, ghandle->file_sfc_index[i] );
			last_file_sfc = MIN( handle->proc_sfc_end, ghandle->file_sfc_index[i+1]-1 );

			for (l = first_file_sfc - ghandle->cache_sfc_begin; 
					l < last_file_sfc - ghandle->cache_sfc_begin + 1; l++) {
				ghandle->sfc_offset_table[l] = cur;
				cur += sizeof(float) * ghandle->num_grid_variables + sizeof(int) * (1
						+ num_levels_per_root_tree[l])
						+ num_octs_per_root_tree[l] * 8 * (sizeof(float)
								* ghandle->num_grid_variables + sizeof(int));
			}

#ifdef ARTIO_MPI
			if ( ghandle->file_sfc_index[i+1] > handle->proc_sfc_end+1 ) {
				MPI_Send( &cur, 1, MPI_LONG_LONG_INT, handle->rank + 1, i, handle->context->comm );
			}
#endif /* ARTIO_MPI */

			/* seek and write our portion of sfc table */
			ret = artio_file_fseek(ghandle->ffh[i], 
					(first_file_sfc - ghandle->file_sfc_index[i]) * sizeof(int64_t), 
					ARTIO_SEEK_SET);
			if ( ret != ARTIO_SUCCESS ) {
				artio_grid_file_destroy(ghandle);
				return ret;
			}

			ret = artio_file_fwrite(ghandle->ffh[i],
					&ghandle->sfc_offset_table[first_file_sfc - ghandle->cache_sfc_begin],
					last_file_sfc - first_file_sfc + 1, ARTIO_TYPE_LONG);
			if ( ret != ARTIO_SUCCESS ) {
				artio_grid_file_destroy(ghandle);
				return ret;
			}
		}
	}

	handle->grid = ghandle;

	artio_parameter_set_long_array(handle, "grid_file_sfc_index",
			ghandle->num_grid_files + 1, ghandle->file_sfc_index);
	artio_parameter_set_int(handle, "grid_max_level", ghandle->file_max_level);
	
	return ARTIO_SUCCESS;
}

artio_grid_file artio_grid_file_allocate(void) {
	artio_grid_file ghandle = 
			(artio_grid_file)malloc(sizeof(struct artio_grid_file_struct));                                          
	if ( ghandle != NULL ) {
		ghandle->ffh = NULL;
		ghandle->num_grid_variables = -1;
		ghandle->num_grid_files = -1;
		ghandle->file_sfc_index = NULL;
		ghandle->cache_sfc_begin = -1;
		ghandle->cache_sfc_end = -1;
		ghandle->sfc_offset_table = NULL;
		ghandle->file_max_level = -1;
		ghandle->cur_file = -1;
		ghandle->cur_num_levels = -1;
		ghandle->cur_level = -1;
		ghandle->cur_octs = -1;
		ghandle->cur_sfc = -1;
		ghandle->octs_per_level = NULL;
    }
	return ghandle;
}

void artio_grid_file_destroy(artio_grid_file ghandle) {
	int i;
	if ( ghandle == NULL ) return;	

	if ( ghandle->ffh != NULL ) {
		for (i = 0; i < ghandle->num_grid_files; i++) {
			if ( ghandle->ffh[i] != NULL ) {
				artio_file_fclose(ghandle->ffh[i]);
			}
		}
		free(ghandle->ffh);
	}

	if ( ghandle->sfc_offset_table != NULL ) free(ghandle->sfc_offset_table);
	if ( ghandle->octs_per_level != NULL ) free(ghandle->octs_per_level);
	if ( ghandle->file_sfc_index != NULL ) free(ghandle->file_sfc_index);
	free(ghandle);
}

int artio_fileset_close_grid(artio_file handle) {
	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if ( !(handle->open_type & ARTIO_OPEN_GRID) || 
			handle->grid == NULL ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	artio_grid_file_destroy(handle->grid);
	handle->grid = NULL;
	return ARTIO_SUCCESS;
}

int artio_grid_cache_sfc_range(artio_file handle, int64_t start, int64_t end) {
	int i;
	int ret;
	int first_file, last_file;
	int64_t first, count, cur;
	artio_grid_file ghandle;

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if (handle->open_mode != ARTIO_FILESET_READ || 
			!(handle->open_type & ARTIO_OPEN_GRID) ||
			handle->grid == NULL ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	ghandle = handle->grid;

	if ( start > end || start < handle->proc_sfc_begin ||
			end > handle->proc_sfc_end ) {
		return ARTIO_ERR_INVALID_SFC_RANGE;
	}


	if (ghandle->sfc_offset_table != NULL) {
		free(ghandle->sfc_offset_table);
	}

	first_file = artio_grid_find_file(ghandle, 0, ghandle->num_grid_files, start);
	last_file = artio_grid_find_file(ghandle, first_file, ghandle->num_grid_files, end);

	ghandle->cache_sfc_begin = start;
	ghandle->cache_sfc_end = end;
	ghandle->sfc_offset_table = (int64_t *)malloc(sizeof(int64_t) * (size_t)(end - start + 1));
	if ( ghandle->sfc_offset_table == NULL ) {
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}

	cur = 0;
	for (i = first_file; i <= last_file; i++) {
		first = MAX( 0, start - ghandle->file_sfc_index[i] );
		count = MIN( ghandle->file_sfc_index[i+1], end+1 )
				- MAX( start, ghandle->file_sfc_index[i]);
		ret = artio_file_fseek(ghandle->ffh[i], 
				sizeof(int64_t) * first, ARTIO_SEEK_SET);
		if ( ret != ARTIO_SUCCESS ) return ret;

		ret = artio_file_fread(ghandle->ffh[i], 
				&ghandle->sfc_offset_table[cur],
				count, ARTIO_TYPE_LONG);
		if ( ret != ARTIO_SUCCESS ) return ret;

		cur += count;
	}

	return ARTIO_SUCCESS;
}

int artio_grid_find_file(artio_grid_file ghandle, int start, int end, int64_t sfc) {
	int j;

	if ( ghandle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

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
		return artio_grid_find_file(ghandle, j, end, sfc);
	} else if (sfc < ghandle->file_sfc_index[j]) {
		return artio_grid_find_file(ghandle, start, j, sfc);
	} else {
		return j;
	}
}

int artio_grid_seek_to_sfc(artio_file handle, int64_t sfc) {
	int64_t offset;
	artio_grid_file ghandle;

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if ( !(handle->open_type & ARTIO_OPEN_GRID) ||
			handle->grid == NULL ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	ghandle = handle->grid;

	if (ghandle->cache_sfc_begin == -1 || 
			sfc < ghandle->cache_sfc_begin || 
			sfc > ghandle->cache_sfc_end) {
		return ARTIO_ERR_INVALID_SFC;
	}

	ghandle->cur_file = artio_grid_find_file(ghandle, 0, ghandle->num_grid_files, sfc);
	offset = ghandle->sfc_offset_table[sfc - ghandle->cache_sfc_begin];
	return artio_file_fseek(ghandle->ffh[ghandle->cur_file], 
			offset, ARTIO_SEEK_SET);
}

int artio_grid_write_root_cell_begin(artio_file handle, int64_t sfc,
		float *variables, int num_oct_levels, int *num_octs_per_level) {
	int i;
	int ret;
	artio_grid_file ghandle;

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if (handle->open_mode != ARTIO_FILESET_WRITE || 
			!(handle->open_type & ARTIO_OPEN_GRID) ||
			handle->grid == NULL ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	ghandle = handle->grid;

	if (num_oct_levels < 0 || num_oct_levels > ghandle->file_max_level) {
		return ARTIO_ERR_INVALID_OCT_LEVELS;
	}	

	ret = artio_grid_seek_to_sfc(handle, sfc);
	if ( ret != ARTIO_SUCCESS ) return ret;

	ret = artio_file_fwrite(ghandle->ffh[ghandle->cur_file], variables, 
			ghandle->num_grid_variables, ARTIO_TYPE_FLOAT);
	if ( ret != ARTIO_SUCCESS ) return ret;

	ret = artio_file_fwrite(ghandle->ffh[ghandle->cur_file], 
			&num_oct_levels, 1, ARTIO_TYPE_INT);
	if ( ret != ARTIO_SUCCESS ) return ret;

	ret = artio_file_fwrite(ghandle->ffh[ghandle->cur_file], 
			num_octs_per_level, num_oct_levels, ARTIO_TYPE_INT);
	if ( ret != ARTIO_SUCCESS ) return ret;

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
	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if (handle->open_mode != ARTIO_FILESET_WRITE || 
			!(handle->open_type & ARTIO_OPEN_GRID) ||
			handle->grid == NULL ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	handle->grid->cur_sfc = -1;
	return ARTIO_SUCCESS;
}

int artio_grid_write_level_begin(artio_file handle, int level) {
	artio_grid_file ghandle;

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if (handle->open_mode != ARTIO_FILESET_WRITE || 
			!(handle->open_type & ARTIO_OPEN_GRID) ||
			handle->grid == NULL ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	ghandle = handle->grid;

	if (ghandle->cur_sfc == -1 ||
			level <= 0 || level > ghandle->cur_num_levels) {
		return ARTIO_ERR_INVALID_STATE;
	}

	ghandle->cur_level = level;
	return ARTIO_SUCCESS;
}

int artio_grid_write_level_end(artio_file handle) {
	artio_grid_file ghandle;

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if (handle->open_mode != ARTIO_FILESET_WRITE || 
			!(handle->open_type & ARTIO_OPEN_GRID) ||
			handle->grid == NULL ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	ghandle = handle->grid;

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
	int i;
	int ret;
	artio_grid_file ghandle;

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if (handle->open_mode != ARTIO_FILESET_WRITE || 
			!(handle->open_type & ARTIO_OPEN_GRID) ||
			handle->grid == NULL ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	ghandle = handle->grid;

	if (ghandle->cur_level == -1 || 
			ghandle->cur_octs >= ghandle->octs_per_level[ghandle->cur_level - 1]) {
		return ARTIO_ERR_INVALID_STATE;
	}

	/* check that no last-level octs have refined cells */
	if ( ghandle->cur_level == ghandle->cur_num_levels ) {
		for ( i = 0; i < 8; i++ ) {
			if ( cellrefined[i] ) {
				return ARTIO_ERR_INVALID_OCT_REFINED;
			}
		}
	}

	ret = artio_file_fwrite(ghandle->ffh[ghandle->cur_file], 
			variables, 8 * ghandle->num_grid_variables,
			ARTIO_TYPE_FLOAT);
	if ( ret != ARTIO_SUCCESS ) return ret;

	ret = artio_file_fwrite(ghandle->ffh[ghandle->cur_file], 
			cellrefined, 8, ARTIO_TYPE_INT);
	if ( ret != ARTIO_SUCCESS ) return ret;

	ghandle->cur_octs++;
	return ARTIO_SUCCESS;
}

/*
 *
 */
int artio_grid_read_root_nocts(artio_file handle, int64_t sfc,
		float *variables, int *num_oct_levels, int *num_octs_per_level) {
	int i;
	int ret;
	artio_grid_file ghandle;

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if (handle->open_mode != ARTIO_FILESET_READ || 
			!(handle->open_type & ARTIO_OPEN_GRID) ||
			handle->grid == NULL ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	ghandle = handle->grid;

	ret = artio_grid_seek_to_sfc(handle, sfc);
	if ( ret != ARTIO_SUCCESS ) return ret;

	ret = artio_file_fread(ghandle->ffh[ghandle->cur_file], 
			variables, ghandle->num_grid_variables,
			ARTIO_TYPE_FLOAT);
	if ( ret != ARTIO_SUCCESS ) return ret;

	ret = artio_file_fread(ghandle->ffh[ghandle->cur_file], 
			num_oct_levels, 1, ARTIO_TYPE_INT);
	if ( ret != ARTIO_SUCCESS ) return ret;

	if ( *num_oct_levels > ghandle->file_max_level ) {
		return ARTIO_ERR_INVALID_OCT_LEVELS;
	}

	if (*num_oct_levels > 0) {
		ret = artio_file_fread(ghandle->ffh[ghandle->cur_file], 
				num_octs_per_level, *num_oct_levels,
				ARTIO_TYPE_INT);
		if ( ret != ARTIO_SUCCESS ) return ret;

		for (i = 0; i < *num_oct_levels; i++) {
			ghandle->octs_per_level[i] = num_octs_per_level[i];
		}
	}

	ghandle->cur_sfc = sfc;
	ghandle->cur_num_levels = *num_oct_levels;
	ghandle->cur_level = -1;

	return ARTIO_SUCCESS;
}
int artio_grid_read_root_cell_begin(artio_file handle, int64_t sfc,
		float *variables, int *num_oct_levels, int *num_octs_per_level) {
	int i;
	int ret;
	artio_grid_file ghandle;

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if (handle->open_mode != ARTIO_FILESET_READ || 
			!(handle->open_type & ARTIO_OPEN_GRID) ||
			handle->grid == NULL ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	ghandle = handle->grid;

	ret = artio_grid_seek_to_sfc(handle, sfc);
	if ( ret != ARTIO_SUCCESS ) return ret;

	ret = artio_file_fread(ghandle->ffh[ghandle->cur_file], 
			variables, ghandle->num_grid_variables,
			ARTIO_TYPE_FLOAT);
	if ( ret != ARTIO_SUCCESS ) return ret;

	ret = artio_file_fread(ghandle->ffh[ghandle->cur_file], 
			num_oct_levels, 1, ARTIO_TYPE_INT);
	if ( ret != ARTIO_SUCCESS ) return ret;

	if ( *num_oct_levels > ghandle->file_max_level ) {
		return ARTIO_ERR_INVALID_OCT_LEVELS;
	}

	if (*num_oct_levels > 0) {
		ret = artio_file_fread(ghandle->ffh[ghandle->cur_file], 
				num_octs_per_level, *num_oct_levels,
				ARTIO_TYPE_INT);
		if ( ret != ARTIO_SUCCESS ) return ret;

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
	int ret;
	artio_grid_file ghandle;

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if (handle->open_mode != ARTIO_FILESET_READ || 
			!(handle->open_type & ARTIO_OPEN_GRID) ||
			handle->grid == NULL ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	ghandle = handle->grid;

	if (ghandle->cur_level == -1 || 
			ghandle->cur_octs > ghandle->octs_per_level[ghandle->cur_level - 1]) {
		return ARTIO_ERR_INVALID_STATE;
	}

	ret = artio_file_fread(ghandle->ffh[ghandle->cur_file], 
			variables, 8 * ghandle->num_grid_variables,
			ARTIO_TYPE_FLOAT);
	if ( ret != ARTIO_SUCCESS ) return ret;

	ret = artio_file_fread(ghandle->ffh[ghandle->cur_file], 
			refined, 8, ARTIO_TYPE_INT);
	if ( ret != ARTIO_SUCCESS ) return ret;

	ghandle->cur_octs++;

	return ARTIO_SUCCESS;
}

/*
 * Description        Obtain data from an appointed level tree node
 */
int artio_grid_read_level_begin(artio_file handle, int level) {
	int i;
	int ret;
	int64_t offset = 0;
	artio_grid_file ghandle;

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if (handle->open_mode != ARTIO_FILESET_READ ||
			!(handle->open_type & ARTIO_OPEN_GRID) ||
			handle->grid == NULL ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	ghandle = handle->grid;

	if ( ghandle->cur_sfc == -1 || level <= 0 || 
			level > ghandle->cur_num_levels ) {
		return ARTIO_ERR_INVALID_STATE;
	}

	offset = ghandle->sfc_offset_table[ghandle->cur_sfc - ghandle->cache_sfc_begin];
	offset += sizeof(float) * ghandle->num_grid_variables + sizeof(int)
			* (ghandle->cur_num_levels + 1);
	for (i = 0; i < level - 1; i++) {
		offset += 8 * (sizeof(float) * ghandle->num_grid_variables + sizeof(int))
				* ghandle->octs_per_level[i];
	}

	ret = artio_file_fseek(ghandle->ffh[ghandle->cur_file], 
			offset, ARTIO_SEEK_SET);
	if ( ret != ARTIO_SUCCESS ) return ret;

	ghandle->cur_level = level;
	ghandle->cur_octs = 0;

	return ARTIO_SUCCESS;
}

/*
 * Description        Do something at the end of each kind of read operation
 */
int artio_grid_read_level_end(artio_file handle) {
	artio_grid_file ghandle;

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if (handle->open_mode != ARTIO_FILESET_READ ||
			!(handle->open_type & ARTIO_OPEN_GRID) ||
			handle->grid == NULL ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	ghandle = handle->grid;

	if (ghandle->cur_level == -1) {
		return ARTIO_ERR_INVALID_STATE;
	}

	ghandle->cur_level = -1;
	ghandle->cur_octs = -1;

	return ARTIO_SUCCESS;
}

int artio_grid_read_root_cell_end(artio_file handle) {
	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if (handle->open_mode != ARTIO_FILESET_READ ||
			!(handle->open_type & ARTIO_OPEN_GRID) ||
			handle->grid == NULL ) {
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
	int ret;
	int *octs_per_level = NULL;
	int refined;
	int oct_refined[8];
	int root_tree_levels;
	float * variables = NULL;

	artio_grid_file ghandle;

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if (handle->open_mode != ARTIO_FILESET_READ ||
			!(handle->open_type & ARTIO_OPEN_GRID) ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	ghandle = handle->grid;

	if ((min_level_to_read < 0) || (min_level_to_read > max_level_to_read)) {
		return ARTIO_ERR_INVALID_LEVEL;
	}

	octs_per_level = (int *)malloc(ghandle->file_max_level * sizeof(int));
	variables = (float *)malloc(8*ghandle->num_grid_variables * sizeof(float));

	if ( octs_per_level == NULL || variables == NULL ) {
		if ( octs_per_level != NULL ) free(octs_per_level);
		if ( variables != NULL ) free(variables);
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}

	ret = artio_grid_cache_sfc_range(handle, sfc1, sfc2);
	if ( ret != ARTIO_SUCCESS ) {
		free(octs_per_level);
		free(variables);
		return ret;
	}
		
	for (sfc = sfc1; sfc <= sfc2; sfc++) {
		ret = artio_grid_read_root_cell_begin(handle, sfc, variables,
				&root_tree_levels, octs_per_level);
		if ( ret != ARTIO_SUCCESS ) {
			free(octs_per_level);
			free(variables);
			return ret;
		}

		if (min_level_to_read == 0 && (options == ARTIO_READ_ALL || 
				(options == ARTIO_READ_REFINED && root_tree_levels > 0) || 
				(options == ARTIO_READ_LEAFS && root_tree_levels == 0)) ) {
			refined = (root_tree_levels > 0) ? 1 : 0;
			callback(variables, 0, refined, sfc);
		}

		for (level = MAX(min_level_to_read,1); 
				level <= MIN(root_tree_levels,max_level_to_read); level++) {
			ret = artio_grid_read_level_begin(handle, level);
			if ( ret != ARTIO_SUCCESS ) {
				free(octs_per_level);
				free(variables);
				return ret;
			}

			for (oct = 0; oct < octs_per_level[level - 1]; oct++) {
				ret = artio_grid_read_oct(handle, variables, oct_refined);
				if ( ret != ARTIO_SUCCESS ) {
					free(octs_per_level);
					free(variables);
					return ret;
				}

				for (j = 0; j < 8; j++) {
					if (options == ARTIO_READ_ALL || 
							(options == ARTIO_READ_REFINED && oct_refined[j]) ||
							(options == ARTIO_READ_LEAFS && !oct_refined[j]) ) {
						callback(&variables[j * ghandle->num_grid_variables],
								level, oct_refined[j], sfc);
					}
				}
			}
			artio_grid_read_level_end(handle);
		}
		artio_grid_read_root_cell_end(handle);
	}

	free(variables);
	free(octs_per_level);

	return ARTIO_SUCCESS;
}
int artio_grid_read_sfc_range_buffer(artio_file handle, 
		int64_t sfc1, int64_t sfc2,
		int min_level_to_read, int max_level_to_read, 
		int options,
                GridCallBackBuffer callback, 
                void *userdata) {
	int64_t sfc;
	int oct, level, j;
	int ret;
	int *octs_per_level = NULL;
	int refined;
	int oct_refined[8];
	int root_tree_levels;
	float * variables = NULL;

	artio_grid_file ghandle;

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if (handle->open_mode != ARTIO_FILESET_READ ||
			!(handle->open_type & ARTIO_OPEN_GRID) ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	ghandle = handle->grid;

	if ((min_level_to_read < 0) || (min_level_to_read > max_level_to_read)) {
		return ARTIO_ERR_INVALID_LEVEL;
	}

	octs_per_level = (int *)malloc(ghandle->file_max_level * sizeof(int));
	variables = (float *)malloc(8*ghandle->num_grid_variables * sizeof(float));

	if ( octs_per_level == NULL || variables == NULL ) {
		if ( octs_per_level != NULL ) free(octs_per_level);
		if ( variables != NULL ) free(variables);
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}

	ret = artio_grid_cache_sfc_range(handle, sfc1, sfc2);
	if ( ret != ARTIO_SUCCESS ) {
		free(octs_per_level);
		free(variables);
		return ret;
	}
		
	for (sfc = sfc1; sfc <= sfc2; sfc++) {
		ret = artio_grid_read_root_cell_begin(handle, sfc, variables,
				&root_tree_levels, octs_per_level);
		if ( ret != ARTIO_SUCCESS ) {
			free(octs_per_level);
			free(variables);
			return ret;
		}

		if (min_level_to_read == 0 && (options == ARTIO_READ_ALL || 
				(options == ARTIO_READ_REFINED && root_tree_levels > 0) || 
				(options == ARTIO_READ_LEAFS && root_tree_levels == 0)) ) {
			refined = (root_tree_levels > 0) ? 1 : 0;
                        j=-1;
			callback(variables, 0, refined, j, userdata);
		}

		for (level = MAX(min_level_to_read,1); 
				level <= MIN(root_tree_levels,max_level_to_read); level++) {
			ret = artio_grid_read_level_begin(handle, level);
			if ( ret != ARTIO_SUCCESS ) {
				free(octs_per_level);
				free(variables);
				return ret;
			}

			for (oct = 0; oct < octs_per_level[level - 1]; oct++) {
				ret = artio_grid_read_oct(handle, variables, oct_refined);
				if ( ret != ARTIO_SUCCESS ) {
					free(octs_per_level);
					free(variables);
					return ret;
				}

				for (j = 0; j < 8; j++) {
					if (options == ARTIO_READ_ALL || 
							(options == ARTIO_READ_REFINED && oct_refined[j]) ||
							(options == ARTIO_READ_LEAFS && !oct_refined[j]) ) {
						callback(&variables[j * ghandle->num_grid_variables],
                                                         level, oct_refined[j], j, userdata);
					}
				}
			}
			artio_grid_read_level_end(handle);
		}
		artio_grid_read_root_cell_end(handle);
	}

	free(variables);
	free(octs_per_level);

	return ARTIO_SUCCESS;
}


/* array which describes how child cells are offset from 
 * the corner of their parent oct */
#define num_children 8
#define min_level 0
const double cell_delta_corner[num_children][nDim] = {
#if nDim == 1
    { 0.0 }, { 1.0 }
#elif nDim == 2
    { 0.0 , 0.0  }, { 1.0, 0.0  }, { 0.0 , 1.0 }, { 1.0, 1.0 }
#elif nDim == 3
    { 0.0 , 0.0 , 0.0  }, {  1.0, 0.0 , 0.0  }, { 0.0 ,  1.0, 0.0  },
    {  1.0,  1.0, 0.0  }, { 0.0 , 0.0 ,  1.0 }, {  1.0, 0.0 ,  1.0 },
    { 0.0 ,  1.0,  1.0 }, {  1.0,  1.0,  1.0 }
#else
#error "No valid cell_delta for that number of dimensions!"
#endif
};
int artio_grid_read_sfc_range_pos(artio_file handle, 
                                  int64_t sfc1, int64_t sfc2,
                                  int min_level_to_read, int max_level_to_read, 
                                  int options,
                                  GridCallBackPos callback, 
                                  void *user_data
    ) {
	int64_t sfc;
	int oct, level, j, i;
	int ret;
	int *num_octs_per_level = NULL;
        int num_level_octs, num_next_level_octs;
        int num_root_levels;

        int coords[nDim];
        double pos[nDim];
        double cell_size;
        double *level_octs_pos_x, *next_level_octs_pos_x;
        double *level_octs_pos_y, *next_level_octs_pos_y;
        double *level_octs_pos_z, *next_level_octs_pos_z;
	int refined;
	int oct_refined[num_children];
	int root_tree_levels;
	float * variables = NULL;

	artio_grid_file ghandle;

	if ( handle == NULL ) {
		return ARTIO_ERR_INVALID_HANDLE;
	}

	if (handle->open_mode != ARTIO_FILESET_READ ||
			!(handle->open_type & ARTIO_OPEN_GRID) ) {
		return ARTIO_ERR_INVALID_FILESET_MODE;
	}

	ghandle = handle->grid;

	if ((min_level_to_read < 0) || (min_level_to_read > max_level_to_read)) {
		return ARTIO_ERR_INVALID_LEVEL;
	}

        num_root_levels = log(handle->num_root_cells)/(3*log(2));
	num_octs_per_level = (int *)malloc(ghandle->file_max_level * sizeof(int));
	variables = (float *)malloc(num_children*ghandle->num_grid_variables * sizeof(float));

	if ( num_octs_per_level == NULL || variables == NULL ) {
		if ( num_octs_per_level != NULL ) free(num_octs_per_level);
		if ( variables != NULL ) free(variables);
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}

	ret = artio_grid_cache_sfc_range(handle, sfc1, sfc2);
	if ( ret != ARTIO_SUCCESS ) {
		free(num_octs_per_level);
		free(variables);
		return ret;
	}
		
	for (sfc = sfc1; sfc <= sfc2; sfc++) {
		ret = artio_grid_read_root_cell_begin(handle, sfc, variables,
				&root_tree_levels, num_octs_per_level);

		if ( ret != ARTIO_SUCCESS ) {
			free(num_octs_per_level);
			free(variables);
			return ret;
		}
                //////////////////////////////////////
                sfc_coords(sfc, coords, num_root_levels); 
                for(i=0; i<nDim; i++){ pos[i] = (double)coords[i]; }
                refined = (root_tree_levels > 0) ? 1 : 0;
                if( refined ){
                    
                    num_next_level_octs = 1;
                    next_level_octs_pos_x = malloc( sizeof(double)*num_next_level_octs );
                    next_level_octs_pos_y = malloc( sizeof(double)*num_next_level_octs );
                    next_level_octs_pos_z = malloc( sizeof(double)*num_next_level_octs );
                    if ( next_level_octs_pos_x == NULL ||
                         next_level_octs_pos_y == NULL ||
                         next_level_octs_pos_z == NULL ){
                        return ARTIO_ERR_MEMORY_ALLOCATION;
                    }
                    next_level_octs_pos_x[0] = pos[0] ;  //need a stupid 2D array
                    next_level_octs_pos_y[0] = pos[1] ;  //need a stupid 2D array
                    next_level_octs_pos_z[0] = pos[2] ;  //need a stupid 2D array
                }else{
                    num_next_level_octs = 0;
                }
                //////////////////////////////////////
                        

		if (min_level_to_read == 0 && 
                    (options == ARTIO_READ_ALL || 
                     (options == ARTIO_READ_REFINED && root_tree_levels > 0) || 
                     (options == ARTIO_READ_LEAFS && root_tree_levels == 0)) ) {
			refined = (root_tree_levels > 0) ? 1 : 0;
			callback(variables, min_level, refined, sfc, pos, user_data);
		}
                // level is the cell_level; current octs live at level-1.
		for (level = min_level+1; level <= MIN(root_tree_levels,max_level_to_read); level++) { 
			ret = artio_grid_read_level_begin(handle, level);
			if ( ret != ARTIO_SUCCESS ) {
				free(num_octs_per_level);
				free(variables);
				return ret;
			}
                        
                        //////////////////////////////////////
                        if( num_next_level_octs > 0 ){
                            num_level_octs = num_next_level_octs;
                            level_octs_pos_x = next_level_octs_pos_x;
                            level_octs_pos_y = next_level_octs_pos_y;
                            level_octs_pos_z = next_level_octs_pos_z;
                            
                            // at max_level-1 the next level can only have cells -- no more octs
                            if(num_level_octs != num_octs_per_level[level-1]){
                                printf("bad oct counting next expected:%d octs per level array %d\n",
                                       num_level_octs, num_octs_per_level[level-1]);
                                       printf("sfc %d \n",sfc);
                            }
                            if ( level < MIN(root_tree_levels,max_level_to_read) && num_octs_per_level[level] > 0 ) {
                                next_level_octs_pos_x = malloc( sizeof(double)*num_children*num_level_octs );
                                next_level_octs_pos_y = malloc( sizeof(double)*num_children*num_level_octs );
                                next_level_octs_pos_z = malloc( sizeof(double)*num_children*num_level_octs );
                                if ( next_level_octs_pos_x == NULL ||
                                     next_level_octs_pos_y == NULL ||
                                     next_level_octs_pos_z == NULL ){
                                    return ARTIO_ERR_MEMORY_ALLOCATION;
                                }
                            }
                        }
                        num_next_level_octs = 0;
                        cell_size = pow(2,-level);
                        //////////////////////////////////////
                        
			for (oct = 0; oct < num_octs_per_level[level - 1]; oct++) {
				ret = artio_grid_read_oct(handle, variables, oct_refined);
				if ( ret != ARTIO_SUCCESS ) {
					free(num_octs_per_level);
					free(variables);
					return ret;
				}

				for (j = 0; j < num_children; j++) {
                                        pos[0] = level_octs_pos_x[oct]+cell_delta_corner[j][0]*cell_size;
                                        pos[1] = level_octs_pos_y[oct]+cell_delta_corner[j][1]*cell_size;
                                        pos[2] = level_octs_pos_z[oct]+cell_delta_corner[j][2]*cell_size;
                                        if(oct_refined[j]){
                                                next_level_octs_pos_x[num_next_level_octs] = pos[0];
                                                next_level_octs_pos_y[num_next_level_octs] = pos[1];
                                                next_level_octs_pos_z[num_next_level_octs] = pos[2];
                                                num_next_level_octs++ ;
                                        }
					if (options == ARTIO_READ_ALL || 
                                            (options == ARTIO_READ_REFINED && oct_refined[j]) ||
                                            (options == ARTIO_READ_LEAFS && !oct_refined[j]) ) {
                                            callback(variables, level, oct_refined[j], sfc, pos, user_data);
								
					}
				}
			}
                        free(level_octs_pos_x);
                        free(level_octs_pos_y);
                        free(level_octs_pos_z);
			artio_grid_read_level_end(handle);
		}
		artio_grid_read_root_cell_end(handle);
	}

	free(variables);
	free(num_octs_per_level);

	return ARTIO_SUCCESS;
}
