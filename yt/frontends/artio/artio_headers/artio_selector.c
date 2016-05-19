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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef _WIN32
typedef __int64 int64_t;
typedef __int32 int32_t;
#else
#include <stdint.h>
#endif

#define ARTIO_SELECTION_LIST_SIZE		1024
#define ARTIO_SELECTION_VOLUME_LIMIT	(1L<<60)

int artio_add_volume_to_selection( artio_fileset *handle, int lcoords[3], int rcoords[3],
            int64_t sfcs[8], artio_selection *selection );

int artio_selection_iterator( artio_selection *selection, 
		 int64_t max_range_size, int64_t *start, int64_t *end ) {

	if ( selection->cursor < 0 ) {
		selection->cursor = 0;
	}

	if ( selection->cursor == selection->num_ranges ) {
		selection->cursor = -1;
		return ARTIO_SELECTION_EXHAUSTED;
	}

	if ( selection->subcycle > 0 ) {
		*start = selection->subcycle+1;
	} else {
		*start = selection->list[2*selection->cursor];
	}

	*end = selection->list[2*selection->cursor+1];
	if ( *end - *start > max_range_size ) {
		*end = *start + max_range_size-1;
		selection->subcycle = *end;
	} else {
		selection->subcycle = -1;
		selection->cursor++;
	}	

	return ARTIO_SUCCESS;
}

int artio_selection_iterator_reset( artio_selection *selection ) {
	selection->cursor = -1;
	selection->subcycle = -1;
	return ARTIO_SUCCESS;
}

int64_t artio_selection_size( artio_selection *selection ) {
	int i;
	int64_t count = 0;
	for ( i = 0; i < selection->num_ranges; i++ ) {
		count += selection->list[2*i+1] - selection->list[2*i] + 1;
	}
	return count;
}

artio_selection *artio_selection_allocate( artio_fileset *handle ) {
	artio_selection *selection = (artio_selection *)malloc(sizeof(artio_selection));
	if ( selection != NULL ) {
		selection->list = (int64_t *)malloc(2*ARTIO_SELECTION_LIST_SIZE*sizeof(int64_t));
		if ( selection->list == NULL ) {
			free(selection);
			return NULL;
		}
	}
	selection->subcycle = -1;
	selection->cursor = -1;
	selection->size = ARTIO_SELECTION_LIST_SIZE;
	selection->num_ranges = 0;
	selection->fileset = handle;
	return selection;
}

int artio_selection_destroy( artio_selection *selection ) {
	if ( selection == NULL ) {
		return ARTIO_ERR_INVALID_SELECTION;
	}

	if ( selection->list != NULL ) {
		free( selection->list );
	}
	free(selection);
	return ARTIO_SUCCESS;
}

int artio_selection_add_range( artio_selection *selection, 
		int64_t start, int64_t end ) {
	int i, j;
	int64_t *new_list;

	if ( selection == NULL ) {
		return ARTIO_ERR_INVALID_SELECTION;
	}

	if ( start < 0 || end >= selection->fileset->num_root_cells || start > end ) {
		return ARTIO_ERR_INVALID_SFC_RANGE;
	}

	for ( i = 0; i < selection->num_ranges; i++ ) {
		if ( (start >= selection->list[2*i] && start <= selection->list[2*i+1] ) ||
			 (end >= selection->list[2*i] && end <= selection->list[2*i+1] ) ) {
			return ARTIO_ERR_INVALID_STATE;
		}
	}

	/* locate page */
	if ( selection->num_ranges == 0 ) {
		selection->list[0] = start;
		selection->list[1] = end;
		selection->num_ranges = 1;
		return ARTIO_SUCCESS;
	} else {
		/* eventually replace with binary search */
		for ( i = 0; i < selection->num_ranges; i++ ) {
			if ( end < selection->list[2*i] ) {
				break;
			}
		}

		if ( ( i == 0 && end < selection->list[2*i]-1 ) ||
				( i == selection->num_ranges && start > selection->list[2*i-1]+1 ) ||
				( end < selection->list[2*i]-1 && start > selection->list[2*i-1]+1 ) ) { 
			if ( selection->num_ranges == selection->size ) {
				new_list = (int64_t *)malloc(4*selection->size*sizeof(int64_t));
				if ( new_list == NULL ) {
					return ARTIO_ERR_MEMORY_ALLOCATION;
				}

				for ( j = 0; j < i; j++ ) {
					new_list[2*j] = selection->list[2*j];
					new_list[2*j+1] = selection->list[2*j+1];
				}
				for ( ; j < selection->num_ranges; j++ ) {
					new_list[2*j+2] = selection->list[2*j];
					new_list[2*j+3] = selection->list[2*j+1];
				}
				selection->size *= 2;
				free( selection->list );
				selection->list = new_list;
			} else {
				for ( j = selection->num_ranges-1; j >= i; j-- ) {
					selection->list[2*j+2] = selection->list[2*j];
					selection->list[2*j+3] = selection->list[2*j+1];
				}
			}

			selection->list[2*i] = start;
			selection->list[2*i+1] = end;
			selection->num_ranges++;
		} else {
			if ( end == selection->list[2*i]-1 ) {
				selection->list[2*i] = start;
			} else if ( start == selection->list[2*i-1]+1 ) {
				selection->list[2*i-1] = end;
			}

			/* merge 2 ranges if necessary */
			if ( selection->list[2*i] == selection->list[2*i-1]+1 ) {
				selection->list[2*i-1] = selection->list[2*i+1];
				selection->num_ranges--;
				for ( ; i < selection->num_ranges; i++ ) {
					selection->list[2*i] = selection->list[2*i+2];
					selection->list[2*i+1] = selection->list[2*i+3];	
				}
			}
		}
	}

	return ARTIO_SUCCESS;	
}

int artio_selection_add_root_cell( artio_selection *selection, int coords[3] ) {
	int i;
	int64_t sfc;

	if ( selection == NULL ) {
		return ARTIO_ERR_INVALID_SELECTION;
	}

	for ( i = 0; i < 3; i++ ) {
		if ( coords[i] < 0 || coords[i] >= selection->fileset->num_grid ) {
			return ARTIO_ERR_INVALID_COORDINATES;
		}
	}

	sfc = artio_sfc_index( selection->fileset, coords );
	return artio_selection_add_range( selection, sfc, sfc );
}

void artio_selection_print( artio_selection *selection ) {
	int i;
	
	for ( i = 0; i < selection->num_ranges; i++ ) {
		printf("%u: %ld %ld\n", i, selection->list[2*i], selection->list[2*i+1] );
	}
}

artio_selection *artio_select_all( artio_fileset *handle ) {
	artio_selection *selection;

	if ( handle == NULL ) {
		return NULL;
	}

	selection = artio_selection_allocate(handle);
	if ( selection == NULL ) {
		return NULL;
	}

	if ( artio_selection_add_range( selection, 0, handle->num_root_cells-1 ) != ARTIO_SUCCESS ) {
		artio_selection_destroy(selection);
		return NULL;
	}

	return selection;
}

artio_selection *artio_select_volume( artio_fileset *handle, double lpos[3], double rpos[3] ) {
	int i;
	int64_t sfc;
	int coords[3];
	int lcoords[3];
	int rcoords[3];
	artio_selection *selection;

	if ( handle == NULL ) {
		return NULL;
	}

	for ( i = 0; i < 3; i++ ) {
		if ( lpos[i] < 0.0 || lpos[i] >= rpos[i] ) {
			return NULL;
		}
	}

	for ( i = 0; i < 3; i++ ) {
		lcoords[i] = (int)lpos[i];
		rcoords[i] = (int)rpos[i];
	}

	selection = artio_selection_allocate( handle );
	if ( selection == NULL ) {
		return NULL;
	}

	for ( coords[0] = lcoords[0]; coords[0] <= rcoords[0]; coords[0]++ ) {
		for ( coords[1] = lcoords[1]; coords[1] <= rcoords[1]; coords[1]++ ) {
			for ( coords[2] = lcoords[2]; coords[2] <= rcoords[2]; coords[2]++ ) {
				sfc = artio_sfc_index( handle, coords );
				if ( artio_selection_add_range( selection, sfc, sfc ) != ARTIO_SUCCESS ) {
					artio_selection_destroy(selection);
					return NULL;
				}
			}
		}
	} 

	return selection;
}

artio_selection *artio_select_cube( artio_fileset *handle, double center[3], double size ) {
	int i, j, k, dx;
    int64_t sfc;
    int coords[3], coords2[3];
	artio_selection *selection;

	if ( handle == NULL ) {
		return NULL;
	}

	if ( size <= 0.0 || size > handle->num_grid/2 ) {
		return NULL;
	}
	dx = (int)(center[0] + 0.5*size) - (int)(center[0] - 0.5*size) + 1;

	for ( i = 0; i < 3; i++ ) {
		if ( center[i] < 0.0 || center[i] >= handle->num_grid ) {
			return NULL;
		}
		coords[i] = (int)(center[i] - 0.5*size + handle->num_grid) % handle->num_grid;
	}

	selection = artio_selection_allocate( handle );
	if ( selection == NULL ) {
		return NULL;
	}

	for ( i = coords[0]-dx; i <= coords[0]+dx; i++ ) {
		coords2[0] = (i + handle->num_grid) % handle->num_grid;
		for ( j = coords[1]-dx; j <= coords[1]+dx; j++ ) {
			coords2[1] = (j + handle->num_grid) % handle->num_grid;
			for ( k = coords[2]-dx; k <= coords[2]+dx; k++ ) {
				coords2[2] = (k + handle->num_grid) % handle->num_grid;
				sfc = artio_sfc_index( handle, coords2 );
				if ( artio_selection_add_range( selection, sfc, sfc ) != ARTIO_SUCCESS ) {
					artio_selection_destroy(selection);
					return NULL;
				}
			}
		}
	}
 			
	return selection;
} 
