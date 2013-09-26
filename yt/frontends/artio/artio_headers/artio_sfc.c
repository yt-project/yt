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

#define rollLeft(x,y,mask) ((x<<y) | (x>>(nDim-y))) & mask
#define rollRight(x,y,mask) ((x>>y) | (x<<(nDim-y))) & mask

/*******************************************************
 * morton_index
 ******************************************************/
int64_t artio_morton_index( artio_fileset *handle, int coords[nDim] ) 
/* purpose: interleaves the bits of the nDim integer
 * 	coordinates, normally called Morton or z-ordering
 *
 * 	Used by the hilbert curve algorithm
 */
{
	int i, d;
	int64_t mortonnumber = 0;
	int64_t bitMask = 1L << (handle->nBitsPerDim - 1);
	
	/* interleave bits of coordinates */
	for ( i = handle->nBitsPerDim; i > 0; i-- ) {
		for ( d = 0; d < nDim; d++ ) {
			mortonnumber |= ( coords[d] & bitMask ) << (((nDim - 1) * i ) - d );
		}
		bitMask >>= 1;
	}

	return mortonnumber;
}

/*******************************************************
 * hilbert_index
 ******************************************************/
int64_t artio_hilbert_index( artio_fileset *handle, int coords[nDim] ) 
/* purpose: calculates the 1-d space-filling-curve index
 * 	corresponding to the nDim set of coordinates
 *
 * 	Uses the Hilbert curve algorithm given in
 * 	Alternative Algorithm for Hilbert's Space-
 * 	Filling Curve, A.R. Butz, IEEE Trans on Comp.,
 * 	p. 424, 1971
 */
{
	int i, j;
	int64_t hilbertnumber;
	int64_t singleMask;
	int64_t dimMask;
	int64_t numberShifts;
	int principal;
	int64_t o;
	int64_t rho;
	int64_t w;
	int64_t interleaved;

	/* begin by transposing bits */
	interleaved = artio_morton_index( handle, coords );

	/* mask out nDim and 1 bit blocks starting 
	 * at highest order bits */
	singleMask = 1L << ((handle->nBitsPerDim - 1) * nDim);
	
	dimMask = singleMask;
	for ( i = 1; i < nDim; i++ ) {
		dimMask |= singleMask << i;
	}

	w = 0;
	numberShifts = 0;
	hilbertnumber = 0;

	while (singleMask) {
		o = (interleaved ^ w) & dimMask;
		o = rollLeft( o, numberShifts, dimMask );

		rho = o;
		for ( j = 1; j < nDim; j++ ) {
			rho ^= (o>>j) & dimMask;
		}

		hilbertnumber |= rho;

		/* break out early (we already have complete number
		 * no need to calculate other numbers) */
		if ( singleMask == 1 ) {
			break;
		}

		/* calculate principal position */
		principal = 0;
		for ( i = 1; i < nDim; i++ ) {
			if ( (hilbertnumber & singleMask) != ((hilbertnumber>>i) & singleMask)) {
				principal = i; 
				break;
			}
		}

		/* complement lowest bit position */
		o ^= singleMask;

		/* force even parity by complementing at principal position if necessary 
		 * Note: lowest order bit of hilbertnumber gives you parity of o at this
		 * point due to xor operations of previous steps */
		if ( !(hilbertnumber & singleMask) ) {
			o ^= singleMask << principal;
		}

		/* rotate o right by numberShifts */
		o = rollRight( o, numberShifts, dimMask );

		/* find next numberShifts */
		numberShifts += (nDim - 1) - principal;
		numberShifts %= nDim;

		w ^= o;
		w >>= nDim;

		singleMask >>= nDim;
		dimMask >>= nDim;
	}

	return hilbertnumber;
}

/*******************************************************
 * hilbert_coords
 ******************************************************/
void artio_hilbert_coords( artio_fileset *handle, int64_t index, int coords[nDim] ) 
/* purpose: performs the inverse of sfc_index,
 * 	taking a 1-d space-filling-curve index
 * 	and transforming it into nDim coordinates
 *
 * returns: the coordinates in coords
 */
{
	int i, j;
	int64_t dimMask;
	int64_t singleMask;
	int64_t sigma;
	int64_t sigma_;
	int64_t tau;
	int64_t tau_;
	int num_shifts;
	int principal;
	int64_t w;
	int64_t x = 0;

	w = 0;
	sigma_ = 0;
	num_shifts = 0;

	singleMask = 1L << ((handle->nBitsPerDim - 1) * nDim);

	dimMask = singleMask;
	for ( i = 1; i < nDim; i++ ) {
		dimMask |= singleMask << i;
	}

	for ( i = 0; i < handle->nBitsPerDim; i++ ) {
		sigma = ((index & dimMask) ^ ( (index & dimMask) >> 1 )) & dimMask;
		sigma_ |= rollRight( sigma, num_shifts, dimMask );

		principal = nDim - 1;
		for ( j = 1; j < nDim; j++ ) {
			if ( (index & singleMask) != ((index >> j) & singleMask) ) {
				principal = nDim - j - 1;
				break;
			}
		}

		/* complement nth bit */
		tau = sigma ^ singleMask;

		/* if even parity, complement in principal bit position */
		if ( !(index & singleMask) ) {
			tau ^= singleMask << ( nDim - principal - 1 ); 
		}

		tau_ = rollRight( tau, num_shifts, dimMask );

		num_shifts += principal;
		num_shifts %= nDim;

		w |= ((w & dimMask) ^ tau_) >> nDim;

		dimMask >>= nDim;
		singleMask >>= nDim;
	}

	x = w ^ sigma_;

	/* undo bit interleaving to get coordinates */
	for ( i = 0; i < nDim; i++ ) {
		coords[i] = 0;

		singleMask = 1L << (nDim*handle->nBitsPerDim - 1 - i);

		for ( j = 0; j < handle->nBitsPerDim; j++ ) {
			if ( x & singleMask ) {
				coords[i] |= 1 << (handle->nBitsPerDim-j-1);
			}
			singleMask >>= nDim;
		}
	}
}

int64_t artio_slab_index( artio_fileset *handle, int coords[nDim], int slab_dim ) {
	int64_t num_grid = 1L << handle->nBitsPerDim;
	int64_t index;

	switch ( slab_dim ) {
		case 0: 
			index = num_grid*num_grid*coords[0] + num_grid*coords[1] + coords[2];
			break;
		case 1: 
			index = num_grid*num_grid*coords[1] + num_grid*coords[0] + coords[2];
			break;
		case 2: 
			index = num_grid*num_grid*coords[2] + num_grid*coords[0] + coords[1];
			break;
		default:
			index = -1;
	}
	return index;
}

void artio_slab_coords( artio_fileset *handle, int64_t index, int coords[nDim], int slab_dim ) {
	int64_t num_grid = 1L << handle->nBitsPerDim;
	switch ( slab_dim ) {
		case 0:
			coords[2] = index % num_grid;
			coords[1] = ((index - coords[2] )/num_grid) % num_grid;
			coords[0] = (index - coords[2] - num_grid*coords[1])/(num_grid*num_grid);
			break;
		case 1:
			coords[2] = index % num_grid;
			coords[0] = ((index - coords[2] )/num_grid) % num_grid;
			coords[1] = (index - coords[2] - num_grid*coords[0])/(num_grid*num_grid);
			break;
		case 2:	
			coords[1] = index % num_grid;
			coords[0] = ((index - coords[1] )/num_grid) % num_grid;
			coords[2] = (index - coords[1] - num_grid*coords[0])/(num_grid*num_grid);
			break;
	}
}

int64_t artio_sfc_index_position( artio_fileset *handle, double position[nDim] ) {
	int i;
	int coords[nDim];

	for ( i = 0; i < nDim; i++ ) {
		coords[i] = (int)position[i];
	}
	return artio_sfc_index(handle, coords);
}

int64_t artio_sfc_index( artio_fileset *handle, int coords[nDim] ) {
	switch ( handle->sfc_type ) {
		case ARTIO_SFC_SLAB_X: return artio_slab_index(handle, coords, 0);
		case ARTIO_SFC_SLAB_Y: return artio_slab_index(handle, coords, 1);
		case ARTIO_SFC_SLAB_Z: return artio_slab_index(handle, coords, 2);
		case ARTIO_SFC_HILBERT: return artio_hilbert_index( handle, coords );
		default: return -1;
	}
}

void artio_sfc_coords( artio_fileset *handle, int64_t index, int coords[nDim] ) {
	int i;

	switch ( handle->sfc_type ) {
		case ARTIO_SFC_SLAB_X: 
			artio_slab_coords( handle, index, coords, 0 );
			break;
		case ARTIO_SFC_SLAB_Y: 
			artio_slab_coords( handle, index, coords, 1 );
			break;
		case ARTIO_SFC_SLAB_Z: 
			artio_slab_coords( handle, index, coords, 2 );
			break;
		case ARTIO_SFC_HILBERT: 
			artio_hilbert_coords( handle, index, coords );	
			break;
		default :
			for ( i = 0; i < nDim; i++ ) {
				coords[i] = -1;
			}
			break;
	}
}
