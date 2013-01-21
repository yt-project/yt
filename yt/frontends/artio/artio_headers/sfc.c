#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "sfc.h"


int nBitsPerDim ;
int nBits;
int max_sfc_index ;
int num_grid ;

/*******************************************************
 * morton_index
 ******************************************************/
int morton_index( int coords[nDim] ) 
/* purpose: interleaves the bits of the nDim integer
 * 	coordinates, normally called Morton or z-ordering
 *
 * 	Used by the hilbert curve algorithm
 */
{
	int i, d;
	int mortonnumber = 0;
	int bitMask = 1 << (nBitsPerDim - 1);
	
	/* interleave bits of coordinates */
	for ( i = nBitsPerDim; i > 0; i-- ) {
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
int hilbert_index( int coords[nDim] ) 
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
	int hilbertnumber;
	int singleMask;
	int dimMask;
	int numberShifts;
	int principal;
	int o;
	int rho;
	int w;
	int interleaved;

	/* begin by transposing bits */
	interleaved = morton_index( coords );

	/* mask out nDim and 1 bit blocks starting 
	 * at highest order bits */
	singleMask = 1 << ((nBitsPerDim - 1) * nDim);
	
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
void hilbert_coords( int index, int coords[nDim] ) 
/* purpose: performs the inverse of sfc_index,
 * 	taking a 1-d space-filling-curve index
 * 	and transforming it into nDim coordinates
 *
 * returns: the coordinates in coords
 */
{
	int i, j;
	int dimMask;
	int singleMask;
	int sigma;
	int sigma_;
	int tau;
	int tau_;
	int num_shifts;
	int principal;
	int w;
	int x = 0;

	if( index < 0 ){ fprintf(stderr, "bad index\n"); exit(-1);}

	w = 0;
	sigma_ = 0;
	num_shifts = 0;

        singleMask = 1 << ((nBitsPerDim - 1) * nDim);

        dimMask = singleMask;
        for ( i = 1; i < nDim; i++ ) {
                dimMask |= singleMask << i;
        }

	for ( i = 0; i < nBitsPerDim; i++ ) {
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

		singleMask = 1 << (nBits - 1 - i);

		for ( j = 0; j < nBitsPerDim; j++ ) {
			if ( x & singleMask ) {
				coords[i] |= 1 << (nBitsPerDim-j-1);
			}
			singleMask >>= nDim;
		}
	}
}

int slab_index( int coords[nDim] ) {
	#if SLAB_DIM == 0
		return num_grid*num_grid*coords[0] + num_grid*coords[1] + coords[2];
	#elif SLAB_DIM == 1
		return num_grid*num_grid*coords[1] + num_grid*coords[0] + coords[2];
	#else
		return num_grid*num_grid*coords[2] + num_grid*coords[0] + coords[1];
	#endif
}

void slab_coords( int index, int coords[nDim] ) {
	#if SLAB_DIM == 0
		coords[2] = index % num_grid;
		coords[1] = ((index - coords[2] )/num_grid) % num_grid;
		coords[0] = (index - coords[2] - num_grid*coords[1])/(num_grid*num_grid);
	#elif SLAB_DIM == 1
		coords[2] = index % num_grid;
		coords[0] = ((index - coords[2] )/num_grid) % num_grid;
		coords[1] = (index - coords[2] - num_grid*coords[0])/(num_grid*num_grid);
	#else
		coords[1] = index % num_grid;
		coords[0] = ((index - coords[1] )/num_grid) % num_grid;
		coords[2] = (index - coords[1] - num_grid*coords[0])/(num_grid*num_grid);
	#endif

	if( slab_index( coords ) != index ){ fprintf(stderr, "bad slab index\n"); exit(-1);}
}

int sfc_index0( int ix,int iy, int iz, int num_root_grid_refinements ) {
    int coords[nDim] = {ix,iy,iz};
        num_grid = pow(2, num_root_grid_refinements );
        nBitsPerDim =    num_root_grid_refinements;
        nBits    =       (nDim * nBitsPerDim);
        max_sfc_index  =  (1<<nBits);
        
	#if SFC == SLAB
		return slab_index( coords );
	#elif SFC == MORTON
		#error "Morton order not supported yet!"
	#else
		return hilbert_index( coords );
	#endif
}

//int sfc_index( int coords[nDim], int num_root_grid_refinements ) {
int sfc_index( int coords[nDim], int num_root_grid_refinements ) {
        num_grid = pow(2, num_root_grid_refinements );
        nBitsPerDim =    num_root_grid_refinements;
        nBits    =       (nDim * nBitsPerDim);
        max_sfc_index  =  (1<<nBits);
        
	#if SFC == SLAB
		return slab_index( coords );
	#elif SFC == MORTON
		#error "Morton order not supported yet!"
	#else
		return hilbert_index( coords );
	#endif
}

void sfc_coords( int index, int coords[nDim], int num_root_grid_refinements ) {
        num_grid = pow(2, num_root_grid_refinements );
        nBitsPerDim =    num_root_grid_refinements;
        nBits    =       (nDim * nBitsPerDim);
        max_sfc_index  =  (1<<nBits);
	#if SFC == SLAB
		slab_coords( index, coords );
	#elif SFC == MORTON
		#error "Morton order not supported yet!"
	#else
		hilbert_coords( index, coords );
	#endif
}
