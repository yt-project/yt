#ifndef __SFC_H__
#define __SFC_H__

#define nDim            3

#define	SLAB		0
#define MORTON		1
#define HILBERT		2

#ifndef SFC
#define SFC		HILBERT
#else
#if SFC == SLAB
#ifndef SLAB_DIM
#define SLAB_DIM	1
#endif
#endif
#endif


#define rollLeft(x,y,mask) ((x<<y) | (x>>(nDim-y))) & mask
#define rollRight(x,y,mask) ((x>>y) | (x<<(nDim-y))) & mask

//int sfc_index( int coords[nDim], int root_level );
int sfc_index0( int ix, int iy, int iz, int root_level );
int sfc_index( int coords[nDim], int root_level );
void sfc_coords( int index, int coords[nDim], int root_level );

#endif
