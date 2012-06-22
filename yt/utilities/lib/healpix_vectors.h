/* -----------------------------------------------------------------------------
 *
 *  Copyright (C) 1997-2010 Krzysztof M. Gorski, Eric Hivon,
 *                          Benjamin D. Wandelt, Anthony J. Banday, 
 *                          Matthias Bartelmann, 
 *                          Reza Ansari & Kenneth M. Ganga 
 *
 *
 *  This file is part of HEALPix.
 *
 *  HEALPix is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  HEALPix is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with HEALPix; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix see http://healpix.jpl.nasa.gov
 *
 *----------------------------------------------------------------------------- */
/* chealpix.h
 *
 */

#ifndef __CHEALPIX_H__
#define __CHEALPIX_H__

#ifdef __cplusplus
extern "C" {
#endif

/* -------------------- */
/* Constant Definitions */
/* -------------------- */

#ifndef HEALPIX_NULLVAL
#define HEALPIX_NULLVAL (-1.6375e30)
#endif /* HEALPIX_NULLVAL */

/* --------------------- */
/* Function Declarations */
/* --------------------- */

/* pixel operations */
/* ---------------- */
void ang2pix_nest(const long nside, double theta, double phi, long *ipix);
void ang2pix_ring(const long nside, double theta, double phi, long *ipix);

void pix2ang_nest(long nside, long ipix, double *theta, double *phi);
void pix2ang_ring(long nside, long ipix, double *theta, double *phi);

void nest2ring(long nside, long ipnest, long *ipring);
void ring2nest(long nside, long ipring, long *ipnest);

void mk_pix2xy(int *pix2x, int *pix2y);
void mk_xy2pix(int *x2pix, int *y2pix);

long nside2npix(const long nside);
long npix2nside(const long pix  );

void ang2vec(double theta, double phi,   double *vec);
void vec2ang(double *vec, double *theta, double *phi);

void vec2pix_nest(const long nside, double *vec, long *ipix);
void vec2pix_ring(const long nside, double *vec, long *ipix);

void pix2vec_nest(long nside, long ipix, double *vec);
void pix2vec_ring(long nside, long ipix, double *vec);


/* FITS operations */
/* --------------- */

void printerror (int) ;

float *read_healpix_map (const char *, long *, char *, char *) ;

int write_healpix_map( float *, long , const char *, char ,char *) ;

long get_fits_size(char *, long *, char * ) ;


/* ------------------ */
/* end of header file */
/* ------------------ */

#ifdef __cplusplus
} /* closing brace for extern "C" */
#endif

#endif /* __CHEALPIX_H__ */

