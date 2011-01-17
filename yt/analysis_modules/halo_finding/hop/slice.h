/* SLICE.H, Daniel Eisenstein, 1997 */
/* Based on a paper by Daniel Eisenstein & Piet Hut, 
"HOP: A New Group-Finding Algorithm for N-body Simulations."
See the included documentation or view it at 
http://www.sns.ias.edu/~eisenste/hop/hop_doc.html */

/* Version 1.0 (12/15/97) -- Original Release */

#ifndef NBODYUTIL_H
#define NBODYUTIL_H

#define RHOCRIT 277.5	/* in h^2 Msun/kpc^3 */
#define HTIME 9.7776	/* in h^-1 Gyr */
#define GNEWT 4.51e-6   /* in kpc^3/Msun/Gyr^2, h cancels */
#define KMS 0.975       /* to convert from kpc/Gyr to km/s */

#define PI 3.141592654
#define ROOT2 1.414213562
#define ROOTPI2 1.253314137     /* sqrt(Pi/2) */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "macros_and_parameters.h"

/* Structure to hold data from one time slice and info about the slice */
/* This needn't hold all of the points; it might just hold some subset */
/* Hence, data vectors are just given as arrays here */
/* Usage: pid[] Array goes from 0 to pid[0].  pid is NULL if the particle
numbering is just the trivial mapping.  Kinematic arrays go from 1 to numlist
and can be NULL if the data hasn't been read.  Numlist==pid[0] often */

typedef struct slicestruct {
#ifdef NOT_USED
	/* First, some generic stuff about the simulation */
	float omega, lambda, curv;
	float h0;	/* H_0 = 100 km/s/Mpc * h0 */
	float specn, gamma;	/* The spectral index and BBKS PS Gamma */
	float sigma8;	/* At z=0 */

	/* Next, some information about this slice in particular */
	float z, a, t, growth; 
		/* a=1 at z=0, g=a at early times */
	float hubb;	/* This is a*H(z), but without h0.  So it's 0.1 km/s/kpc
		   redshifted appropriately.  This is used to relate 
		   comoving positions and peculiar velocities */

	/* Now some information about the data */
	int numblocks;	/* Number of blocks in the data file */
	int numperblock;  /* Number of particles per block */
	float masspart;	/* Mass per particle in h^-1 Msun */
	float boxsize;	/* Comoving size of box in h^-1 kpc */
	float physsize;	/* Physical Size in h^-1 kpc */
	float velscale;	/* To turn raw velocity data into peculiar vel in km/s*/
#endif

	int numpart;    /* Total number of particles in the simulation */
	/* Now the data itself */
	int *pid;	/* The id number of the particle */
			/* pid[0] holds the number of particles in the list */
	int offset;	/* If pid==NULL, then the arrays are consecutively
				numbered, starting from offset+1 */
	int numlist;	/* Length of arrays below, set when allocated */
	float *px, *py, *pz, *vx, *vy, *vz;	/* The kinematic information */

	/* And here's the group tag information */
	int *ntag;	/* Only stored for the numlist above */
        //int *ID;       /* The real, true ID of the particle. S Skory */
	int numgroups;	/* The number of groups read out of the tag file */
} Slice;	/* Type Slice is defined */

/* Prototypes */
Slice *newslice();
void free_tags(Slice *s);
void free_data(Slice *s);
void free_slice(Slice *s);
int f77write(FILE *f, void *p, int len);
int f77read(FILE *f, void *p, int len);
float *vector(long nl, long nh);
int *ivector(long nl, long nh);
void free_vector(float *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);

void myerror(char *message);
void mywarn(char *message);

int read_header(FILE *f, Slice *s);
void normalizedata(Slice *s, int conp, int conv);
int read_alldata(FILE *f, FILE *ftag, Slice *s, int conp, int conv);
int read_partdata(FILE *f, FILE *ftag, Slice *s);

int readtag(FILE *f, int numread, int *ntag);
int skiptagheader(FILE *f, Slice *s);
int readalltags(FILE *f, Slice *s);


#endif 		/* NBODYUTIL_H */
