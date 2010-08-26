/* KD.H */
/* This was written by Joachim Stadel and the NASA HPCC ESS at
the University of Washington Department of Astronomy as part of
the SMOOTH program, v2.0.1.
URL: http://www-hpcc.astro.washington.edu/tools/SMOOTH */

/* DJE--I have made a few alterations to the PARTICLE structure
in order to reduce memory consumption. */

/* HOP Version 1.0 (12/15/97) -- Original Release */

/* GLB--set different masses on */

#define DIFFERENT_MASSES

//#include "macros_and_parameters.h"

#ifndef KD_HINCLUDED
#define KD_HINCLUDED

#include "Python.h"
#include "numpy/ndarrayobject.h"

#define ROOT		1
#define LOWER(i)	(i<<1)
#define UPPER(i)	((i<<1)+1)
#define PARENT(i)	(i>>1)
#define SIBLING(i) 	((i&1)?i-1:i+1)
#define SETNEXT(i)\
{\
	while (i&1) i=i>>1;\
	++i;\
	}

#define DARK	1
#define GAS		2
#define STAR	4

typedef struct Particle {
    int np_index;
    int iHop;
	int iOrder;
#if 0
	float r[3];
	float fDensity;
	// int iID;  /* the real ID of the particle S. Skory */
	int iHop;	/* DJE: The number of the highest-density neighbor;
				Later, the group number. */
#ifdef DIFFERENT_MASSES
	float fMass;
#endif
	/* DJE: The following are unused and cost too much memory to keep */
	/* float v[3]; */
	/* float fMass; */
	/* int iMark; */
	/* float vMean[3]; */
	/* float fVelDisp2; */
#endif
	} PARTICLE;

typedef struct bndBound {
	float fMin[3];
	float fMax[3];
	} BND;

typedef struct kdNode {
	float fSplit;
	BND bnd;
	int iDim;
	int pLower;
	int pUpper;
	} KDN;

typedef struct kdContext {
	int nBucket;
	int nParticles;
	int nDark;
	int nGas;
	int nStar;
	int bDark;
	int bGas;
	int bStar;
	int nActive;
	float fTime;
	BND bnd;
	int nLevels;
	int nNodes;
	int nSplit;
	float fMass;	/* DJE: If all particles have the same mass */
	PARTICLE *p;
	KDN *kdNodes;
	int uSecond;
	int uMicro;
    npy_float64 *np_densities;
    npy_float64 *np_pos[3];
    npy_float64 *np_masses;
    float totalmass;
	} * KD;


#define INTERSECT(c,cp,fBall2,lx,ly,lz,x,y,z,sx,sy,sz)\
{\
	float dx,dy,dz,dx1,dy1,dz1,fDist2;\
	dx = c[cp].bnd.fMin[0]-x;\
	dx1 = x-c[cp].bnd.fMax[0];\
	dy = c[cp].bnd.fMin[1]-y;\
	dy1 = y-c[cp].bnd.fMax[1];\
	dz = c[cp].bnd.fMin[2]-z;\
	dz1 = z-c[cp].bnd.fMax[2];\
	if (dx > 0.0) {\
		dx1 += lx;\
		if (dx1 < dx) {\
			fDist2 = dx1*dx1;\
			sx = x+lx;\
			}\
		else {\
			fDist2 = dx*dx;\
			sx = x;\
			}\
		if (fDist2 > fBall2) goto GetNextCell;\
		}\
	else if (dx1 > 0.0) {\
		dx += lx;\
		if (dx < dx1) {\
			fDist2 = dx*dx;\
			sx = x-lx;\
			}\
		else {\
			fDist2 = dx1*dx1;\
			sx = x;\
			}\
		if (fDist2 > fBall2) goto GetNextCell;\
		}\
	else {\
		fDist2 = 0.0;\
		sx = x;\
		}\
	if (dy > 0.0) {\
		dy1 += ly;\
		if (dy1 < dy) {\
			fDist2 += dy1*dy1;\
			sy = y+ly;\
			}\
		else {\
			fDist2 += dy*dy;\
			sy = y;\
			}\
		if (fDist2 > fBall2) goto GetNextCell;\
		}\
	else if (dy1 > 0.0) {\
		dy += ly;\
		if (dy < dy1) {\
			fDist2 += dy*dy;\
			sy = y-ly;\
			}\
		else {\
			fDist2 += dy1*dy1;\
			sy = y;\
			}\
		if (fDist2 > fBall2) goto GetNextCell;\
		}\
	else {\
		sy = y;\
		}\
	if (dz > 0.0) {\
		dz1 += lz;\
		if (dz1 < dz) {\
			fDist2 += dz1*dz1;\
			sz = z+lz;\
			}\
		else {\
			fDist2 += dz*dz;\
			sz = z;\
			}\
		if (fDist2 > fBall2) goto GetNextCell;\
		}\
	else if (dz1 > 0.0) {\
		dz += lz;\
		if (dz < dz1) {\
			fDist2 += dz*dz;\
			sz = z-lz;\
			}\
		else {\
			fDist2 += dz1*dz1;\
			sz = z;\
			}\
		if (fDist2 > fBall2) goto GetNextCell;\
		}\
	else {\
		sz = z;\
		}\
	}


void kdTime(KD,int *,int *);
int kdInit(KD *,int);
int kdReadTipsy(KD,FILE *,int,int,int);
void kdInMark(KD,char *);
int kdBuildTree(KD);
void kdOrder(KD);
void kdFinish(KD);

#endif
