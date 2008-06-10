/* SMOOTH.H */
/* This was written by Joachim Stadel and the NASA HPCC ESS at
the University of Washington Department of Astronomy as part of
the SMOOTH program, v2.0.1.
URL: http://www-hpcc.astro.washington.edu/tools/SMOOTH */

/* DJE--I have made a few additions to the SMX structure
in order to store information necessary for HOP.  I have also
added the Boundary structure. */

/* HOP Version 1.0 (12/15/97) -- Original Release */

#ifndef SMOOTH_HINCLUDED
#define SMOOTH_HINCLUDED

#include "kd.h"
//#include "macros_and_parameters.h"

#define RESMOOTH_SAFE  30

/* DJE: Define this structure to hold the boundary data. */
typedef struct boundarystruct {
	int nGroup1, nGroup2;	/* The two groups involved, ordered such
				that nGroup1<nGroup2 */
	float fDensity;		/* The highest average density of two
				boundary particles */
	} Boundary;	/* Type Boundary is defined */


typedef struct pqNode {
	float fKey;
	struct pqNode *pqLoser;
	struct pqNode *pqFromInt;
	struct pqNode *pqFromExt;
	struct pqNode *pqWinner;	/* Only used when building initial tree */
	int p;
	float ax;
	float ay;
	float az;
	} PQ;


typedef struct smContext {
	KD kd;
	int nSmooth;
	float fPeriod[3];
	PQ *pq;
	PQ *pqHead;
	float *pfBall2;
	char *iMark;
	int nListSize;
	float *fList;
	int *pList;
	/* DJE -- Added the following fields to SMX */
	int nDens;	/* The number of neighbors for calculating density */
	int nHop;	/* How many neighbors to consider when hopping*/
	int nMerge;		/* The # of neighbors to look at when merging */
	int nGroups;		/* Groups number 0 to nGroups - 1 */
	int *nmembers;		/* The number of particles in the groups */
				/* nmembers[nGroups] is the # not in groups */
	int *densestingroup;	/* The ID # of the densest particle in the gr */
	int nHashLength;	/* The length of the hash table */
	Boundary *hash;		/* The hash table for boundaries */
	float fDensThresh;	/* Density Threshold for group finding */
	} * SMX;

#define PQ_STATIC	PQ *PQ_t,*PQ_lt;int PQ_j,PQ_i


#define PQ_INIT(pq,n)\
{\
	for (PQ_j=0;PQ_j<(n);++PQ_j) {\
		if (PQ_j < 2) (pq)[PQ_j].pqFromInt = NULL;\
		else (pq)[PQ_j].pqFromInt = &(pq)[PQ_j>>1];\
		(pq)[PQ_j].pqFromExt = &(pq)[(PQ_j+(n))>>1];\
		}\
	}


#define PQ_BUILD(pq,n,q)\
{\
	for (PQ_j=(n)-1;PQ_j>0;--PQ_j) {\
		PQ_i = (PQ_j<<1);\
		if (PQ_i < (n)) PQ_t = (pq)[PQ_i].pqWinner;\
		else PQ_t = &(pq)[PQ_i-(n)];\
		++PQ_i;\
		if (PQ_i < (n)) PQ_lt = (pq)[PQ_i].pqWinner;\
		else PQ_lt = &(pq)[PQ_i-(n)];\
		if (PQ_t->fKey < PQ_lt->fKey) {\
			(pq)[PQ_j].pqLoser = PQ_t;\
			(pq)[PQ_j].pqWinner = PQ_lt;\
			}\
		else {\
			(pq)[PQ_j].pqLoser = PQ_lt;\
			(pq)[PQ_j].pqWinner = PQ_t;\
			}\
		}\
	(q) = (pq)[1].pqWinner;\
	}


#define PQ_REPLACE(q)\
{\
	PQ_t = (q)->pqFromExt;\
	while (PQ_t) {\
		if (PQ_t->pqLoser->fKey > (q)->fKey) {\
			PQ_lt = PQ_t->pqLoser;\
			PQ_t->pqLoser = (q);\
			(q) = PQ_lt;\
			}\
		PQ_t = PQ_t->pqFromInt;\
		}\
	}



int smInit(SMX *,KD,int,float *);
void smFinish(SMX);
void smBallSearch(SMX,float,float *);
int  smBallGather(SMX,float,float *);
void smSmooth(SMX,void (*)(SMX,int,int,int *,float *));
void smReSmooth(SMX,void (*)(SMX,int,int,int *,float *));
void smDensity(SMX,int,int,int *,float *);
void smDensitySym(SMX,int,int,int *,float *);
void smMeanVel(SMX,int,int,int *,float *);
void smMeanVelSym(SMX,int,int,int *,float *);
void smVelDisp(SMX,int,int,int *,float *);
void smVelDispSym(SMX,int,int,int *,float *);
void smNull(SMX,int,int,int *,float *);
void smOutDensity(SMX,FILE *);
void smOutMeanVel(SMX,FILE *);
void smOutVelDisp(SMX,FILE *);
void smOutPhase(SMX,FILE *);
void smOutMach(SMX,FILE *);
void smOutSpeed(SMX,FILE *);

#endif



