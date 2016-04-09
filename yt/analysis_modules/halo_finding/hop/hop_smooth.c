/* SMOOTH.C */
/* This was written by Joachim Stadel and the NASA HPCC ESS at
the University of Washington Department of Astronomy as part of
the SMOOTH program, v2.0.1.
URL: http://www-hpcc.astro.washington.edu/tools/SMOOTH */
 
/* DJE--I have removed unneeded subroutines, notably those having
to do with velocity field reconstructions (because they refer to
particle data that I chose not to store) and output routines
(because I wanted binary output).  Also, the density subroutine
was slightly customized to reduce memory consumption in
the case of equal mass particles. */
 
/* HOP Version 1.0 (12/15/97) -- Original Release */
 
#include <stdio.h>
#include <stdlib.h>
#ifdef _WIN32
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#include <assert.h>
#include "smooth.h"
#include "kd.h"
#include "hop_numpy.h"

#define ISYM "d"
#define GSYM "g"
 
//#include "macros_and_parameters.h"
 
#define IMARK 1		/* All particles are marked to be included */
 
int smInit(SMX *psmx,KD kd,int nSmooth,float *fPeriod)
{
	SMX smx;
	PQ_STATIC;
	int pi,j;
    fprintf(stderr,"nSmooth = %d kd->nActive = %d\n", nSmooth, kd->nActive);
	assert(nSmooth <= kd->nActive);
	smx = (SMX)malloc(sizeof(struct smContext));
	assert(smx != NULL);
    smx->kd = NULL;
    
	smx->kd = kd;
	smx->nSmooth = nSmooth;
	smx->pq = (PQ *)malloc(nSmooth*sizeof(PQ));
	assert(smx->pq != NULL);
	PQ_INIT(smx->pq,nSmooth);
	smx->pfBall2 = (float *)malloc((kd->nActive+1)*sizeof(int));
	assert(smx->pfBall2 != NULL);
	smx->iMark = (char *)malloc(kd->nActive*sizeof(char));
	assert(smx->iMark);
	smx->nListSize = smx->nSmooth+RESMOOTH_SAFE;
	smx->fList = (float *)malloc(smx->nListSize*sizeof(float));
	assert(smx->fList != NULL);
	smx->pList = (int *)malloc(smx->nListSize*sizeof(int));
	assert(smx->pList != NULL);
	/*
	 ** Set for Periodic Boundary Conditions.
	 */
	for (j=0;j<3;++j) smx->fPeriod[j] = fPeriod[j];
	/*
	 ** Initialize arrays for calculated quantities.--DJE
	 */
	for (pi=0;pi<smx->kd->nActive;++pi) {
        NP_DENS(smx->kd, pi) = 0.0;
		smx->kd->p[pi].iHop = 0;
		}
	*psmx = smx;	
	return(1);
	}
 
 
void smFinish(SMX smx)
{
	free(smx->pfBall2);
	free(smx->iMark);
	free(smx->pq);
	free(smx);
	}
 
 
void smBallSearch(SMX smx,float fBall2,float *ri)
{
	KDN *c;
	PARTICLE *p;
	int cell,cp,ct,pj;
	float fDist2,dx,dy,dz,lx,ly,lz,sx,sy,sz,x,y,z;
	PQ *pq;
	PQ_STATIC;
 
	c = smx->kd->kdNodes;
	p = smx->kd->p;
	pq = smx->pqHead;
	x = ri[0];
	y = ri[1];
	z = ri[2];
	lx = smx->fPeriod[0];
	ly = smx->fPeriod[1];
	lz = smx->fPeriod[2];
	cell = ROOT;
	/*
	 ** First find the "local" Bucket.
	 ** This could mearly be the closest bucket to ri[3].
	 */
	while (cell < smx->kd->nSplit) {
		if (ri[c[cell].iDim] < c[cell].fSplit) cell = LOWER(cell);
		else cell = UPPER(cell);
		}
	/*
	 ** Now start the search from the bucket given by cell!
	 */
	for (pj=c[cell].pLower;pj<=c[cell].pUpper;++pj) {
		dx = x - NP_POS(smx->kd, pj, 0);
		dy = y - NP_POS(smx->kd, pj, 1);
		dz = z - NP_POS(smx->kd, pj, 2);
		fDist2 = dx*dx + dy*dy + dz*dz;
		if (fDist2 < fBall2) {
			if (smx->iMark[pj]) continue;
			smx->iMark[pq->p] = 0;
			smx->iMark[pj] = 1;
			pq->fKey = fDist2;
			pq->p = pj;
			pq->ax = 0.0;
			pq->ay = 0.0;
			pq->az = 0.0;
			PQ_REPLACE(pq);
			fBall2 = pq->fKey;
			}
		}
	while (cell != ROOT) {
		cp = SIBLING(cell);
		ct = cp;
		SETNEXT(ct);
		while (1) {
			INTERSECT(c,cp,fBall2,lx,ly,lz,x,y,z,sx,sy,sz);
			/*
			 ** We have an intersection to test.
			 */
			if (cp < smx->kd->nSplit) {
				cp = LOWER(cp);
				continue;
				}
			else {
				for (pj=c[cp].pLower;pj<=c[cp].pUpper;++pj) {
                    dx = sx - NP_POS(smx->kd, pj, 0);
                    dy = sy - NP_POS(smx->kd, pj, 1);
                    dz = sz - NP_POS(smx->kd, pj, 2);
					fDist2 = dx*dx + dy*dy + dz*dz;
					if (fDist2 < fBall2) {
						if (smx->iMark[pj]) continue;
						smx->iMark[pq->p] = 0;
						smx->iMark[pj] = 1;
						pq->fKey = fDist2;
						pq->p = pj;
						pq->ax = sx - x;
						pq->ay = sy - y;
						pq->az = sz - z;
						PQ_REPLACE(pq);
						fBall2 = pq->fKey;
						}
					}
				}
		GetNextCell:
			SETNEXT(cp);
			if (cp == ct) break;
			}
		cell = PARENT(cell);
		}
	smx->pqHead = pq;
	}
 
 
int smBallGather(SMX smx,float fBall2,float *ri)
{
	KDN *c;
	PARTICLE *p;
	int pj,nCnt,cp,nSplit;
	float dx,dy,dz,x,y,z,lx,ly,lz,sx,sy,sz,fDist2;
 
	c = smx->kd->kdNodes;
	p = smx->kd->p;
	nSplit = smx->kd->nSplit;
	lx = smx->fPeriod[0];
	ly = smx->fPeriod[1];
	lz = smx->fPeriod[2];
	x = ri[0];
	y = ri[1];
	z = ri[2];
	nCnt = 0;
	cp = ROOT;
	while (1) {
		INTERSECT(c,cp,fBall2,lx,ly,lz,x,y,z,sx,sy,sz);
		/*
		 ** We have an intersection to test.
		 */
		if (cp < nSplit) {
			cp = LOWER(cp);
			continue;
			}
		else {
			for (pj=c[cp].pLower;pj<=c[cp].pUpper;++pj) {
                dx = sx - NP_POS(smx->kd, pj, 0);
                dy = sy - NP_POS(smx->kd, pj, 1);
                dz = sz - NP_POS(smx->kd, pj, 2);
				fDist2 = dx*dx + dy*dy + dz*dz;
				if (fDist2 < fBall2) {
					smx->fList[nCnt] = fDist2;
					smx->pList[nCnt++] = pj;
					/* Insert debugging flag here */
					if (nCnt > smx->nListSize) {
					    fprintf(stderr,"nCnt too big.\n");
					    }
					}
				}
			}
	GetNextCell:
		SETNEXT(cp);
		if (cp == ROOT) break;
		}
	assert(nCnt <= smx->nListSize);
	return(nCnt);
	}
 
 
void smSmooth(SMX smx,void (*fncSmooth)(SMX,int,int,int *,float *))
{
	KDN *c;
	PARTICLE *p;
    PQ *pq,*pqLast;
	PQ_STATIC;
	int cell;
	int pi,pin,pj,pNext,nCnt,nSmooth;
	float dx,dy,dz,x,y,z,h2,ax,ay,az;
    float temp_ri[3];
 
 
	for (pi=0;pi<smx->kd->nActive;++pi) {
		if (IMARK) smx->pfBall2[pi] = -1.0;
		else smx->pfBall2[pi] = 1.0;	/* pretend it is already done! */
		}
	smx->pfBall2[smx->kd->nActive] = -1.0; /* stop condition */
	for (pi=0;pi<smx->kd->nActive;++pi) {
		smx->iMark[pi] = 0;
		}
	pqLast = &smx->pq[smx->nSmooth-1];
	c = smx->kd->kdNodes;
	p = smx->kd->p;
	nSmooth = smx->nSmooth;
	/*
	 ** Initialize Priority Queue.
	 */
	pin = 0;
	pNext = 1;
	ax = 0.0;
	ay = 0.0;
	az = 0.0;
	for (pq=smx->pq,pj=0;pq<=pqLast;++pq,++pj) {
		smx->iMark[pj] = 1;
		pq->p = pj;
		pq->ax = ax;
		pq->ay = ay;
		pq->az = az;
		}
	while (1) {
		if (smx->pfBall2[pin] >= 0) {
			/*
			 ** Find next particle which is not done, and load the
			 ** priority queue with nSmooth number of particles.
			 */
			while (smx->pfBall2[pNext] >= 0) ++pNext;
			/*
			 ** Check if we are really finished.
			 */
			if (pNext == smx->kd->nActive) break;
			pi = pNext;
			++pNext;
			x = NP_POS(smx->kd, pi, 0);
			y = NP_POS(smx->kd, pi, 1);
			z = NP_POS(smx->kd, pi, 2);
			/* printf("%"ISYM": %"GSYM" %"GSYM" %"GSYM"\n", pi, x, y, z); */
			/*
			 ** First find the "local" Bucket.
			 ** This could mearly be the closest bucket to ri[3].
			 */
			cell = ROOT;
			while (cell < smx->kd->nSplit) {
                if (NP_POS(smx->kd, pi, c[cell].iDim) < c[cell].fSplit)
					cell = LOWER(cell);
				else
					cell = UPPER(cell);
				}
			/*
			 ** Remove everything from the queue.
			 */
			smx->pqHead = NULL;
			for (pq=smx->pq;pq<=pqLast;++pq) smx->iMark[pq->p] = 0;
			/*
			 ** Add everything from pj up to and including pj+nSmooth-1.
			 */
			pj = c[cell].pLower;
			if (pj > smx->kd->nActive - nSmooth)
				pj = smx->kd->nActive - nSmooth;
			for (pq=smx->pq;pq<=pqLast;++pq) {
				smx->iMark[pj] = 1;
                dx = x - NP_POS(smx->kd, pj, 0);
                dy = y - NP_POS(smx->kd, pj, 1);
                dz = z - NP_POS(smx->kd, pj, 2);
				pq->fKey = dx*dx + dy*dy + dz*dz;
				pq->p = pj++;
				pq->ax = 0.0;
				pq->ay = 0.0;
				pq->az = 0.0;
				}
			PQ_BUILD(smx->pq,nSmooth,smx->pqHead);
			}
		else {
			/*
			 ** Calculate the priority queue using the previous particles!
			 */
			pi = pin;
			x = NP_POS(smx->kd, pi, 0);
			y = NP_POS(smx->kd, pi, 1);
			z = NP_POS(smx->kd, pi, 2);
			smx->pqHead = NULL;
			for (pq=smx->pq;pq<=pqLast;++pq) {
				pq->ax -= ax;
				pq->ay -= ay;
				pq->az -= az;
				dx = x + pq->ax - NP_POS(smx->kd, pq->p, 0);
				dy = y + pq->ay - NP_POS(smx->kd, pq->p, 1);
				dz = z + pq->az - NP_POS(smx->kd, pq->p, 2);
				pq->fKey = dx*dx + dy*dy + dz*dz;
				}
			PQ_BUILD(smx->pq,nSmooth,smx->pqHead);
			ax = 0.0;
			ay = 0.0;
			az = 0.0;
			}
        temp_ri[0] = NP_POS(smx->kd, pi, 0);
        temp_ri[1] = NP_POS(smx->kd, pi, 1);
        temp_ri[2] = NP_POS(smx->kd, pi, 2);
		smBallSearch(smx,smx->pqHead->fKey,temp_ri);
		smx->pfBall2[pi] = smx->pqHead->fKey;
		/*
		 ** Pick next particle, 'pin'.
		 ** Create fList and pList for function 'fncSmooth'.
		 */
		pin = pi;
		nCnt = 0;
		h2 = smx->pqHead->fKey;
		for (pq=smx->pq;pq<=pqLast;++pq) {
			if (pq == smx->pqHead) continue;
			smx->pList[nCnt] = pq->p;
			smx->fList[nCnt++] = pq->fKey;
			if (smx->pfBall2[pq->p] >= 0) continue;
			if (pq->fKey < h2) {
				pin = pq->p;
				h2 = pq->fKey;
				ax = pq->ax;
				ay = pq->ay;
				az = pq->az;
				}
			}
		(*fncSmooth)(smx,pi,nCnt,smx->pList,smx->fList);
		}
	}
 
 
void smReSmooth(SMX smx,void (*fncSmooth)(SMX,int,int,int *,float *))
{
	PARTICLE *p;
	int pi,nSmooth;
    float temp_ri[3];
 
	p = smx->kd->p;
	for (pi=0;pi<smx->kd->nActive;++pi) {
		if (IMARK == 0) continue;
		/*
		 ** Do a Ball Gather at the radius of the most distant particle
		 ** which is smDensity sets in smx->pBall[pi].
		 */
        temp_ri[0] = NP_POS(smx->kd, pi, 0);
        temp_ri[1] = NP_POS(smx->kd, pi, 1);
        temp_ri[2] = NP_POS(smx->kd, pi, 2);
		nSmooth = smBallGather(smx,smx->pfBall2[pi],temp_ri);
		(*fncSmooth)(smx,pi,nSmooth,smx->pList,smx->fList);
		}
 	}
 
 
void smDensity(SMX smx,int pi,int nSmooth,int *pList,float *fList)
{
	float ih2,r2,rs,fDensity;
	int i,pj;
 
	ih2 = 4.0/smx->pfBall2[pi];
	fDensity = 0.0;
	for (i=0;i<nSmooth;++i) {
		pj = pList[i];
		r2 = fList[i]*ih2;
		rs = 2.0 - sqrt(r2);
		if (r2 < 1.0) rs = (1.0 - 0.75*rs*r2);
		else rs = 0.25*rs*rs*rs;
#ifdef DIFFERENT_MASSES
        fDensity += rs*NP_MASS(smx->kd, pj);
#else
		fDensity += rs*smx->kd->fMass;
#endif
		}
    NP_DENS(smx->kd, pi) = M_1_PI*sqrt(ih2)*ih2*fDensity;
	}
 
 
void smDensitySym(SMX smx,int pi,int nSmooth,int *pList,float *fList)
{
	float fNorm,ih2,r2,rs;
	int i,pj;
 
	ih2 = 4.0/smx->pfBall2[pi];
	fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
	for (i=0;i<nSmooth;++i) {
		pj = pList[i];
		r2 = fList[i]*ih2;
		rs = 2.0 - sqrt(r2);
		if (r2 < 1.0) rs = (1.0 - 0.75*rs*r2);
		else rs = 0.25*rs*rs*rs;
		rs *= fNorm;
#ifdef DIFFERENT_MASSES
        NP_DENS(smx->kd, pi) += rs*NP_MASS(smx->kd, pj);
        NP_DENS(smx->kd, pj) += rs*NP_MASS(smx->kd, pi);
#else
		smx->kd->p[pi].fDensity += rs*smx->kd->fMass;
		smx->kd->p[pj].fDensity += rs*smx->kd->fMass;
#endif
		}
	}
 
/* I'm not using the following function, but I left it here in case someone
wants the densities outputted in Tipsy format.  But you're probably better
off just fetching the smooth() program from the HPCC web site... */
 
void smOutDensity(SMX smx,FILE *fp)
{
  int i,iCnt;

  fprintf(fp,"%"ISYM"\n",smx->kd->nParticles);
  iCnt = 0;
  for (i=0;i<smx->kd->nGas;++i) {
    if (smx->kd->bGas) {
      if (IMARK)
        fprintf(fp,"%.8"GSYM"\n",NP_DENS(smx->kd, iCnt));
      else fprintf(fp,"0\n");
      ++iCnt;
    }
    else fprintf(fp,"0\n");
  }
  for (i=0;i<smx->kd->nDark;++i) {
    if (smx->kd->bDark) {
      if (IMARK)
        fprintf(fp,"%.8"GSYM"\n",NP_DENS(smx->kd, iCnt));
      else fprintf(fp,"0\n");
      ++iCnt;
    }
    else fprintf(fp,"0\n");
  }
  for (i=0;i<smx->kd->nStar;++i) {
    if (smx->kd->bStar) {
      if (IMARK)
        fprintf(fp,"%.8"GSYM"\n",NP_DENS(smx->kd, iCnt));
      else fprintf(fp,"0\n");
      ++iCnt;
    }
    else fprintf(fp,"0\n");
  }
}


