/* KD.C */
/* This was written by Joachim Stadel and the NASA HPCC ESS at
the University of Washington Department of Astronomy as part of
the SMOOTH program, v2.0.1.
URL: http://www-hpcc.astro.washington.edu/tools/SMOOTH */
 
/* DJE--I have removed all the subroutines not used by HOP, notably
the input and output routines. */
 
/* HOP Version 1.0 (12/15/97) -- Original Release */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <assert.h>
#include "kd.h"
//#include "macros_and_parameters.h"
/* #include "tipsydefs.h" */ /* Don't need this, since I removed kdReadTipsy()*/
 
 
#define MAX_ROOT_ITTR	32
 
 
void kdTime(KD kd,int *puSecond,int *puMicro)
{
	struct rusage ru;
 
	getrusage(0,&ru);
	*puMicro = ru.ru_utime.tv_usec - kd->uMicro;
	*puSecond = ru.ru_utime.tv_sec - kd->uSecond;
	if (*puMicro < 0) {
		*puMicro += 1000000;
		*puSecond -= 1;
		}
	kd->uSecond = ru.ru_utime.tv_sec;
	kd->uMicro = ru.ru_utime.tv_usec;
	}
 
 
int kdInit(KD *pkd,int nBucket)
{
	KD kd;
 
	kd = (KD)malloc(sizeof(struct kdContext));
	assert(kd != NULL);
	kd->nBucket = nBucket;
	*pkd = kd;
	return(1);
	}
 
/*
 ** JST's Median Algorithm
 */
int kdMedianJst(KD kd,int d,int l,int u)
{
	float fm;
    int i,k,m;
    PARTICLE *p,t;
 
	p = kd->p;
    k = (l+u)/2;
	m = k;
    while (l < u) {
		m = (l+u)/2;
		fm = p[m].r[d];
		t = p[m];
		p[m] = p[u];
		p[u] = t;
		i = u-1;
		m = l;
		while (p[m].r[d] < fm) ++m;
		while (m < i) {
			while (p[i].r[d] >= fm) if (--i == m) break;
			/*
			 ** Swap
			 */
			t = p[m];
			p[m] = p[i];
			p[i] = t;
			--i;
			while (p[m].r[d] < fm) ++m;
			}
		t = p[m];
		p[m] = p[u];
		p[u] = t;
        if (k <= m) u = m-1;
        if (k >= m) l = m+1;
        }
    return(m);
    }
 
 
void kdCombine(KDN *p1,KDN *p2,KDN *pOut)
{
	int j;
 
	/*
	 ** Combine the bounds.
	 */
	for (j=0;j<3;++j) {
		if (p2->bnd.fMin[j] < p1->bnd.fMin[j])
			pOut->bnd.fMin[j] = p2->bnd.fMin[j];
		else
			pOut->bnd.fMin[j] = p1->bnd.fMin[j];
		if (p2->bnd.fMax[j] > p1->bnd.fMax[j])
			pOut->bnd.fMax[j] = p2->bnd.fMax[j];
		else
			pOut->bnd.fMax[j] = p1->bnd.fMax[j];
		}
	}
 
 
void kdUpPass(KD kd,int iCell)
{
	KDN *c;
	int l,u,pj,j;
 
	c = kd->kdNodes;
	if (c[iCell].iDim != -1) {
		l = LOWER(iCell);
		u = UPPER(iCell);
		kdUpPass(kd,l);
		kdUpPass(kd,u);
		kdCombine(&c[l],&c[u],&c[iCell]);
		}
	else {
		l = c[iCell].pLower;
		u = c[iCell].pUpper;
		for (j=0;j<3;++j) {
			c[iCell].bnd.fMin[j] = kd->p[u].r[j];
			c[iCell].bnd.fMax[j] = kd->p[u].r[j];
			}
		for (pj=l;pj<u;++pj) {
			for (j=0;j<3;++j) {
				if (kd->p[pj].r[j] < c[iCell].bnd.fMin[j])
					c[iCell].bnd.fMin[j] = kd->p[pj].r[j];
				if (kd->p[pj].r[j] > c[iCell].bnd.fMax[j])
					c[iCell].bnd.fMax[j] = kd->p[pj].r[j];
				}
			}
		}
	}
 
int kdBuildTree(KD kd)
{
	int l,n,i,d,m,j,ct;
	KDN *c;
 
	n = kd->nActive;
	kd->nLevels = 1;
	l = 1;
	while (n > kd->nBucket) {
		n = n>>1;
		l = l<<1;
		++kd->nLevels;
		}
	kd->nSplit = l;
	kd->nNodes = l<<1;
	kd->kdNodes = (KDN *)malloc(kd->nNodes*sizeof(KDN));
	assert(kd->kdNodes != NULL);
	/*
	 ** Set up ROOT node
	 */
	c = kd->kdNodes;
	c[ROOT].pLower = 0;
	c[ROOT].pUpper = kd->nActive-1;
	c[ROOT].bnd = kd->bnd;
	i = ROOT;
	ct = ROOT;
	SETNEXT(ct);
	while (1) {
		if (i < kd->nSplit) {
			d = 0;
			for (j=1;j<3;++j) {
				if (c[i].bnd.fMax[j]-c[i].bnd.fMin[j] >
					c[i].bnd.fMax[d]-c[i].bnd.fMin[d]) d = j;
				}
			c[i].iDim = d;
			m = kdMedianJst(kd,d,c[i].pLower,c[i].pUpper);
			c[i].fSplit = kd->p[m].r[d];
			c[LOWER(i)].bnd = c[i].bnd;
			c[LOWER(i)].bnd.fMax[d] = c[i].fSplit;
			c[LOWER(i)].pLower = c[i].pLower;
			c[LOWER(i)].pUpper = m-1;
			c[UPPER(i)].bnd = c[i].bnd;
			c[UPPER(i)].bnd.fMin[d] = c[i].fSplit;
			c[UPPER(i)].pLower = m;
			c[UPPER(i)].pUpper = c[i].pUpper;
			i = LOWER(i);
			}
		else {
			c[i].iDim = -1;
			SETNEXT(i);
			if (i == ct) break;
			}
		}
	kdUpPass(kd,ROOT);
	return(1);
	}
 
 
int cmpParticles(const void *v1,const void *v2)
{
	PARTICLE *p1=(PARTICLE *)v1,*p2=(PARTICLE *)v2;
	
	return(p1->iOrder - p2->iOrder);
	}
 
 
void kdOrder(KD kd)
{
	qsort(kd->p,kd->nActive,sizeof(PARTICLE),cmpParticles);
	}
 
void kdFinish(KD kd)
{
	free(kd->p);
	free(kd->kdNodes);
	free(kd);
	}
 
