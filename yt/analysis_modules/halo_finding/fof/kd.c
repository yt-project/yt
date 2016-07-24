#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <sys/resource.h>
#include <sys/time.h>
#endif
#include <assert.h>
#include "kd.h"
#include "tipsydefs.h"


void kdTimeFoF(KDFOF kd,int *puSecond,int *puMicro)
{

#ifdef _WIN32
        int secs, usecs;
        HANDLE hProcess = GetCurrentProcess();
	FILETIME ftCreation, ftExit, ftKernel, ftUser;
	SYSTEMTIME stUser;
	GetProcessTimes(hProcess, &ftCreation, &ftExit, 
			&ftKernel, &ftUser);
	FileTimeToSystemTime(&ftUser, &stUser);
	secs = (int)((double)stUser.wHour*3600.0 +
			  (double)stUser.wMinute*60.0 +
			  (double)stUser.wSecond);
	usecs = (int)((double)stUser.wMilliseconds/1000.0);
	*puMicro = usecs;
	*puSecond = secs;
	if (*puMicro < 0) {
	  *puMicro += 1000000;
	  *puSecond -= 1;
	}
	kd->uSecond = secs;
	kd->uMicro = usecs;
#else
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
#endif
}

int kdInitFoF(KDFOF *pkd,int nBucket,float *fPeriod)
{
	KDFOF kd;
	int j;

	kd = (KDFOF)malloc(sizeof(struct kdContext));
	assert(kd != NULL);
	kd->nBucket = nBucket;
	for (j=0;j<3;++j) kd->fPeriod[j] = fPeriod[j];
	kd->p = NULL;
	kd->kdNodes = NULL;
	*pkd = kd;
	return(1);
	}


void kdReadTipsyFoF(KDFOF kd,FILE *fp,int bDark,int bGas,int bStar)
{
	int i,j,nCnt;
	struct dump h;
	struct gas_particle gp;
	struct dark_particle dp;
	struct star_particle sp;

	fread(&h,sizeof(struct dump),1,fp);
	kd->nParticles = h.nbodies;
	kd->nDark = h.ndark;
	kd->nGas = h.nsph;
	kd->nStar = h.nstar;
	kd->fTime = h.time;
	kd->nActive = 0;
	if (bDark) kd->nActive += kd->nDark;
	if (bGas) kd->nActive += kd->nGas;
	if (bStar) kd->nActive += kd->nStar;
	kd->bDark = bDark;
	kd->bGas = bGas;
	kd->bStar = bStar;
	/*
	 ** Allocate particles.
	 */
	kd->p = (PARTICLEFOF *)malloc(kd->nActive*sizeof(PARTICLEFOF));
	assert(kd->p != NULL);
	/*
	 ** Read Stuff!
	 */
	nCnt = 0;
	for (i=0;i<h.nsph;++i) {
		fread(&gp,sizeof(struct gas_particle),1,fp);
		if (bGas) {
			kd->p[nCnt].iOrder = nCnt;
			for (j=0;j<3;++j) kd->p[nCnt].r[j] = gp.pos[j];
			++nCnt;
			}
		}
	for (i=0;i<h.ndark;++i) {
		fread(&dp,sizeof(struct dark_particle),1,fp);
		if (bDark) {
			kd->p[nCnt].iOrder = nCnt;
			for (j=0;j<3;++j) kd->p[nCnt].r[j] = dp.pos[j];
			++nCnt;
			}
		}
	for (i=0;i<h.nstar;++i) {
		fread(&sp,sizeof(struct star_particle),1,fp);
		if (bStar) {
			kd->p[nCnt].iOrder = nCnt;
			for (j=0;j<3;++j) kd->p[nCnt].r[j] = sp.pos[j];
			++nCnt;
			}
		}
	}


void kdSelectFoF(KDFOF kd,int d,int k,int l,int r)
{
	PARTICLEFOF *p,t;
	double v;
	int i,j;

	p = kd->p;
	while (r > l) {
		v = p[k].r[d];
		t = p[r];
		p[r] = p[k];
		p[k] = t;
		i = l - 1;
		j = r;
		while (1) {
			while (i < j) if (p[++i].r[d] >= v) break;
			while (i < j) if (p[--j].r[d] <= v) break;
			t = p[i];
			p[i] = p[j];
			p[j] = t;
			if (j <= i) break;
			}
		p[j] = p[i];
		p[i] = p[r];
		p[r] = t;
		if (i >= k) r = i - 1;
		if (i <= k) l = i + 1;
		}
	}


void kdCombineFoF(KDNFOF *p1,KDNFOF *p2,KDNFOF *pOut)
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


void kdUpPassFoF(KDFOF kd,int iCell)
{
	KDNFOF *c;
	int l,u,pj,j;

	c = kd->kdNodes;
	if (c[iCell].iDim != -1) {
		l = LOWERFOF(iCell);
		u = UPPERFOF(iCell);
		kdUpPassFoF(kd,l);
		kdUpPassFoF(kd,u);
		kdCombineFoF(&c[l],&c[u],&c[iCell]);
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

void kdBuildTreeFoF(KDFOF kd)
{
	int l,n,i,d,m,j,diff;
	KDNFOF *c;
	BNDFOF bnd;

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
	if (kd->kdNodes != NULL) free(kd->kdNodes);
	kd->kdNodes = (KDNFOF *)malloc(kd->nNodes*sizeof(KDNFOF));
	assert(kd->kdNodes != NULL);
	/*
	 ** Calculate Bounds.
	 */
	for (j=0;j<3;++j) {
		bnd.fMin[j] = kd->p[0].r[j];
		bnd.fMax[j] = kd->p[0].r[j];
		}
	for (i=1;i<kd->nActive;++i) {
		for (j=0;j<3;++j) {
			if (bnd.fMin[j] > kd->p[i].r[j]) 
				bnd.fMin[j] = kd->p[i].r[j];
			else if (bnd.fMax[j] < kd->p[i].r[j])
				bnd.fMax[j] = kd->p[i].r[j];
			}
		}
	/*
	 ** Set up ROOTFOF node
	 */
	c = kd->kdNodes;
	c[ROOTFOF].pLower = 0;
	c[ROOTFOF].pUpper = kd->nActive-1;
	c[ROOTFOF].bnd = bnd;
	i = ROOTFOF;
	while (1) {
		assert(c[i].pUpper - c[i].pLower + 1 > 0);
		if (i < kd->nSplit && (c[i].pUpper - c[i].pLower) > 0) {
			d = 0;
			for (j=1;j<3;++j) {
				if (c[i].bnd.fMax[j]-c[i].bnd.fMin[j] > 
					c[i].bnd.fMax[d]-c[i].bnd.fMin[d]) d = j;
				}
			c[i].iDim = d;

			m = (c[i].pLower + c[i].pUpper)/2;
			kdSelectFoF(kd,d,m,c[i].pLower,c[i].pUpper);

			c[i].fSplit = kd->p[m].r[d];
			c[LOWERFOF(i)].bnd = c[i].bnd;
			c[LOWERFOF(i)].bnd.fMax[d] = c[i].fSplit;
			c[LOWERFOF(i)].pLower = c[i].pLower;
			c[LOWERFOF(i)].pUpper = m;
			c[UPPERFOF(i)].bnd = c[i].bnd;
			c[UPPERFOF(i)].bnd.fMin[d] = c[i].fSplit;
			c[UPPERFOF(i)].pLower = m+1;
			c[UPPERFOF(i)].pUpper = c[i].pUpper;
			diff = (m-c[i].pLower+1)-(c[i].pUpper-m);
			assert(diff == 0 || diff == 1);
			i = LOWERFOF(i);
			}
		else {
			c[i].iDim = -1;
			SETNEXTFOF(i);
			if (i == ROOTFOF) break;
			}
		}
	kdUpPassFoF(kd,ROOTFOF);
	}


int kdFoF(KDFOF kd,float fEps)
{
	PARTICLEFOF *p;
	KDNFOF *c;
	int pi,pj,pn,cp;

	int iGroup;

	int *Fifo,iHead,iTail,nFifo;
	float fEps2;
	float dx,dy,dz,x,y,z,lx,ly,lz,sx,sy,sz,fDist2;

	p = kd->p;
	c = kd->kdNodes;
	lx = kd->fPeriod[0];
	ly = kd->fPeriod[1];
	lz = kd->fPeriod[2];
	fEps2 = fEps*fEps;
	for (pn=0;pn<kd->nActive;++pn) p[pn].iGroup = 0;
	nFifo = kd->nActive;
	Fifo = (int *)malloc(nFifo*sizeof(int));
	assert(Fifo != NULL);
	iHead = 0;
	iTail = 0;
	iGroup = 0;
	for (pn=0;pn<kd->nActive;++pn) {
		if (p[pn].iGroup) continue;
		++iGroup;
		/*
		 ** Mark it and add to the do-fifo.
		 */
		p[pn].iGroup = iGroup;
		Fifo[iTail++] = pn;
		if (iTail == nFifo) iTail = 0;
		while (iHead != iTail) {
			pi = Fifo[iHead++];
			if (iHead == nFifo) iHead = 0;
			/*
			 ** Now do an fEps-Ball Gather!
			 */
			x = p[pi].r[0];
			y = p[pi].r[1];
			z = p[pi].r[2];
			cp = ROOTFOF;
			while (1) {
				INTERSECTFOF(c,cp,fEps2,lx,ly,lz,x,y,z,sx,sy,sz);
				/*
				 ** We have an intersection to test.
				 */
				if (c[cp].iDim >= 0) {
					cp = LOWERFOF(cp);
					continue;
					}
				else {
					for (pj=c[cp].pLower;pj<=c[cp].pUpper;++pj) {
						if (p[pj].iGroup) continue;
						dx = sx - p[pj].r[0];
						dy = sy - p[pj].r[1];
						dz = sz - p[pj].r[2];
						fDist2 = dx*dx + dy*dy + dz*dz;
						if (fDist2 < fEps2) {
							/*
							 ** Mark it and add to the do-fifo.
							 */
							p[pj].iGroup = iGroup;
							Fifo[iTail++] = pj;
							if (iTail == nFifo) iTail = 0;
							}
						}
					SETNEXTFOF(cp);
					if (cp == ROOTFOF) break;
					continue;
					}
			ContainedCell:
				for (pj=c[cp].pLower;pj<=c[cp].pUpper;++pj) {
					if (p[pj].iGroup) continue;
					/*
					 ** Mark it and add to the do-fifo.
					 */
					p[pj].iGroup = iGroup;
					Fifo[iTail++] = pj;
					if (iTail == nFifo) iTail = 0;
					}
			GetNextCell:
				SETNEXTFOF(cp);
				if (cp == ROOTFOF) break;
				}
			}
		}
	free(Fifo);
	kd->nGroup = iGroup+1;
	return(kd->nGroup-1);
	}


int kdTooSmallFoF(KDFOF kd,int nMembers)
{
	int *pnMembers,*pMap;
	int i,pi,nGroup;

	pnMembers = (int *)malloc(kd->nGroup*sizeof(int));
	assert(pnMembers != NULL);
	pMap = (int *)malloc(kd->nGroup*sizeof(int));
	assert(pMap != NULL);
	for (i=0;i<kd->nGroup;++i) pnMembers[i] = 0;
	for (pi=0;pi<kd->nActive;++pi) {
		++pnMembers[kd->p[pi].iGroup];
		}
	for (i=1;i<kd->nGroup;++i) {
		if (pnMembers[i] < nMembers) {
			pnMembers[i] = 0;
			}
		}
	/*
	 ** Create a remapping!
	 */
	pMap[0] = 0;
	nGroup = 1;
	for (i=1;i<kd->nGroup;++i) {
		pMap[i] = nGroup;
		if (pnMembers[i] == 0) {
			pMap[i] = -1; /* was 0 */
			}
		else {
			++nGroup;
			}
		}
	/*
	 ** Remap the groups.
	 */
	for (pi=0;pi<kd->nActive;++pi) {
		kd->p[pi].iGroup = pMap[kd->p[pi].iGroup];
		}
	free(pMap);
	free(pnMembers);
	kd->nGroup = nGroup;
	return(nGroup-1);
	}


int CmpParticlesFoF(const void *v1,const void *v2)
{
	PARTICLEFOF *p1 = (PARTICLEFOF *)v1;
	PARTICLEFOF *p2 = (PARTICLEFOF *)v2;
	return(p1->iOrder - p2->iOrder);
	}

void kdOrderFoF(KDFOF kd)
{
	qsort(kd->p,kd->nActive,sizeof(PARTICLEFOF),CmpParticlesFoF);
	}


void kdOutGroupFoF(KDFOF kd,char *pszFile)
{
	FILE *fp;
	int i,iCnt;

	fp = fopen(pszFile,"w");
	assert(fp != NULL);
	fprintf(fp,"%d\n",kd->nParticles);
	iCnt = 0;
	for (i=0;i<kd->nGas;++i) {
		if (kd->bGas) fprintf(fp,"%d\n",kd->p[iCnt++].iGroup);
		else fprintf(fp,"0\n");
		}
	for (i=0;i<kd->nDark;++i) {
		if (kd->bDark) fprintf(fp,"%d\n",kd->p[iCnt++].iGroup);
		else fprintf(fp,"0\n");
		}
	for (i=0;i<kd->nStar;++i) {
		if (kd->bStar) fprintf(fp,"%d\n",kd->p[iCnt++].iGroup);
		else fprintf(fp,"0\n");
		}
	fclose(fp);
	}


void kdFinishFoF(KDFOF kd)
{
	free(kd->p);
	free(kd->kdNodes);
	free(kd);
	}

