#ifndef KD_HINCLUDED
#define KD_HINCLUDED

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

#define KD_ORDERTEMP	256

typedef struct Particle {
	float r[3];
	int iGroup;
	int iOrder;
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
	float fPeriod[3];
	int nLevels;
	int nNodes;
	int nSplit;
	PARTICLE *p;
	KDN *kdNodes;
	int nGroup;
	int uSecond;
	int uMicro;
	} * KD;


#define INTERSECT(c,cp,fBall2,lx,ly,lz,x,y,z,sx,sy,sz)\
{\
	float dx,dy,dz,dx1,dy1,dz1,fDist2,fMax2;\
	dx = c[cp].bnd.fMin[0]-x;\
	dx1 = x-c[cp].bnd.fMax[0];\
	dy = c[cp].bnd.fMin[1]-y;\
	dy1 = y-c[cp].bnd.fMax[1];\
	dz = c[cp].bnd.fMin[2]-z;\
	dz1 = z-c[cp].bnd.fMax[2];\
	if (dx > 0.0) {\
		if (dx1+lx < dx) {\
			dx1 += lx;\
			dx -= lx;\
			sx = x+lx;\
			fDist2 = dx1*dx1;\
			fMax2 = dx*dx;\
			}\
		else {\
			sx = x;\
			fDist2 = dx*dx;\
			fMax2 = dx1*dx1;\
			}\
		if (fDist2 > fBall2) goto GetNextCell;\
		}\
	else if (dx1 > 0.0) {\
		if (dx+lx < dx1) {\
		    dx += lx;\
			dx1 -= lx;\
			sx = x-lx;\
			fDist2 = dx*dx;\
			fMax2 = dx1*dx1;\
			}\
		else {\
			sx = x;\
			fDist2 = dx1*dx1;\
			fMax2 = dx*dx;\
			}\
		if (fDist2 > fBall2) goto GetNextCell;\
		}\
	else {\
		sx = x;\
		fDist2 = 0.0;\
		if (dx < dx1) fMax2 = dx*dx;\
		else fMax2 = dx1*dx1;\
		}\
	if (dy > 0.0) {\
		if (dy1+ly < dy) {\
		    dy1 += ly;\
			dy -= ly;\
			sy = y+ly;\
			fDist2 += dy1*dy1;\
			fMax2 += dy*dy;\
			}\
		else {\
			sy = y;\
			fDist2 += dy*dy;\
			fMax2 += dy1*dy1;\
			}\
		if (fDist2 > fBall2) goto GetNextCell;\
		}\
	else if (dy1 > 0.0) {\
		if (dy+ly < dy1) {\
		    dy += ly;\
			dy1 -= ly;\
			sy = y-ly;\
			fDist2 += dy*dy;\
			fMax2 += dy1*dy1;\
			}\
		else {\
			sy = y;\
			fDist2 += dy1*dy1;\
			fMax2 += dy*dy;\
			}\
		if (fDist2 > fBall2) goto GetNextCell;\
		}\
	else {\
		sy = y;\
		if (dy < dy1) fMax2 += dy*dy;\
		else fMax2 += dy1*dy1;\
		}\
	if (dz > 0.0) {\
		if (dz1+lz < dz) {\
		    dz1 += lz;\
            dz -= lz;\
			sz = z+lz;\
			fDist2 += dz1*dz1;\
			fMax2 += dz*dz;\
			}\
		else {\
			sz = z;\
			fDist2 += dz*dz;\
			fMax2 += dz1*dz1;\
			}\
		if (fDist2 > fBall2) goto GetNextCell;\
		}\
	else if (dz1 > 0.0) {\
		if (dz+lz < dz1) {\
			dz += lz;\
		    dz1 -= lz;\
			sz = z-lz;\
			fDist2 += dz*dz;\
			fMax2 += dz1*dz1;\
			}\
		else {\
			sz = z;\
			fDist2 += dz1*dz1;\
			fMax2 += dz*dz;\
			}\
		if (fDist2 > fBall2) goto GetNextCell;\
		}\
	else {\
		sz = z;\
		if (dz < dz1) fMax2 += dz*dz;\
		else fMax2 += dz1*dz1;\
		}\
	if (fMax2 < fBall2) goto ContainedCell;\
	}


void kdTime(KD,int *,int *);
int kdInit(KD *,int,float *);
void kdReadTipsy(KD,FILE *,int,int,int);
void kdBuildTree(KD);
int kdFoF(KD,float);
int kdTooSmall(KD,int);
void kdOrder(KD);
void kdOutGroup(KD,char *);
void kdFinish(KD);

#endif











