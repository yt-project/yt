#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include "kd.h"


void usage(void)
{
	fprintf(stderr,"USAGE:\n");
	fprintf(stderr,"fof -e <Linking Length>\n");
	fprintf(stderr,"   [-m <nMinMembers>] [-dgs] [-v]\n");
	fprintf(stderr,"   [-o <Output Name>] [-p <xyzPeriod>]\n");
	fprintf(stderr,"   [-px <xPeriod>] [-py <yPeriod>] [-pz <zPeriod>]\n");
	fprintf(stderr,"Input taken from stdin in tipsy binary format.\n");
	fprintf(stderr,"SEE MAN PAGE: fof(1) for more information.\n");
	exit(1);
	}

void main(int argc,char **argv)
{
	KD kd;
	int nBucket,i,j;
	char ach[80];
	float fPeriod[3],fEps;
	int bDark,bGas,bStar;
	int nMembers,nGroup,bVerbose;
	int sec,usec;
	char *p;
	
   	nBucket = 16;
	nMembers = 8;
	bDark = 1;
	bGas = 1;
	bStar = 1;
	bVerbose = 0;
	strcpy(ach,"fof");
	i = 1;
	for (j=0;j<3;++j) fPeriod[j] = HUGE;
	while (i < argc) {
		if (!strcmp(argv[i],"-e")) {
			++i;
			fEps = atof(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-m")) {
			++i;
			nMembers = atoi(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-o")) {
			++i;
			strcpy(ach,argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-p")) {
			++i;
			fPeriod[0] = atof(argv[i]);
			fPeriod[1] = atof(argv[i]);
			fPeriod[2] = atof(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-px")) {
			++i;
			fPeriod[0] = atof(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-py")) {
			++i;
			fPeriod[1] = atof(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-pz")) {
			++i;
		    fPeriod[2] = atof(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-v")) {
			bVerbose = 1;
			++i;
			}
		else if (*argv[i] == '-') {
			p = argv[i];
			++p;
			if (*p == 'd' || *p == 'g' || *p == 's') {
				bDark = 0;
				bGas = 0;
				bStar = 0;
				}
			else usage();
			while (isalpha(*p)) {
				switch (*p) {
				case 'd':
					bDark = 1;
					break;
				case 'g':
					bGas = 1;
					break;
				case 's':
					bStar = 1;
					break;
				default:
					usage();
					}
				++p;
				}
			++i;
			}
		else usage();
		}
	kdInit(&kd,nBucket,fPeriod);
	kdReadTipsy(kd,stdin,bDark,bGas,bStar);
	kdBuildTree(kd);
	kdTime(kd,&sec,&usec);
	nGroup = kdFoF(kd,fEps);
	kdTime(kd,&sec,&usec);
	if (bVerbose) printf("Number of initial groups:%d\n",nGroup);
	nGroup = kdTooSmall(kd,nMembers);
	if (bVerbose) {
		printf("Number of groups:%d\n",nGroup);
		printf("FOF CPU TIME: %d.%06d secs\n",sec,usec);
		}
	kdOrder(kd);
	strcat(ach,".grp");
	kdOutGroup(kd,ach);
	kdFinish(kd);
	}





