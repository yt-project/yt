// This code was originally written by Bridget Falck and Mark Neyrinck.
// They have agreed to release it under the terms of the BSD license.
#include "origami_tags.h"

int isneg(int h) {
  return (int)(h < 0);
}

int par(int i, int j, int k, int ng) {
  return i + (j + k*ng)*ng;
}

int compute_tags(int ng, double boxsize, double **r, int np,
                 unsigned char *m) {
  /* Note that the particles must be fed in according to the order specified in
   * the README file */
  double predict, xmin,xmax,ymin,ymax,zmin,zmax;
  
  double negb2, b2;
  int ng2,ng4, h, i,i2, x,y,z,nhalo,nhalo0,nhalo1,nhalo2,nhaloany;
  unsigned char *m0,*m1,*m2, mn,m0n,m1n,m2n; /*Morphology tag */
  
  double dx,d1,d2;

  b2 = boxsize/2.;
  negb2 = -boxsize/2.;
  ng2=ng/2;
  ng4=ng/4;

  /* Boxsize should be the range in r, yielding a range 0-1 */
  printf("%d particles\n",np);fflush(stdout);
  xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;
  m0 = (unsigned char *)malloc(np*sizeof(unsigned char)); /* for the diagonals */
  m1 = (unsigned char *)malloc(np*sizeof(unsigned char));
  m2 = (unsigned char *)malloc(np*sizeof(unsigned char));
  for (i=0; i<np;i++) {
    if (r[0][i]<xmin) xmin = r[0][i]; if (r[0][i]>xmax) xmax = r[0][i];
    if (r[1][i]<ymin) ymin = r[1][i]; if (r[1][i]>ymax) ymax = r[1][i];
    if (r[2][i]<zmin) zmin = r[2][i]; if (r[2][i]>zmax) zmax = r[2][i];
   
    m[i] = 1;
    m0[i] = 1;
    m1[i] = 1;
    m2[i] = 1;
  }

  if (m==NULL) {
    printf("Morphology array cannot be allocated.\n");
    return 1;
  }
  //  printf("np: %d, x: %f,%f; y: %f,%f; z: %f,%f\n",np,xmin,xmax, ymin,ymax, zmin,zmax); fflush(stdout);

  printf("Calculating ORIGAMI morphology.\n");

  for (x=0; x<ng; x++){
    //    printf("%d\n",x);fflush(stdout);
    for (y=0; y<ng; y++) {
      for (z=0; z<ng; z++) {
	i = par(x,y,z,ng);
	/* First just along the Cartesian axes */
	/* x-direction */
	  for (h=1; h<ng4; h++) {
	    i2 = par((x+h)%ng,y,z,ng);
	    dx = r[0][i2]-r[0][i];
	    if (dx < negb2) dx += boxsize;
	    if (dx > b2) dx -= boxsize;
	    if (dx < 0.) {
	      /*printf("x:%d %d %d %d %f\n",x,y,z,h,dx);*/
	    if (m[i] % 2 > 0) {
	      m[i] *= 2;
	    }
	      if (m[i2] % 2 > 0){
		m[i2] *= 2;
	      }
	    break;
	    }
	  }
	  for (h=1; h<ng4; h++) {
	    i2 = par(x,(y+h)%ng,z,ng);
	    dx = r[1][i2]-r[1][i];
	    if (dx < negb2) dx += boxsize;
	    if (dx > b2) dx -= boxsize;
	    if (dx < 0.) {
	      /*printf("y:%d %d %d %d %f\n",x,y,z,h,dx);*/
	    if (m[i] % 3 > 0) {
	      m[i] *= 3;
	    }
	      if (m[i2] % 3 > 0){
		m[i2] *= 3;
	      }
	      break;
	    }
	  }
	  for (h=1; h<ng4; h++) {
	    i2 = par(x,y,(z+h)%ng,ng);
	    dx = r[2][i2]-r[2][i];
	    if (dx < negb2) dx += boxsize;
	    if (dx > b2) dx -= boxsize;
	    if (dx < 0.) {
	      /*printf("z:%d %d %d %d %f\n",x,y,z,h,dx);*/
	    if (m[i] % 5 > 0) {
	      m[i] *= 5;
	    }
	      if (m[i2] % 5 > 0){
		m[i2] *= 5;
	      }
	      break;
	    }
	  }
	// Now do diagonal directions 
	for (h=1; h<ng4; h = -h + isneg(h)) {
	  i2 = par(x,goodmod(y+h,ng),goodmod(z+h,ng),ng);
	  d1 = r[1][i2]-r[1][i];
	  d2 = r[2][i2]-r[2][i];
	  if (d1 < negb2) d1 += boxsize;
	  if (d1 > b2) d1 -= boxsize;
	  if (d2 < negb2) d2 += boxsize;
	  if (d2 > b2) d2 -= boxsize;
	  if ((d1 + d2)*h < 0.) {
	    m0[i] *= 2;
	    break;
	  }
	}
	for (h=1; h<ng4; h = -h + isneg(h)) {
	  i2 = par(x,goodmod(y+h,ng),goodmod(z-h,ng),ng);
	  d1 = r[1][i2]-r[1][i];
	  d2 = r[2][i2]-r[2][i];
	  if (d1 < negb2) d1 += boxsize;
	  if (d1 > b2) d1 -= boxsize;
	  if (d2 < negb2) d2 += boxsize;
	  if (d2 > b2) d2 -= boxsize;
	  if ((d1 - d2)*h < 0.) {
	    m0[i] *= 3;
	    break;
	  }
	}
	// y
	for (h=1; h<ng4; h = -h + isneg(h)) {
	  i2 = par(goodmod(x+h,ng),y,goodmod(z+h,ng),ng);
	  d1 = r[0][i2]-r[0][i];
	  d2 = r[2][i2]-r[2][i];
	  if (d1 < negb2) d1 += boxsize;
	  if (d1 > b2) d1 -= boxsize;
	  if (d2 < negb2) d2 += boxsize;
	  if (d2 > b2) d2 -= boxsize;
	  if ((d1 + d2)*h < 0.) {
	    m1[i] *= 2;
	    break;
	  }
	}
	for (h=1; h<ng4; h = -h + isneg(h)) {
	  i2 = par(goodmod(x+h,ng),y,goodmod(z-h,ng),ng);
	  d1 = r[0][i2]-r[0][i];
	  d2 = r[2][i2]-r[2][i];
	  if (d1 < negb2) d1 += boxsize;
	  if (d1 > b2) d1 -= boxsize;
	  if (d2 < negb2) d2 += boxsize;
	  if (d2 > b2) d2 -= boxsize;
	  if ((d1 - d2)*h < 0.) {
	    m1[i] *= 3;
	    break;
	  }
	}
	// z
	for (h=1; h<ng4; h = -h + isneg(h)) {
	  i2 = par(goodmod(x+h,ng),goodmod(y+h,ng),z,ng);
	  d1 = r[0][i2]-r[0][i];
	  d2 = r[1][i2]-r[1][i];
	  if (d1 < negb2) d1 += boxsize;
	  if (d1 > b2) d1 -= boxsize;
	  if (d2 < negb2) d2 += boxsize;
	  if (d2 > b2) d2 -= boxsize;
	  if ((d1 + d2)*h < 0.) {
	    m2[i] *=2;
	    break;
	  }
	}
	for (h=1; h<ng4; h = -h + isneg(h)) {
	  i2 = par(goodmod(x+h,ng),goodmod(y-h,ng),z,ng);
	  d1 = r[0][i2]-r[0][i];
	  d2 = r[1][i2]-r[1][i];
	  if (d1 < negb2) d1 += boxsize;
	  if (d1 > b2) d1 -= boxsize;
	  if (d2 < negb2) d2 += boxsize;
	  if (d2 > b2) d2 -= boxsize;
	  if ((d1 - d2)*h < 0.) {
	    m2[i] *= 3;
	    break;
	  }
	  }
      }
    }
  }

  nhalo = 0;
  nhalo0 = 0;
  nhalo1 = 0;
  nhalo2 = 0;
  nhaloany = 0;
  for (i=0;i<np;i++){
    mn = (m[i]%2 == 0) + (m[i]%3 == 0) + (m[i]%5 == 0);
    m0n = (unsigned char)(m[i]%2 == 0) + (unsigned char)(m0[i]%2 == 0) + (unsigned char)(m0[i]%3 == 0);
    m1n = (unsigned char)(m[i]%3 == 0) + (unsigned char)(m1[i]%2 == 0) + (unsigned char)(m1[i]%3 == 0);
    m2n = (unsigned char)(m[i]%5 == 0) + (unsigned char)(m2[i]%2 == 0) + (unsigned char)(m2[i]%3 == 0);
    m[i] = max(mn,max(m0n,max(m1n,m2n)));
    if (mn == 3) nhalo++;
    if (m0n == 3) nhalo0++;
    if (m1n == 3) nhalo1++;
    if (m2n == 3) nhalo2++;
    if (m[i] == 3) {
      nhaloany++;
    }
  }
  printf("nhalo=%d,%d,%d,%d,%d\n",nhalo,nhalo0,nhalo1,nhalo2,nhaloany);
  free(m0);
  free(m1);
  free(m2);

  /* Output m */
  return 0;
}
