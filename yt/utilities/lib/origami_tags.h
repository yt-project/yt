#ifndef __ORIGAMI_TAGS_H__
#define __ORIGAMI_TAGS_H__
#include <stdio.h>
#include <stdlib.h>

#define BF 1e30
#define max(A,B) (((A)>(B)) ? (A):(B))
#define goodmod(A,B) (((A) >= (B)) ? (A-B):(((A) < 0) ? (A+B):(A)))

int isneg(int h);
int par(int i, int j, int k, int ng);
int compute_tags(int ng, double boxsize, double **r, int npart,
                 unsigned char *m);
#endif // __ORIGAMI_TAGS_H__
