#include <vector>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <stdint.h>
#include <float.h>
#include <cstdlib>

bool isEqual(double f1, double f2);
double* max_pts(double *pts, uint64_t n, uint32_t m);
double* min_pts(double *pts, uint64_t n, uint32_t m);
uint64_t argmin_pts_dim(double *pts, uint64_t *idx,
                        uint32_t m, uint32_t d,
                        uint64_t Lidx, uint64_t Ridx);
uint64_t argmax_pts_dim(double *pts, uint64_t *idx,
                        uint32_t m, uint32_t d,
                        uint64_t Lidx, uint64_t Ridx);
void quickSort(double *pts, uint64_t *idx,
               uint32_t ndim, uint32_t d,
               int64_t l, int64_t r);
void insertSort(double *pts, uint64_t *idx,
                uint32_t ndim, uint32_t d,
                int64_t l, int64_t r);
int64_t pivot(double *pts, uint64_t *idx,
              uint32_t ndim, uint32_t d,
              int64_t l, int64_t r);
int64_t partition_given_pivot(double *pts, uint64_t *idx,
                              uint32_t ndim, uint32_t d,
                              int64_t l, int64_t r, double pivot);
int64_t partition(double *pts, uint64_t *idx,
                  uint32_t ndim, uint32_t d,
                  int64_t l, int64_t r, int64_t p);
int64_t select(double *pts, uint64_t *idx,
               uint32_t ndim, uint32_t d,
               int64_t l, int64_t r, int64_t n);
uint32_t split(double *all_pts, uint64_t *all_idx,
               uint64_t Lidx, uint64_t n, uint32_t ndim,
               double *mins, double *maxes,
               int64_t &split_idx, double &split_val,
	       bool use_sliding_midpoint = false);
