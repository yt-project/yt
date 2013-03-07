#include "artio_endian.h"

#include <stdint.h>

void artio_int_swap(int32_t *src, int count) {
	int i;
	union {
		int32_t f;
		unsigned char c[4];
	} d1, d2;

	for ( i = 0; i < count; i++ ) {
		d1.f = src[i];
		d2.c[0] = d1.c[3];
		d2.c[1] = d1.c[2];
		d2.c[2] = d1.c[1];
		d2.c[3] = d1.c[0];
		src[i] = d2.f;
	}
}

void artio_float_swap(float *src, int count) {
	int i;
	union {
		float f;
		unsigned char c[4];
	} d1, d2;

	for ( i = 0; i < count; i++ ) {
		d1.f = src[i];
		d2.c[0] = d1.c[3];
		d2.c[1] = d1.c[2];
		d2.c[2] = d1.c[1];
		d2.c[3] = d1.c[0];
		src[i] = d2.f;
	}
}

void artio_double_swap(double *src, int count) {	
	int i;
	union
	{
		double d;
		unsigned char c[8];
	} d1, d2;

	for ( i = 0; i < count; i++ ) {
		d1.d = src[i];
		d2.c[0] = d1.c[7];
		d2.c[1] = d1.c[6];
		d2.c[2] = d1.c[5];
		d2.c[3] = d1.c[4];
		d2.c[4] = d1.c[3];
		d2.c[5] = d1.c[2];
		d2.c[6] = d1.c[1];
		d2.c[7] = d1.c[0];
		src[i] = d2.d;
	}
}

void artio_long_swap(int64_t *src, int count) {
	int i;
	union
	{
		int64_t d;
		unsigned char c[8];
	} d1, d2;

	for ( i = 0; i < count; i++ ) {
		d1.d = src[i];
		d2.c[0] = d1.c[7];
		d2.c[1] = d1.c[6];
		d2.c[2] = d1.c[5];
		d2.c[3] = d1.c[4];
		d2.c[4] = d1.c[3];
		d2.c[5] = d1.c[2];
		d2.c[6] = d1.c[1];
		d2.c[7] = d1.c[0];
		src[i] = d2.d;
	}
}
