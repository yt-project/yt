/**********************************************************************
 * Copyright (c) 2012-2013, Douglas H. Rudd
 * All rights reserved.
 *
 * This file is part of the artio library.
 *
 * artio is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * artio is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * Copies of the GNU Lesser General Public License and the GNU General
 * Public License are available in the file LICENSE, included with this
 * distribution.  If you failed to receive a copy of this file, see
 * <http://www.gnu.org/licenses/>
 **********************************************************************/

#include "artio_endian.h"

#ifdef _WIN32
typedef __int64 int64_t;
typedef __int32 int32_t;
#else
#include <stdint.h>
#endif

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
