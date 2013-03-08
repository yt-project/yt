#ifndef __ARTIO_EDIAN_H__
#define __ARTIO_EDIAN_H__

#include <stdint.h>

void artio_int_swap(int32_t *src, int count);
void artio_float_swap(float *src, int count);
void artio_double_swap(double *src, int count);
void artio_long_swap(int64_t *src, int count);

#endif /* __ARTIO_ENDIAN_H__ */
