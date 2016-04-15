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

#include "artio.h"
#include "artio_internal.h"

#ifndef ARTIO_MPI 

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#ifdef _WIN32
typedef __int64 int64_t;
typedef __int32 int32_t;
#else
#include <stdint.h>
#endif

struct ARTIO_FH {
	FILE *fh;
	int mode;
	char *data;                                                                                                  
	int bfptr;
	int bfsize;
	int bfend;
};

#ifdef _WIN32
#define FOPEN_FLAGS "b"
#define fseek _fseeki64
#else
#define FOPEN_FLAGS ""
#endif

artio_context artio_context_global_struct = { 0 };
const artio_context *artio_context_global = &artio_context_global_struct;

artio_fh *artio_file_fopen_i( char * filename, int mode, const artio_context *not_used ) {
	artio_fh *ffh;
	/* check for invalid combination of mode parameter */
	if ( ( mode & ARTIO_MODE_READ && mode & ARTIO_MODE_WRITE ) ||
			!( mode & ARTIO_MODE_READ || mode & ARTIO_MODE_WRITE ) ) {
		return NULL;
	}

	ffh = (artio_fh *)malloc(sizeof(artio_fh));
	if ( ffh == NULL ) {
		return NULL;
	}
	
	ffh->mode = mode;
	ffh->bfsize = -1;
	ffh->bfend = -1;
	ffh->bfptr = -1;
	ffh->data = NULL;

	if ( mode & ARTIO_MODE_ACCESS ) {
		ffh->fh = fopen( filename, ( mode & ARTIO_MODE_WRITE ) ? "w"FOPEN_FLAGS : "r"FOPEN_FLAGS );
		if ( ffh->fh == NULL ) {
			free( ffh );
			return NULL;
		}
	} 

	return ffh;
}

int artio_file_attach_buffer_i( artio_fh *handle, void *buf, int buf_size ) {
	if ( !(handle->mode & ARTIO_MODE_ACCESS ) ) {
		return ARTIO_ERR_INVALID_FILE_MODE;
	}

	if ( handle->data != NULL ) {
		return ARTIO_ERR_BUFFER_EXISTS;
	}
	
	handle->bfsize = buf_size;
	handle->bfend = -1;
	handle->bfptr = 0;
	handle->data = (char *)buf;

	return ARTIO_SUCCESS;
}

int artio_file_detach_buffer_i( artio_fh *handle ) {
	int ret;
	ret = artio_file_fflush(handle);
	if ( ret != ARTIO_SUCCESS ) return ret;

	handle->data = NULL;
	handle->bfsize = -1;
    handle->bfend = -1;
    handle->bfptr = -1;

	return ARTIO_SUCCESS;	
}

int artio_file_fwrite_i(artio_fh *handle, const void *buf, int64_t count, int type ) {
	size_t size;
	int64_t remain;
	char *p;
	int size32;

	if ( !(handle->mode & ARTIO_MODE_WRITE) ||
			!(handle->mode & ARTIO_MODE_ACCESS) ) {
		return ARTIO_ERR_INVALID_FILE_MODE;
	}

	size = artio_type_size( type );
	if ( size == (size_t)-1 ) {
		return ARTIO_ERR_INVALID_DATATYPE;
	}

	if ( count > ARTIO_INT64_MAX / size ) {
		return ARTIO_ERR_IO_OVERFLOW;
	}

	remain = count*size;
	p = (char *)buf;

	if ( handle->data == NULL ) {
		/* force writes to 32-bit sizes */
		while ( remain > 0 ) {
			size32 = MIN( ARTIO_IO_MAX, remain );		
			if ( fwrite( p, 1, size32, handle->fh ) != size32 ) {
				return ARTIO_ERR_IO_WRITE;
			}
			remain -= size32;
			p += size32;
		}
	} else if ( remain < handle->bfsize - handle->bfptr ) {
		memcpy( handle->data + handle->bfptr, p, (size_t)remain );
		handle->bfptr += remain;
	} else {
		size32 = handle->bfsize - handle->bfptr;
		memcpy( handle->data + handle->bfptr, p, size32 );
		if ( fwrite( handle->data, 1, handle->bfsize, 
				handle->fh ) != handle->bfsize ) {
			return ARTIO_ERR_IO_WRITE;
		}
		p += size32;
		remain -= size32;

		while ( remain > handle->bfsize ) {
			/* write directly to file-handle in unbuffered case */
			if ( fwrite( p, 1, handle->bfsize, 
					handle->fh ) != handle->bfsize ) {
				return ARTIO_ERR_IO_WRITE;
			}
			remain -= handle->bfsize;
			p += handle->bfsize;
		}

		memcpy( handle->data, p, (size_t)remain);
		handle->bfptr = remain;
	}

	return ARTIO_SUCCESS;
}

int artio_file_fflush_i(artio_fh *handle) {
	if ( !(handle->mode & ARTIO_MODE_ACCESS) ) {
		return ARTIO_ERR_INVALID_FILE_MODE;
	}

    if ( handle->mode & ARTIO_MODE_WRITE ) {
		if ( handle->bfptr > 0 ) {
			if ( fwrite( handle->data, 1, handle->bfptr, 
					handle->fh ) != handle->bfptr ) {
				return ARTIO_ERR_IO_WRITE;
			}
			handle->bfptr = 0;
		}
	} else if ( handle->mode & ARTIO_MODE_READ ) {
		handle->bfend = -1;
		handle->bfptr = 0;
	} else {
        return ARTIO_ERR_INVALID_FILE_MODE;
    }

	return ARTIO_SUCCESS;
}

int artio_file_fread_i(artio_fh *handle, void *buf, int64_t count, int type ) {
	size_t size, avail, remain;
	int size32;
	char *p;

	if ( !(handle->mode & ARTIO_MODE_READ) ) {
		return ARTIO_ERR_INVALID_FILE_MODE;
	}

	size = artio_type_size( type );
	if ( size == (size_t)-1 ) {
		return ARTIO_ERR_INVALID_DATATYPE;
	}

	if ( count > ARTIO_INT64_MAX / size ) {
		return ARTIO_ERR_IO_OVERFLOW;
	}

	remain = size*count;
	p = (char *)buf;

	if ( handle->data == NULL ) {
		while ( remain > 0 ) {
			size32 = MIN( ARTIO_IO_MAX, remain );
			if ( fread( p, 1, size32, handle->fh) != size32 ) {
				return ARTIO_ERR_INSUFFICIENT_DATA;
			}
			remain -= size32;
			p += size32;
		}
	} else {
		if ( handle->bfend == -1 ) {
			/* load initial data into buffer */
			handle->bfend = fread( handle->data, 1, handle->bfsize, handle->fh );
			handle->bfptr = 0;
		}

		/* read from buffer */
		while ( remain > 0 && 
				handle->bfend > 0 && 
				handle->bfptr + remain >= handle->bfend ) {
			avail = handle->bfend - handle->bfptr;
			memcpy( p, handle->data + handle->bfptr, avail );
			p += avail;
			remain -= avail;

			/* refill buffer */
			handle->bfend = fread( handle->data, 1, handle->bfsize, handle->fh );
            handle->bfptr = 0;
		}

		if ( remain > 0 ) {
			if ( handle->bfend == 0 ) {
				/* ran out of data, eof */
				return ARTIO_ERR_INSUFFICIENT_DATA;
			}

			memcpy( p, handle->data + handle->bfptr, (size_t)remain );
			handle->bfptr += (int)remain;
		}
	}

	if(handle->mode & ARTIO_MODE_ENDIAN_SWAP) {
		switch (type) {
			case ARTIO_TYPE_INT :
				artio_int_swap( (int32_t *)buf, count );
				break;
			case ARTIO_TYPE_FLOAT :
				artio_float_swap( (float *)buf, count );
				break;
			case ARTIO_TYPE_DOUBLE :
				artio_double_swap( (double *)buf, count );
				break;
			case ARTIO_TYPE_LONG :
				artio_long_swap( (int64_t *)buf, count );
				break;
			default :
				return ARTIO_ERR_INVALID_DATATYPE;
		}
	}

    return ARTIO_SUCCESS;
}

int artio_file_ftell_i( artio_fh *handle, int64_t *offset ) {
	size_t current = ftell( handle->fh );

	if ( handle->bfend > 0 ) {
		current -= handle->bfend;
	} 
	if ( handle->bfptr > 0 ) {
		current += handle->bfptr;
	}
	*offset = (int64_t)current;

    return ARTIO_SUCCESS;
}

int artio_file_fseek_i(artio_fh *handle, int64_t offset, int whence ) {
	size_t current;

	if ( handle->mode & ARTIO_MODE_ACCESS ) {
		if ( whence == ARTIO_SEEK_CUR ) {
			if ( offset == 0 ) {
				return ARTIO_SUCCESS;
			} else if ( handle->mode & ARTIO_MODE_READ &&
					handle->bfend > 0 &&
					handle->bfptr + offset >= 0 && 
					handle->bfptr + offset < handle->bfend ) {
				handle->bfptr += offset;
				return ARTIO_SUCCESS;
            } else {
				/* modify offset due to offset in buffer */
				if ( handle->bfptr > 0 ) {
					current = offset - handle->bfend + handle->bfptr;
				} else {
					current = offset;
				}
                artio_file_fflush( handle );
				fseek( handle->fh, (size_t)current, SEEK_CUR );
			}
		} else if ( whence == ARTIO_SEEK_SET ) {
			current = ftell( handle->fh );
			if ( handle->mode & ARTIO_MODE_WRITE &&
					current <= offset && 
					offset < current + handle->bfsize &&                                      
					handle->bfptr == offset - current ) {
				return ARTIO_SUCCESS;
			} else if ( handle->mode & ARTIO_MODE_READ &&
					handle->bfptr > 0 &&
					handle->bfend > 0 &&
					handle->bfptr < handle->bfend &&
					offset >= current - handle->bfend &&
					offset < current ) {
				handle->bfptr = offset - current + handle->bfend;
			} else {
				artio_file_fflush( handle );
				fseek( handle->fh, (size_t)offset, SEEK_SET );
			}
		} else if ( whence == ARTIO_SEEK_END ) {
			artio_file_fflush( handle );
			fseek( handle->fh, (size_t)offset, SEEK_END );
		} else {
			/* unknown whence */
			return ARTIO_ERR_INVALID_SEEK;
		}
	} else {
		return ARTIO_ERR_INVALID_FILE_MODE;
	}

	return ARTIO_SUCCESS;
}

int artio_file_fclose_i(artio_fh *handle) {
	if ( handle->mode & ARTIO_MODE_ACCESS ) {
		artio_file_fflush(handle);
		fclose(handle->fh);
	}
	free(handle);

	return ARTIO_SUCCESS;
}

void artio_file_set_endian_swap_tag_i(artio_fh *handle) {
	handle->mode |= ARTIO_MODE_ENDIAN_SWAP;
}

#endif /* ifndef ARTIO_MPI */
