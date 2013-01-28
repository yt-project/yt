/*
 * artio_posix.c
 *
 *  Created on: Apr 9, 2010
 *      Author: Yongen Yu
 *  Modified: Nov 18, 2010 - Douglas Rudd
 */

#include "artio.h"
#include "artio_internal.h"

#ifndef ARTIO_MPI 

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

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

#define FH_BUFFERSIZE		16384

artio_context artio_context_global_struct = { 0 };
const artio_context *artio_context_global = &artio_context_global_struct;

artio_fh *artio_file_fopen( char * filename, int mode, const artio_context *not_used ) {
	int flag;

	artio_fh *ffh = (artio_fh *)malloc(sizeof(artio_fh));
	if ( ffh == NULL ) {
		return NULL;
	}
	
	ffh->mode = mode;
	flag = mode & ARTIO_MODE_ACCESS;

	if ( flag ) {
		if ( mode & ARTIO_MODE_DIRECT) {
			ffh->data = NULL;
			ffh->bfsize = 0;
			ffh->bfend = 0;
			ffh->bfptr = 0;
		} else {
			ffh->data = (char *)malloc(FH_BUFFERSIZE);
			if ( ffh->data == NULL ) {
				free(ffh);
				return NULL;
			}
			memset(ffh->data, 0, FH_BUFFERSIZE);
			ffh->bfptr = 0;
			ffh->bfsize = FH_BUFFERSIZE;
			ffh->bfend = -1;
		}

		ffh->fh = fopen( filename, ( mode & ARTIO_MODE_WRITE ) ? "w"FOPEN_FLAGS : "r"FOPEN_FLAGS );
		if ( ffh->fh == NULL ) {
			free( ffh->data );
			free( ffh );
			return NULL;
		}
	} 

	return ffh;
}

int artio_file_fwrite(artio_fh *handle, void *buf, int64_t count, int type ) {
	size_t size, part, remain;

	if ( !(handle->mode & ARTIO_MODE_WRITE) ||
			!(handle->mode & ARTIO_MODE_ACCESS) ) {
		return ARTIO_ERR_INVALID_FILE_MODE;
	}

	size = artio_type_size( type );
	if ( size < 0 ) {
		return ARTIO_ERR_INVALID_DATATYPE;
	}
	size *= count;

	if ( handle->mode & ARTIO_MODE_DIRECT ) {	
		fwrite( buf, size, (size_t)count, handle->fh );
	} else {
		if ( handle->bfptr + size < handle->bfsize) {
			memcpy( handle->data + handle->bfptr, buf, size );
			handle->bfptr += size;
		} else {
			part = handle->bfsize - handle->bfptr;
			remain = handle->bfptr + size - handle->bfsize;

			memcpy( handle->data + handle->bfptr, buf, part);
			fwrite( handle->data, 1, handle->bfsize, handle->fh );

			while ( remain > handle->bfsize ) {
				fwrite( (char *)buf + part, 1, handle->bfsize, handle->fh );
				/* memcpy( handle->data, (char *)buf + part, handle->bfsize);
				fwrite( handle->data, handle->bfsize, 1, handle->fh );
				*/
				remain -= handle->bfsize;
				part += handle->bfsize;
			}

			memcpy( handle->data, (char *)buf + part, remain);
			handle->bfptr = remain;
		}
	}

	return ARTIO_SUCCESS;
}

int artio_file_fflush(artio_fh *handle) {
	if ( !(handle->mode & ARTIO_MODE_ACCESS) ) {
		return ARTIO_ERR_INVALID_FILE_MODE;
	}

    if ( handle->mode & ARTIO_MODE_WRITE ) {
		if ( handle->bfptr > 0 ) {
			fwrite( handle->data, 1, handle->bfptr, handle->fh );
			handle->bfptr = 0;
		}
	} else if ( handle->mode & ARTIO_MODE_READ ) {
		handle->bfptr = 0;
		handle->bfend = -1;
	} else {
        return ARTIO_ERR_INVALID_FILE_MODE;
    }

	/* match semantics of MPI call */
	/* fflush(handle->fh); */
	return ARTIO_SUCCESS;
}

int artio_file_fread(artio_fh *handle, void *buf, int64_t count, int type ) {
	size_t size, avail, remain, size_read;
	char *curbuf;

	if ( !(handle->mode & ARTIO_MODE_READ) ) {
		return ARTIO_ERR_INVALID_FILE_MODE;
	}

	size = artio_type_size( type );
	if ( size < 0 ) {
		return ARTIO_ERR_INVALID_DATATYPE;
	}

	if ( handle->mode & ARTIO_MODE_DIRECT ) {
		size_read = fread( buf, size, (size_t)count, handle->fh);
		if ( size_read != count ) {
			return ARTIO_ERR_INSUFFICIENT_DATA;
		}
	} else {
		if ( handle->bfend == -1 ) {
			/* load initial data into buffer */
			size_read = fread( handle->data, 1, handle->bfsize, handle->fh );
			handle->bfend = size_read;
			handle->bfptr = 0;
		}

		/* read from buffer */
		curbuf = (char *)buf;
		remain = size*count;
		while ( remain > 0 && 
				handle->bfend > 0 && 
				handle->bfptr + remain >= handle->bfend ) {

			avail = handle->bfend - handle->bfptr;
			memcpy( curbuf, (char *)handle->data + handle->bfptr, avail );
			curbuf += avail;
			remain -= avail;

			/* refill buffer */
			size_read = fread( handle->data, 1, handle->bfsize, handle->fh );
            handle->bfend = size_read;
            handle->bfptr = 0;
		}

		if ( remain > 0 ) {
			if ( handle->bfend == 0 ) {
				/* ran out of data, eof */
				return ARTIO_ERR_INSUFFICIENT_DATA;
			}

			memcpy( curbuf, (char *)handle->data + handle->bfptr, remain );
			handle->bfptr += remain;
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

int artio_file_ftell( artio_fh *handle, int64_t *offset ) {
	size_t current = ftell( handle->fh );

	if ( handle->bfend < 0 ) {
		*offset = current;
	} else {
        *offset = current - handle->bfend + handle->bfptr;
    }

    return ARTIO_SUCCESS;
}

int artio_file_fseek(artio_fh *handle, int64_t offset, int whence ) {
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
					handle->bfptr < handle->bfend &&
					handle->bfend > 0 &&
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

int artio_file_fclose(artio_fh *handle) {
	if ( handle->mode & ARTIO_MODE_ACCESS ) {
		artio_file_fflush(handle);
		fclose(handle->fh);

		if ( handle->data != NULL ) {
			free( handle->data );
		}
	}
	free(handle);

	return ARTIO_SUCCESS;
}

void artio_set_endian_swap_tag(artio_fh *handle) {
	handle->mode |= ARTIO_MODE_ENDIAN_SWAP;
}

#endif /* ifndef ARTIO_MPI */
