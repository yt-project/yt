/*
 * artio_posix.c
 *
 *  Created on: Apr 9, 2010
 *      Author: Yongen Yu
 *  Modified: Nov 18, 2010 - Douglas Rudd
 */

#include "artio.h"
#include "artio_internal.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#ifndef ARTIO_MPI

struct ARTIO_FH {
	FILE *fh;
	int mode;
};


#ifdef _WIN32
#define FOPEN_FLAGS "b"
#define fseek _fseeki64
#else
#define FOPEN_FLAGS ""
#endif

struct artio_context_struct artio_context_global_struct = { 0 };
artio_context artio_context_global = &artio_context_global_struct;

artio_fh artio_file_fopen( char * filename, int mode, artio_context not_used ) {
	int flag;

	artio_fh ffh = (artio_fh)malloc(sizeof(struct ARTIO_FH));
	
	ffh->mode = mode;
	flag = mode & ARTIO_MODE_ACCESS;

	if ( flag ) {
		ffh->fh = fopen( filename, ( mode & ARTIO_MODE_WRITE ) ? "w"FOPEN_FLAGS : "r"FOPEN_FLAGS );
		if ( ffh->fh == NULL ) {
			free( ffh );
			return NULL;
		}
	} 

	return ffh;
}

int artio_file_fwrite(artio_fh handle, void *buf, int64_t count, int type ) {
	int size;

	if ( !(handle->mode & ARTIO_MODE_WRITE) ||
			!(handle->mode & ARTIO_MODE_ACCESS) ) {
		return ARTIO_ERR_INVALID_FILE_MODE;
	}

	size = artio_type_size( type );
	fwrite( buf, (size_t)count, size, handle->fh );
	return ARTIO_SUCCESS;
}

int artio_file_fflush(artio_fh handle) {
	if ( !(handle->mode & ARTIO_MODE_ACCESS) ) {
		return ARTIO_ERR_INVALID_FILE_MODE;
	}

	fflush(handle->fh);
	return ARTIO_SUCCESS;
}

int artio_file_fread(artio_fh handle, void *buf, int64_t count, int type ) {
	int size;

	if ( !(handle->mode & ARTIO_MODE_READ) ) {
		return ARTIO_ERR_INVALID_FILE_MODE;
	}

	size = artio_type_size( type );
	if ( fread( buf, size, (size_t)count, handle->fh ) == count ) {
		return ARTIO_SUCCESS;
	} else {
		return ARTIO_ERR_INSUFFICIENT_DATA;
	}
}

int artio_file_fseek(artio_fh handle, int64_t offset, int whence ) {
	int64_t current;

	if (  handle->mode & ARTIO_MODE_ACCESS ) {
		if ( whence == ARTIO_SEEK_CUR ) {
			if ( offset == 0 ) {
				return ARTIO_SUCCESS;
			} else { 
				fseek( handle->fh, offset, SEEK_CUR );
			}
		} else if ( whence == ARTIO_SEEK_SET ) {
			current = ftell( handle->fh );

			if ( handle->mode & ARTIO_MODE_WRITE && 
					current == offset ) {
				return ARTIO_SUCCESS;
			} else {
				artio_file_fflush( handle );
				fseek( handle->fh, offset, SEEK_SET );
			}
		} else {
			/* unknown whence */
			return ARTIO_ERR_INVALID_SEEK;
		}
	} else {
		return ARTIO_ERR_INVALID_FILE_MODE;
	}

	return ARTIO_SUCCESS;
}

int artio_file_fclose(artio_fh handle) {
	if ( handle->mode & ARTIO_MODE_ACCESS ) {
		artio_file_fflush(handle);
		fclose(handle->fh);
	}
	return ARTIO_SUCCESS;
}

void artio_set_endian_swap_tag(artio_fh handle) {
	handle->mode |= ARTIO_MODE_ENDIAN_SWAP;
}

#endif /* ifndef ARTIO_MPI */
