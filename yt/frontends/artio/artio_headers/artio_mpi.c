/*
 * artio_mpi.c
 *
 *  Created on: Apr 9, 2010
 *      Author: Yongen Yu
 *  Modified: Nov 18, 2010 - Douglas Rudd
 */

#include "artio.h"
#include "artio_internal.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#ifdef ARTIO_MPI

#include "artio_mpi.h"

#define MPI_FH_BUFFERSIZE		16384

struct ARTIO_FH {
	MPI_File fh;
	MPI_Comm comm;
	char *data;
	int mode;
	int bfptr;
	int bfsize;
	int bfend;
};

artio_fh artio_file_fopen( char * filename, int mode, artio_context context) {
	int status;
	int flag;
	int rank;
	int amode;

	if ( mode & ARTIO_MODE_WRITE &&
			mode & ARTIO_MODE_READ ) {
		return NULL;
	} else if ( mode & ARTIO_MODE_WRITE ) {
		amode = MPI_MODE_CREATE | MPI_MODE_WRONLY;
	} else if ( mode & ARTIO_MODE_READ ) {
		amode = MPI_MODE_RDONLY;
	} else {
		return NULL;
	}

	artio_fh ffh = (artio_fh)malloc(sizeof(struct ARTIO_FH));

	ffh->mode = mode;
	flag = mode & ARTIO_MODE_ACCESS;

	MPI_Comm_rank( context->comm, &rank );
	MPI_Comm_split( context->comm, flag, rank, &ffh->comm );

	if ( flag ) {
		status = MPI_File_open( ffh->comm, filename, amode, MPI_INFO_NULL, &ffh->fh);
		if (status != MPI_SUCCESS) {
			return NULL;
		}

		/* truncate the file on write */
		if ( mode & ARTIO_MODE_WRITE ) {
			MPI_File_set_size(ffh->fh, 0);
		}

		if ( mode & ARTIO_MODE_DIRECT) {
			ffh->bfsize = 0;
			ffh->bfend = 0;
			ffh->bfptr = 0;
		} else {
		  ffh->data = (char *)malloc(MPI_FH_BUFFERSIZE);
			memset(ffh->data, 0, MPI_FH_BUFFERSIZE);
			ffh->bfptr = 0;
			ffh->bfsize = MPI_FH_BUFFERSIZE;
			ffh->bfend = -1;
		}
	}

	return ffh;
}

int artio_file_fwrite(artio_fh handle, void *buf, int count, int type ) {
	int size, part, remain;

	if ( !(handle->mode & ARTIO_MODE_WRITE) ||
			!(handle->mode & ARTIO_MODE_ACCESS) ) {
		return ARTIO_ERR_INVALID_FILE_MODE;
	}

	size = count*artio_type_size( type );

	if ( handle->mode & ARTIO_MODE_DIRECT ) {
		MPI_File_write( handle->fh, buf, size, 
				MPI_BYTE, MPI_STATUS_IGNORE );
	} else {
		if ( handle->bfptr + size < handle->bfsize) {
			memcpy( handle->data + handle->bfptr, buf, size );
			handle->bfptr += size;
		} else {
			part = handle->bfsize - handle->bfptr;
			remain = handle->bfptr + size - handle->bfsize;

			memcpy( handle->data + handle->bfptr, buf, part);
			MPI_File_write(handle->fh, handle->data, handle->bfsize, 
					MPI_BYTE, MPI_STATUS_IGNORE );

			while ( remain > handle->bfsize ) {
				memcpy( handle->data, (char *)buf + part, handle->bfsize);
				MPI_File_write(handle->fh, handle->data, handle->bfsize, 
						MPI_BYTE, MPI_STATUS_IGNORE );
				remain -= handle->bfsize;
				part += handle->bfsize;
			}

			memcpy( handle->data, (char *)buf + part, remain);
			handle->bfptr = remain;
		}
	}

	return ARTIO_SUCCESS;
}

int artio_file_fflush(artio_fh handle) {
	if ( !(handle->mode & ARTIO_MODE_ACCESS) ) {
		return ARTIO_ERR_INVALID_FILE_MODE;
	}

	if ( handle->mode & ARTIO_MODE_WRITE ) {
		if ( handle->bfptr > 0 ) {
			MPI_File_write(handle->fh, handle->data, handle->bfptr, 
					MPI_BYTE, MPI_STATUS_IGNORE );
			handle->bfptr = 0;
		}
	} else if ( handle->mode & ARTIO_MODE_READ ) {
		handle->bfptr = 0;
		handle->bfend = -1;
	} else {
		return ARTIO_ERR_INVALID_FILE_MODE;
	}

	return ARTIO_SUCCESS;
}

int artio_file_fread(artio_fh handle, void *buf, int count, int type ) {
	MPI_Status status;

	int size, avail, remain, size_read;
	char *curbuf;

	if ( !(handle->mode & ARTIO_MODE_READ) ) {
		return ARTIO_ERR_INVALID_FILE_MODE;
	}

	size = count*artio_type_size( type );

	if ( handle->mode & ARTIO_MODE_DIRECT ) {
		MPI_File_read(handle->fh, buf, size, MPI_BYTE, &status );
		MPI_Get_count( &status, MPI_BYTE, &size_read );
		if ( size_read != size ) {
			return ARTIO_ERR_INSUFFICIENT_DATA;
		}
	} else {
		if ( handle->bfend == -1 ) {
			/* load initial data into buffer */
			MPI_File_read(handle->fh, handle->data, handle->bfsize, MPI_BYTE, &status);
			MPI_Get_count(&status, MPI_BYTE, &handle->bfend);
			handle->bfptr = 0;
		}

		/* read from buffer */
		curbuf = (char *)buf;
		remain = size;
		while ( remain > 0 && 
				handle->bfend > 0 && 
				handle->bfptr + remain >= handle->bfend ) {
			avail = handle->bfend - handle->bfptr;
			memcpy( curbuf, (char *)handle->data + handle->bfptr, avail );
			curbuf += avail;
			remain -= avail;

			/* refill buffer */
			MPI_File_read(handle->fh, handle->data, handle->bfsize, 
					MPI_BYTE, &status );
			MPI_Get_count(&status, MPI_BYTE, &handle->bfend );
			handle->bfptr = 0;
		}

		if ( handle->bfend == 0 ) {
			/* ran out of data, eof */
			return ARTIO_ERR_INSUFFICIENT_DATA;
		}

		memcpy( curbuf, (char *)handle->data + handle->bfptr, remain );
		handle->bfptr += remain;
	}

	if(handle->mode & ARTIO_MODE_ENDIAN_SWAP){
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

int artio_file_fseek(artio_fh handle, int64_t offset, int whence ) {
	MPI_Offset current;

	if (  handle->mode & ARTIO_MODE_ACCESS ) {
		if ( whence == ARTIO_SEEK_CUR ) {
			if ( offset == 0 ) {
				return ARTIO_SUCCESS;
			} else if ( handle->mode & ARTIO_MODE_READ &&
                	handle->bfptr + offset >= 0 && 
					handle->bfptr + offset < handle->bfsize ) {
				handle->bfptr += offset;
				return ARTIO_SUCCESS;
			} else { 
				artio_file_fflush( handle );
				current = (MPI_Offset)offset;
				MPI_File_seek( handle->fh, current, MPI_SEEK_CUR );
			}
		} else if ( whence == ARTIO_SEEK_SET ) {
			MPI_File_get_position( handle->fh, &current );
			if (handle->mode & ARTIO_MODE_WRITE &&
					current<=offset && offset<(current + handle->bfsize) && 
					handle->bfptr==(offset - current)) {
				return ARTIO_SUCCESS;
			} else if ( handle->mode & ARTIO_MODE_READ && 
					handle->bfend > 0 &&
					offset >= (current - handle->bfend) &&
					offset < current ) {
				handle->bfptr = offset - current + handle->bfend;
			} else {
				artio_file_fflush( handle );
				current = (MPI_Offset)offset;
				MPI_File_seek( handle->fh, current, MPI_SEEK_SET );
			}
		} else {
			/* unknown whence */
			return ARTIO_ERR_INVALID_SEEK;
		}
	} else {
		/* seek on non-active file handle */
		return ARTIO_ERR_INVALID_FILE_MODE;
	}

	return ARTIO_SUCCESS;
}

int artio_file_fclose(artio_fh handle) {
	if ( handle->mode & ARTIO_MODE_ACCESS ) {
		artio_file_fflush(handle);
		MPI_File_close(&handle->fh);

		if ( handle->bfsize > 0 ) {
			free(handle->data);
		}
	}
	MPI_Comm_free(&handle->comm);
	return ARTIO_SUCCESS;
}

void artio_set_endian_swap_tag(artio_fh handle) {
	handle->mode |= ARTIO_MODE_ENDIAN_SWAP;
}

#endif /* MPI */
