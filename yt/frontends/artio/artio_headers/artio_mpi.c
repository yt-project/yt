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

#ifdef ARTIO_MPI

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef _WIN32
typedef __int64 int64_t;
typedef __int32 int32_t;
#else
#include <stdint.h>
#endif

artio_context artio_context_global_struct = { MPI_COMM_WORLD };
const artio_context *artio_context_global = &artio_context_global_struct;

struct ARTIO_FH {
	MPI_File fh;
	MPI_Comm comm;
	int mode;
	char *data;
	int bfsize;
	int bfptr;
	int bfend;
};

artio_fh *artio_file_fopen_i( char * filename, int mode, const artio_context *context) {
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

	artio_fh *ffh = (artio_fh *)malloc(sizeof(artio_fh));
	if ( ffh == NULL ) {
		return NULL;
	}

	ffh->mode = mode;
	ffh->data = NULL;
	ffh->bfsize = -1;
	ffh->bfend = -1;
	ffh->bfptr = -1;

	flag = mode & ARTIO_MODE_ACCESS;

	MPI_Comm_rank( context->comm, &rank );
	MPI_Comm_split( context->comm, flag, rank, &ffh->comm );

	if ( flag ) {
		status = MPI_File_open( ffh->comm, filename, amode, MPI_INFO_NULL, &ffh->fh);
		if (status != MPI_SUCCESS) {
			MPI_Comm_free(&ffh->comm);
			free( ffh );
			return NULL;
		}

		/* truncate the file on write */
		if ( mode & ARTIO_MODE_WRITE ) {
			MPI_File_set_size(ffh->fh, 0);
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

int artio_file_fwrite_i( artio_fh *handle, const void *buf, int64_t count, int type ) {
	size_t size;
	int64_t remain;
	int size32;
	char *p;

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

	remain = size*count;
	p = (char *)buf;

	if ( handle->data == NULL ) {
		while ( remain > 0 ) {
			size32 = MIN( ARTIO_IO_MAX, remain );
			if ( MPI_File_write( handle->fh, p, size32, 
					MPI_BYTE, MPI_STATUS_IGNORE ) != MPI_SUCCESS ) {
				return ARTIO_ERR_IO_WRITE;
			}
			remain -= size32;
			p += size32;
		}
	} else if ( remain < handle->bfsize - handle->bfptr ) {
		memcpy( handle->data + handle->bfptr, p, (size_t)remain );
		handle->bfptr += remain;
	} else {
		/* complete buffer */
		size32 = handle->bfsize - handle->bfptr;
		memcpy( handle->data + handle->bfptr, p, size32 );
		if ( MPI_File_write(handle->fh, handle->data, handle->bfsize, 
					MPI_BYTE, MPI_STATUS_IGNORE ) != MPI_SUCCESS ) {
			return ARTIO_ERR_IO_WRITE;
		}
		p += size32;
		remain -= size32;

		while ( remain > handle->bfsize ) {
			if ( MPI_File_write(handle->fh, p, handle->bfsize, 
						MPI_BYTE, MPI_STATUS_IGNORE ) != MPI_SUCCESS ) {
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
			if ( MPI_File_write(handle->fh, handle->data, handle->bfptr, 
					MPI_BYTE, MPI_STATUS_IGNORE ) != MPI_SUCCESS ) {
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
	MPI_Status status;
	size_t size, avail, remain;
	int size_read, size32;
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
			if ( MPI_File_read(handle->fh, p, size32, 
					MPI_BYTE, &status ) != MPI_SUCCESS ) {
				return ARTIO_ERR_IO_READ;
			}
			MPI_Get_count( &status, MPI_BYTE, &size_read );
			if ( size_read != size32 ) {
				return ARTIO_ERR_INSUFFICIENT_DATA;
			}
			remain -= size32;
			p += size32;
		}
	} else {
		if ( handle->bfend == -1 ) {
			/* load initial data into buffer */
			if ( MPI_File_read(handle->fh, handle->data, 
					handle->bfsize, MPI_BYTE, &status) != MPI_SUCCESS ) {
				return ARTIO_ERR_IO_READ;
			}
			MPI_Get_count(&status, MPI_BYTE, &handle->bfend);
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
			if ( MPI_File_read(handle->fh, handle->data, handle->bfsize, 
					MPI_BYTE, &status ) != MPI_SUCCESS ) {
				return ARTIO_ERR_IO_READ;
			}
			MPI_Get_count(&status, MPI_BYTE, &handle->bfend );
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

int artio_file_ftell_i(artio_fh *handle, int64_t *offset) {
	MPI_Offset current;
	MPI_File_get_position( handle->fh, &current );

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
	MPI_Offset current;

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
				if ( handle->bfptr > 0 ) {
					current = (MPI_Offset)offset - handle->bfend + handle->bfptr;
				} else {
					current = (MPI_Offset)offset;
				}

				artio_file_fflush( handle );
				MPI_File_seek( handle->fh, current, MPI_SEEK_CUR );
			}
		} else if ( whence == ARTIO_SEEK_SET ) {
			MPI_File_get_position( handle->fh, &current );
			if (handle->mode & ARTIO_MODE_WRITE &&
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
				MPI_File_seek( handle->fh, (MPI_Offset)offset, MPI_SEEK_SET );
			}
		} else if ( whence == ARTIO_SEEK_END ) {
			artio_file_fflush(handle);
			MPI_File_seek( handle->fh, (MPI_Offset)offset, MPI_SEEK_END );	
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

int artio_file_fclose_i(artio_fh *handle) {
	if ( handle->mode & ARTIO_MODE_ACCESS ) {
		artio_file_fflush(handle);
		MPI_File_close(&handle->fh);
	}
	MPI_Comm_free(&handle->comm);
	free(handle);
	return ARTIO_SUCCESS;
}

void artio_file_set_endian_swap_tag_i(artio_fh *handle) {
	handle->mode |= ARTIO_MODE_ENDIAN_SWAP;
}

#endif /* MPI */
