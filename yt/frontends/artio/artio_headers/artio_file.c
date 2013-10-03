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

#include <stdlib.h>
#include <stdio.h>

artio_fh *artio_file_fopen( char * filename, int mode, const artio_context *context) {
	artio_fh *fh;
#ifdef ARTIO_DEBUG
	printf( "artio_file_fopen( filename=%s, mode=%d, context=%p )\n", 
			filename, mode, context ); fflush(stdout);
#endif /* ARTIO_DEBUG */
	fh = artio_file_fopen_i(filename,mode,context);
#ifdef ARTIO_DEBUG
	printf(" artio_file_fopen = %p\n", fh ); fflush(stdout);
#endif /* ARTIO_DEBUG */
	return fh;
}

int artio_file_attach_buffer( artio_fh *handle, void *buf, int buf_size ) {
    int status;
#ifdef ARTIO_DEBUG
	printf( "artio_file_attach_buffer( handle=%p, buf=%p, buf_size = %d )\n",
			handle, buf, buf_size ); fflush(stdout);
#endif /* ARTIO_DEBUG */
    status = artio_file_attach_buffer_i(handle,buf,buf_size);
#ifdef ARTIO_DEBUG
	if ( status != ARTIO_SUCCESS ) {
		printf(" artio_file_attach_buffer(%p) = %d\n", handle, status ); fflush(stdout);
	}
#endif /* ARTIO_DEBUG */
    return status;
}

int artio_file_detach_buffer( artio_fh *handle ) {
	int status;
#ifdef ARTIO_DEBUG
	printf( "artio_file_detach_buffer( handle=%p )\n", handle ); fflush(stdout);
#endif /* ARTIO_DEBUG */
	status = artio_file_detach_buffer_i(handle);
#ifdef ARTIO_DEBUG
	if ( status != ARTIO_SUCCESS ) {
		printf( "artio_file_detach_buffer(%p) = %d\n", handle, status ); fflush(stdout);
	}
#endif /* ARTIO_DEBUG */
	return status;
}

int artio_file_fwrite( artio_fh *handle, const void *buf, int64_t count, int type ) {
	int status;
#ifdef ARTIO_DEBUG
	printf( "artio_file_fwrite( handle=%p, buf=%p, count=%ld, type=%d )\n",
			handle, buf, count, type ); fflush(stdout);
#endif /* ARTIO_DEBUG */
	status = artio_file_fwrite_i(handle,buf,count,type);
#ifdef ARTIO_DEBUG
	if ( status != ARTIO_SUCCESS ) {
		printf( "artio_file_fwrite(%p) = %d", handle, status ); fflush(stdout);
	}
#endif /* ARTIO_DEBUG */
	return status;
}

int artio_file_fflush(artio_fh *handle) {
	int status;
#ifdef ARTIO_DEBUG
	printf( "artio_file_fflush( handle=%p )\n", handle ); fflush(stdout);
#endif /* ARTIO_DEBUG */
	status = artio_file_fflush_i(handle);
#ifdef ARTIO_DEBUG
	if ( status != ARTIO_SUCCESS ) {
		printf( "artio_file_fflush(%p) = %d\n", handle, status ); fflush(stdout);
	}
#endif /* ARTIO_DEBUG */
	return status;
}

int artio_file_fread(artio_fh *handle, void *buf, int64_t count, int type ) {
	int status;
#ifdef ARTIO_DEBUG
	printf( "artio_file_fread( handle=%p, buf=%p, count=%ld, type=%d )\n",
			handle, buf, count, type ); fflush(stdout);
#endif /* ARTIO_DEBUG */
	status = artio_file_fread_i(handle,buf,count,type);
#ifdef ARTIO_DEBUG
	if ( status != ARTIO_SUCCESS ) {
		printf( "artio_file_fread(%p) = %d", handle, status );
	}
#endif /* ARTIO_DEBUG */
	return status;
}

int artio_file_ftell(artio_fh *handle, int64_t *offset) {
	int status;
#ifdef ARTIO_DEBUG
	printf( "artio_file_ftell( handle=%p, offset=%p )\n",
		handle, offset ); fflush(stdout);
#endif /* ARTIO_DEBUG */
	status = artio_file_ftell_i(handle,offset);
#ifdef ARTIO_DEBUG
	if ( status != ARTIO_SUCCESS ) {
		printf("artio_file_ftell(%p) = %d\n", handle, status ); fflush(stdout);
	}
#endif /* ARTIO_DEBUG */
	return status;
}

int artio_file_fseek(artio_fh *handle, int64_t offset, int whence ) {
	int status;
#ifdef ARTIO_DEBUG
    printf( "artio_file_fseek( handle=%p, offset=%ld, whence=%d )\n",
        handle, offset, whence ); fflush(stdout);
#endif /* ARTIO_DEBUG */
	status = artio_file_fseek_i(handle,offset,whence);
#ifdef ARTIO_DEBUG
	if ( status != ARTIO_SUCCESS ) {
		printf( "artio_file_fseek(%p) = %d\n", handle, status ); fflush(stdout);
	}
#endif /* ARTIO_DEBUG */
	return status;
}

int artio_file_fclose(artio_fh *handle) {
	int status;
#ifdef ARTIO_DEBUG
	printf( "artio_file_fclose( handle=%p )\n", handle ); fflush(stdout);
#endif /* ARTIO_DEBUG */
	status = artio_file_fclose_i(handle);
#ifdef ARTIO_DEBUG
	if ( status != ARTIO_SUCCESS ) {
		printf( "artio_file_fclose(%p) = %d\n", handle, status ); fflush(stdout);
	}
#endif /* ARTIO_DEBUG */
	return status;
}

void artio_file_set_endian_swap_tag(artio_fh *handle) {
	artio_file_set_endian_swap_tag_i(handle);
}
