/*
 * artio_internal.h
 *
 *  Created on: Apr 9, 2010
 *      Author: Yongen Yu
 *  Renamed/Modified: Nov 18, 2010 - Douglas Rudd
 */

#ifndef __ARTIO_INTERNAL_H__
#define __ARTIO_INTERNAL_H__

#include <stdlib.h>
#include <stdint.h>

#include "artio.h"
#include "artio_endian.h"

#ifndef MIN
#define MIN(x,y)        (((x) < (y)) ? (x): (y))
#endif
#ifndef MAX
#define MAX(x,y)        (((x) > (y)) ? (x): (y))
#endif

#ifdef ARTIO_MPI
#include <mpi.h>
#endif

typedef struct ARTIO_FH *artio_fh;

typedef struct artio_particle_file_struct
{
	artio_fh *ffh;
	int num_particle_files;
	int64_t *file_sfc_index;
	int64_t cache_sfc_begin;
	int64_t cache_sfc_end;
	int64_t *sfc_offset_table;

	/* maintained for consistency and user-error detection */
	int num_species;
	int cur_file;
	int cur_species;
	int cur_particle;
	int64_t cur_sfc;
	int *num_primary_variables;
	int *num_secondary_variables;
	int *num_particles_per_species;
} *artio_particle_file;

typedef struct artio_grid_file_struct
{
	artio_fh *ffh;
	int num_grid_variables;
	int num_grid_files;
	int64_t *file_sfc_index;
	int64_t cache_sfc_begin;
	int64_t cache_sfc_end;
	int64_t *sfc_offset_table;

	int file_max_level;
	/* maintained for consistency and user-error detection */
	int cur_file;
	int cur_num_levels;
	int cur_level;
	int cur_octs;
	int64_t cur_sfc;
	int *octs_per_level;
} *artio_grid_file;

typedef struct param_struct
{
	int key_length;
	char key[64];
	int val_length;
	int type;
	char *value;
	struct param_struct * next;
} param;

typedef struct list_struct
{
	param * head;
	param * tail;
	param * cursor;
	int iterate_flag;
} list;

struct artio_file_struct
{
	char file_prefix[256];
	int endian_swap;
	int open_type;
	int open_mode;
	int rank;
	int num_procs;
	artio_context context;

	int64_t *proc_sfc_index;
	int64_t proc_sfc_begin;
	int64_t proc_sfc_end;
	int64_t num_root_cells;

	list param_list;
	artio_grid_file grid;
	artio_particle_file particle;
};

#define ARTIO_FILESET_READ                  0                                                                       
#define ARTIO_FILESET_WRITE                 1

#define ARTIO_MODE_READ         1
#define ARTIO_MODE_WRITE        2
#define ARTIO_MODE_DIRECT       4
#define ARTIO_MODE_ACCESS       8
#define ARTIO_MODE_ENDIAN_SWAP 16

#define ARTIO_SEEK_SET          0
#define ARTIO_SEEK_CUR          1

artio_fh artio_file_fopen( char * filename, int amode, artio_context context );
int artio_file_fwrite(artio_fh handle, void *buf, int64_t count, int type );
int artio_file_fflush(artio_fh handle);
int artio_file_fseek(artio_fh ffh, int64_t offset, int whence);
int artio_file_fread(artio_fh handle, void *buf, int64_t count, int type );
int artio_file_fclose(artio_fh handle);
void artio_set_endian_swap_tag(artio_fh handle);


#define ARTIO_ENDIAN_MAGIC	0x1234

void artio_parameter_list_init( list *param_list );

param * artio_parameter_list_search(list * param_list, char * key);

int param_array_length( param *item );

int artio_parameter_list_insert(list * param_list, char * key, int length,
		void * value, int type);

int artio_parameter_read(artio_fh handle, list * param_list);

int artio_parameter_write(artio_fh handle, list * param_list);

int artio_parameter_list_print(list * param_list);
int artio_parameter_free_list(list * param_list);
int artio_parameter_list_print(list * param_list);

size_t artio_type_size( int type );
int artio_parameter_list_unpack(list *param_list, char *key, int length, void *value, int type );

#endif /* __ARTIO_INTERNAL_H__ */
