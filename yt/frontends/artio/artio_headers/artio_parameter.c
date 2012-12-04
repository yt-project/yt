/*
 * artio_parameter.c
 *
 *  Created on: Jun 8, 2010
 *      Author: Yongen Yu
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

#include "artio.h"
#include "artio_internal.h"

size_t artio_type_size(int type) {
	size_t t_len=0;

	switch (type) {
		case ARTIO_TYPE_STRING:
		case ARTIO_TYPE_CHAR:
			t_len = 1;
			break;
	case ARTIO_TYPE_INT:
			t_len = sizeof(int32_t);
			break;
	case ARTIO_TYPE_FLOAT:
			t_len = sizeof(float);
			break;
	case ARTIO_TYPE_DOUBLE:
			t_len = sizeof(double);
			break;
	case ARTIO_TYPE_LONG:
			t_len = sizeof(int64_t);
			break;
	default:
			t_len = (size_t)-1;
			break;
	}

	return t_len;
}

void artio_parameter_list_init(list *param_list) {
	param_list->head = NULL;
	param_list->tail = NULL;
	param_list->cursor = NULL;
	param_list->iterate_flag = 0;
}

int artio_parameter_read(artio_fh handle, list * param_list) {
	param * param_item;
	int i;
	int length, re;
	int t_len;
	int32_t endian_tag;

    /* endian check */
	re = artio_file_fread(handle, &endian_tag, 1, ARTIO_TYPE_INT);
	if ( re != ARTIO_SUCCESS ) {
		return ARTIO_ERR_PARAM_CORRUPTED;
	}

	if ( endian_tag != ARTIO_ENDIAN_MAGIC ) {
		artio_int_swap( &endian_tag, 1 );
		if ( endian_tag == ARTIO_ENDIAN_MAGIC ) {
			artio_set_endian_swap_tag(handle);
		} else {
			return ARTIO_ERR_PARAM_CORRUPTED_MAGIC;
		}
	}

	re = artio_file_fread(handle, &length, 1, ARTIO_TYPE_INT);
	if ( re != ARTIO_SUCCESS ) {
		return ARTIO_ERR_PARAM_CORRUPTED;
	}

	for ( i = 0; i < length; i++ ) {
		param_item = (param *)malloc(sizeof(param));
		artio_file_fread(handle, &param_item->key_length, 1, ARTIO_TYPE_INT);
		artio_file_fread(handle, param_item->key, param_item->key_length, ARTIO_TYPE_CHAR);
		param_item->key[param_item->key_length] = 0;

		artio_file_fread(handle, &param_item->val_length, 1, ARTIO_TYPE_INT);
		artio_file_fread(handle, &param_item->type, 1, ARTIO_TYPE_INT);

		t_len = artio_type_size(param_item->type);
		param_item->value = (char *)malloc(param_item->val_length * t_len);

		re = artio_file_fread(handle, param_item->value, 
				param_item->val_length, param_item->type);
		if ( re != ARTIO_SUCCESS ) {
			return ARTIO_ERR_PARAM_CORRUPTED;
		}

		param_item->next = NULL;
		if (NULL == param_list->tail) {
			param_list->tail = param_item;
			param_list->head = param_item;
		} else {
			param_list->tail->next = param_item;
			param_list->tail = param_item;
		}
	}

	return ARTIO_SUCCESS;
}

int artio_parameter_write(artio_fh handle, list * param_list) {
	param * item;

	/* retain a number for endian check */
	int32_t endian_tag = ARTIO_ENDIAN_MAGIC;
	int32_t length = 0;

	item = param_list->head;
	while (NULL != item) {
		length++;
		item = item->next;
	}

	artio_file_fwrite(handle, &endian_tag, 
			1, ARTIO_TYPE_INT);
	artio_file_fwrite(handle, &length,
			1, ARTIO_TYPE_INT);

	item = param_list->head;
	while (NULL != item) {
		artio_file_fwrite(handle, &item->key_length, 
				1, ARTIO_TYPE_INT);
		artio_file_fwrite(handle, item->key, 
				item->key_length, ARTIO_TYPE_CHAR);
		artio_file_fwrite(handle, &item->val_length, 
				1, ARTIO_TYPE_INT);
		artio_file_fwrite(handle, &item->type, 
				1, ARTIO_TYPE_INT);
		artio_file_fwrite(handle, item->value, 
				item->val_length, item->type);
		item = item->next;
	}

	return ARTIO_SUCCESS;
}

int artio_parameter_iterate( artio_file handle, char *key, int *type, int *length ) {
	param *item;
	list *param_list = &handle->param_list;	

	if ( param_list->iterate_flag == 0 ) {
		param_list->cursor = param_list->head;	
		param_list->iterate_flag = 1;
	}

	if ( param_list->cursor == NULL ) {
		param_list->iterate_flag = 0;
		return ARTIO_PARAMETER_EXHAUSTED;
	}
 
	item = param_list->cursor;
	strncpy( key, item->key, 64 );
	*type = item->type;
	*length = param_array_length(item);

	param_list->cursor = item->next;
	return ARTIO_SUCCESS;
}

param *artio_parameter_list_search(list * param_list, char *key ) {
	param * item = param_list->head;
	while ( NULL != item && strcmp(item->key, key) ) {
		item = item->next;
	}
	return item;
}

int artio_parameter_list_insert(list * param_list, char * key, 
		int length, void *value, int type) {
	int key_len, val_len = 0;
	param * item;

	if ( length <= 0 ) {
		return ARTIO_ERR_PARAM_LENGTH_INVALID;
	}

	item = artio_parameter_list_search(param_list, key);
	if (NULL != item) {
		return ARTIO_ERR_PARAM_DUPLICATE;
	}

	/* create the list node */
	item = (param *)malloc(sizeof(param));
	key_len = strlen(key);
	item->key_length = key_len;
	strcpy(item->key, key);
	item->val_length = length;
	item->type = type;
	val_len = artio_type_size(type);
	item->value = (char *)malloc(length * val_len);

	memcpy(item->value, value, length * val_len);

	item->next = NULL;
	/* add to the list */
	if (NULL == param_list->tail) {
		param_list->tail = item;
		param_list->head = item;
	} else {
		param_list->tail->next = item;
		param_list->tail = item;
	}

	return ARTIO_SUCCESS;
}

int artio_parameter_list_unpack(list *param_list, char *key, int length,
		void *value, int type) {
	size_t t_len;

	param *item = artio_parameter_list_search(param_list, key);

	if (item != NULL) {
		if (length != item->val_length ) {
			return ARTIO_ERR_PARAM_LENGTH_MISMATCH;
		} else if ( type != item->type ) {
			return ARTIO_ERR_PARAM_TYPE_MISMATCH;
		} else {
			t_len = artio_type_size(type);
			memcpy(value, item->value, item->val_length * t_len);
		}
	} else {
		return ARTIO_ERR_PARAM_NOT_FOUND;
	}

	return ARTIO_SUCCESS;
}

int param_array_length( param *item ) {
	int i, length;

	if ( item->type == ARTIO_TYPE_STRING ) {
		length = 0;
		for ( i = 0; i < item->val_length; i++ ) {
			if ( item->value[i] == '\0' ) {
				length++;
			}
		}
	} else {
		length = item->val_length;
	}

	return length;
}

int artio_parameter_get_array_length(artio_file handle, char * key, int *length) {
	param * item = artio_parameter_list_search(&handle->param_list, key);

	if (item != NULL) {
		*length = param_array_length(item);
	} else {
		return ARTIO_ERR_PARAM_NOT_FOUND;
	}

	return ARTIO_SUCCESS;
}

int artio_parameter_free_list(list * param_list) {
	param * tmp;
	param * param_item = param_list->head;

	if ( param_list == NULL ) {
		return ARTIO_ERR_INVALID_PARAM_LIST;
	}

	while (NULL != param_item) {
		tmp = param_item;
		param_item = param_item->next;
		free(tmp->value);
		free(tmp);
	}

	param_list->head = NULL;
	param_list->tail = NULL;

	return ARTIO_SUCCESS;
}

int artio_parameter_list_print(list * param_list) {
	int32_t a;
	float b;
	double c;
	int64_t d;

	param * item = param_list->head;
	while (NULL != item) {
		switch ( item->type ) {
			case ARTIO_TYPE_STRING:
				printf("string: key %s %s\n", item->key, item->value);
				break;
			case ARTIO_TYPE_CHAR:
				printf("char: key %s %c\n", item->key, *item->value);
				break;
			case ARTIO_TYPE_INT:
				memcpy(&a, item->value, sizeof(int32_t));
				printf("int: key %s %d\n", item->key, a);
				break;
			case ARTIO_TYPE_FLOAT:
				memcpy(&b, item->value, sizeof(float));
				printf("float: key %s %f\n", item->key, b);
				break;
			case ARTIO_TYPE_DOUBLE:
				memcpy(&c, item->value, sizeof(double));
				printf("double: key %s %f\n", item->key, c);
				break;
			case ARTIO_TYPE_LONG:
				memcpy(&d, item->value, sizeof(int64_t));
				printf("long: %ld\n", d);
				break;
			default:
				printf("unrecognized type %d\n", item->type);
		}
		item = item->next;
	}

	return ARTIO_SUCCESS;
}

void artio_parameter_set_int(artio_file handle, char * key, int32_t value) {
	int32_t tmp = value;
	artio_parameter_set_int_array(handle, key, 1, &tmp);
}

int artio_parameter_get_int(artio_file handle, char * key, int32_t * value) {
	return artio_parameter_get_int_array(handle, key, 1, value);
}

void artio_parameter_set_int_array(artio_file handle, char * key, int length,
		int32_t * value) {
	artio_parameter_list_insert(&handle->param_list, key, length, value,
			ARTIO_TYPE_INT);
}

int artio_parameter_get_int_array(artio_file handle, char * key, int length,
		int32_t * value) {
	return artio_parameter_list_unpack(&handle->param_list, key, length, 
			value, ARTIO_TYPE_INT);
}

void artio_parameter_set_float(artio_file handle, char * key, float value) {
	float tmp = value;
	artio_parameter_set_float_array(handle, key, 1, &tmp);
}

int artio_parameter_get_float(artio_file handle, char * key, float * value) {
	return artio_parameter_get_float_array(handle, key, 1, value);
}

void artio_parameter_set_float_array(artio_file handle, char * key,
		int length, float * value) {
	artio_parameter_list_insert(&handle->param_list, key, length, value,
			ARTIO_TYPE_FLOAT);
}

int artio_parameter_get_float_array(artio_file handle, char * key,
		int length, float * value) {
	return artio_parameter_list_unpack(&handle->param_list, key, length, 
			value, ARTIO_TYPE_FLOAT);
}

void artio_parameter_set_double(artio_file handle, char * key, double value) {
	double tmp = value;
	artio_parameter_set_double_array(handle, key, 1, &tmp);
}

int artio_parameter_get_double(artio_file handle, char * key, double * value) {
	return artio_parameter_get_double_array(handle, key, 1, value);
}

void artio_parameter_set_double_array(artio_file handle, char * key,
		int length, double * value) {
	artio_parameter_list_insert(&handle->param_list, key, length, value,
			ARTIO_TYPE_DOUBLE);
}

int artio_parameter_get_double_array(artio_file handle, char * key,
		int length, double * value) {
	return artio_parameter_list_unpack(&handle->param_list, key, length, 
			value, ARTIO_TYPE_DOUBLE);
}

void artio_parameter_set_long(artio_file handle, char * key, int64_t value) {
	int64_t tmp = value;
	artio_parameter_set_long_array(handle, key, 1, &tmp);
}

int artio_parameter_get_long(artio_file handle, char * key, int64_t * value) {
	return artio_parameter_get_long_array(handle, key, 1, value);
}

void artio_parameter_set_long_array(artio_file handle, char * key,
		int length, int64_t * value) {
	artio_parameter_list_insert(&handle->param_list, key, length, value,
			ARTIO_TYPE_LONG);
}

int artio_parameter_get_long_array(artio_file handle, char * key,
		int length, int64_t * value) {
	return artio_parameter_list_unpack(&handle->param_list, key, length, 
			value, ARTIO_TYPE_LONG);
}

void artio_parameter_set_string(artio_file handle, char *key, char *value) {
	artio_parameter_set_string_array(handle, key, 1, &value);
}

void artio_parameter_set_string_array(artio_file handle, char *key,
		int length, char **value) {
	int i;
	int loc_length;
	char *loc_value;
	char *p;

	for (i = 0, loc_length = 0; i < length; i++) {
		loc_length += strlen(value[i]) + 1;
	}

	loc_value = (char *)malloc(loc_length * sizeof(char));

	for (i = 0, p = loc_value; i < length; i++) {
		strcpy(p, value[i]);
		p += strlen(value[i]) + 1;
	}

	artio_parameter_list_insert(&handle->param_list, key, 
			loc_length, loc_value, ARTIO_TYPE_STRING);
	free(loc_value);
}

int artio_parameter_get_string(artio_file handle, char *key, char *value, int max_length) {
	return artio_parameter_get_string_array(handle, key, 1, &value, max_length);
}

int artio_parameter_get_string_array(artio_file handle, char *key,
		int length, char **value, int max_length) {
	int i;
	char *p;
	int count;

	param *item = artio_parameter_list_search(&handle->param_list, key);

	if (item != NULL) {
		/* count string items in item->value */
		count = 0;
		p = item->value;
		while (p < item->value + item->val_length) {
			p += strlen(p) + 1;
			count++;
		}

		if (count != length) {
			return ARTIO_ERR_PARAM_LENGTH_MISMATCH;
		}

		for (i = 0, p = item->value; i < length; i++) {
			strncpy(value[i], p, max_length-1);
			value[i][max_length-1] = 0;
			p += strlen(p) + 1;
		}
	} else {
		return ARTIO_ERR_PARAM_NOT_FOUND;
	}

	return ARTIO_SUCCESS;
}
