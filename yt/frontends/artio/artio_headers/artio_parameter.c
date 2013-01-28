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

parameter_list *artio_parameter_list_init() {
	parameter_list *parameters = (parameter_list *)malloc(sizeof(parameter_list));
	if ( parameters != NULL ) {
		parameters->head = NULL;
		parameters->tail = NULL;
		parameters->cursor = NULL;
		parameters->iterate_flag = 0;
	} 
	return parameters;
}

int artio_parameter_read(artio_fh *handle, parameter_list *parameters) {
	parameter * item;
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
		item = (parameter *)malloc(sizeof(parameter));
		if ( item == NULL ) {
			return ARTIO_ERR_MEMORY_ALLOCATION;
		}

		artio_file_fread(handle, &item->key_length, 1, ARTIO_TYPE_INT);
		artio_file_fread(handle, item->key, item->key_length, ARTIO_TYPE_CHAR);
		item->key[item->key_length] = 0;

		artio_file_fread(handle, &item->val_length, 1, ARTIO_TYPE_INT);
		artio_file_fread(handle, &item->type, 1, ARTIO_TYPE_INT);

		t_len = artio_type_size(item->type);
		item->value = (char *)malloc(item->val_length * t_len);

		re = artio_file_fread(handle, item->value, 
				item->val_length, item->type);
		if ( re != ARTIO_SUCCESS ) {
			return ARTIO_ERR_PARAM_CORRUPTED;
		}

		item->next = NULL;
		if (NULL == parameters->tail) {
			parameters->tail = item;
			parameters->head = item;
		} else {
			parameters->tail->next = item;
			parameters->tail = item;
		}
	}

	return ARTIO_SUCCESS;
}

int artio_parameter_write(artio_fh *handle, parameter_list *parameters) {
	parameter * item;

	/* retain a number for endian check */
	int32_t endian_tag = ARTIO_ENDIAN_MAGIC;
	int32_t length = 0;

	item = parameters->head;
	while (NULL != item) {
		length++;
		item = item->next;
	}

	artio_file_fwrite(handle, &endian_tag, 
			1, ARTIO_TYPE_INT);
	artio_file_fwrite(handle, &length,
			1, ARTIO_TYPE_INT);

	item = parameters->head;
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

int artio_parameter_iterate( artio_fileset *handle, 
		char *key, int *type, int *length ) {
	parameter *item;
	parameter_list *parameters = handle->parameters;

	if ( parameters->iterate_flag == 0 ) {
		parameters->cursor = parameters->head;	
		parameters->iterate_flag = 1;
	}

	if ( parameters->cursor == NULL ) {
		parameters->iterate_flag = 0;
		return ARTIO_PARAMETER_EXHAUSTED;
	}
 
	item = parameters->cursor;
	strncpy( key, item->key, 64 );
	*type = item->type;
	*length = artio_parameter_array_length(item);

	parameters->cursor = item->next;
	return ARTIO_SUCCESS;
}

parameter *artio_parameter_list_search(parameter_list * parameters, const char *key ) {
	parameter * item = parameters->head;
	while ( NULL != item && strcmp(item->key, key) ) {
		item = item->next;
	}
	return item;
}

int artio_parameter_list_insert(parameter_list * parameters, const char * key, 
		int length, void *value, int type) {
	int key_len, val_len = 0;
	parameter * item;

	if ( length <= 0 ) {
		return ARTIO_ERR_PARAM_LENGTH_INVALID;
	}

	item = artio_parameter_list_search(parameters, key);
	if (NULL != item) {
		return ARTIO_ERR_PARAM_DUPLICATE;
	}

	/* create the list node */
	item = (parameter *)malloc(sizeof(parameter));
	if ( item == NULL ) {
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}

	key_len = strlen(key);
	item->key_length = key_len;
	strcpy(item->key, key);
	item->val_length = length;
	item->type = type;
	val_len = artio_type_size(type);
	item->value = (char *)malloc(length * val_len);
	if ( item->value == NULL ) {
		free(item);
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}

	memcpy(item->value, value, length * val_len);

	item->next = NULL;
	/* add to the list */
	if (NULL == parameters->tail) {
		parameters->tail = item;
		parameters->head = item;
	} else {
		parameters->tail->next = item;
		parameters->tail = item;
	}

	return ARTIO_SUCCESS;
}

int artio_parameter_list_unpack(parameter_list *parameters, 
		const char *key, int length,
		void *value, int type) {
	size_t t_len;
	parameter *item = artio_parameter_list_search(parameters, key);

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

int artio_parameter_array_length( parameter *item ) {
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

int artio_parameter_get_array_length(artio_fileset *handle, const char * key, int *length) {
	parameter *item = artio_parameter_list_search(handle->parameters, key);

	if (item != NULL) {
		*length = artio_parameter_array_length(item);
	} else {
		return ARTIO_ERR_PARAM_NOT_FOUND;
	}

	return ARTIO_SUCCESS;
}

int artio_parameter_list_free(parameter_list * parameters) {
	parameter * tmp;
	parameter * item;

	if ( parameters != NULL ) {
		item  = parameters->head;
		while (NULL != item) {
			tmp = item;
			item = item->next;
			free(tmp->value);
			free(tmp);
		}

		parameters->head = NULL;
		parameters->tail = NULL;
	
		free( parameters );
	}

	return ARTIO_SUCCESS;
}

int artio_parameter_list_print(parameter_list * parameters) {
	int32_t a;
	float b;
	double c;
	int64_t d;

	parameter * item = parameters->head;
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

int artio_parameter_set_int(artio_fileset *handle, const char * key, int32_t value) {
	int32_t tmp = value;
	return artio_parameter_set_int_array(handle, key, 1, &tmp);
}

int artio_parameter_get_int(artio_fileset *handle, const char * key, int32_t * value) {
	return artio_parameter_get_int_array(handle, key, 1, value);
}

int artio_parameter_set_int_array(artio_fileset *handle, const char * key, int length,
		int32_t * value) {
	return artio_parameter_list_insert(handle->parameters, key, length, value,
			ARTIO_TYPE_INT);
}

int artio_parameter_get_int_array(artio_fileset *handle, const char * key, int length,
		int32_t * value) {
	return artio_parameter_list_unpack(handle->parameters, key, length, 
			value, ARTIO_TYPE_INT);
}

int artio_parameter_set_float(artio_fileset *handle, const char * key, float value) {
	float tmp = value;
	return artio_parameter_set_float_array(handle, key, 1, &tmp);
}

int artio_parameter_get_float(artio_fileset *handle, const char * key, float * value) {
	return artio_parameter_get_float_array(handle, key, 1, value);
}

int artio_parameter_set_float_array(artio_fileset *handle, const char * key,
		int length, float * value) {
	return artio_parameter_list_insert(handle->parameters, key, length, value,
			ARTIO_TYPE_FLOAT);
}

int artio_parameter_get_float_array(artio_fileset *handle, const char * key,
		int length, float * value) {
	return artio_parameter_list_unpack(handle->parameters, key, length, 
			value, ARTIO_TYPE_FLOAT);
}

int artio_parameter_set_double(artio_fileset *handle, const char * key, double value) {
	double tmp = value;
	return artio_parameter_set_double_array(handle, key, 1, &tmp);
}

int artio_parameter_get_double(artio_fileset *handle, const char * key, double * value) {
	return artio_parameter_get_double_array(handle, key, 1, value);
}

int artio_parameter_set_double_array(artio_fileset *handle, const char * key,
		int length, double * value) {
	return artio_parameter_list_insert(handle->parameters, key, length, value,
			ARTIO_TYPE_DOUBLE);
}

int artio_parameter_get_double_array(artio_fileset *handle, const char * key,
		int length, double * value) {
	return artio_parameter_list_unpack(handle->parameters, key, length, 
			value, ARTIO_TYPE_DOUBLE);
}

int artio_parameter_set_long(artio_fileset *handle, const char * key, int64_t value) {
	int64_t tmp = value;
	return artio_parameter_set_long_array(handle, key, 1, &tmp);
}

int artio_parameter_get_long(artio_fileset *handle, const char * key, int64_t * value) {
	return artio_parameter_get_long_array(handle, key, 1, value);
}

int artio_parameter_set_long_array(artio_fileset *handle, const char * key,
		int length, int64_t * value) {
	return artio_parameter_list_insert(handle->parameters, key, length, value,
			ARTIO_TYPE_LONG);
}

int artio_parameter_get_long_array(artio_fileset *handle, const char * key,
		int length, int64_t * value) {
	return artio_parameter_list_unpack(handle->parameters, key, length, 
			value, ARTIO_TYPE_LONG);
}

int artio_parameter_set_string(artio_fileset *handle, const char *key, char *value) {
	return artio_parameter_set_string_array(handle, key, 1, &value);
}

int artio_parameter_set_string_array(artio_fileset *handle, const char *key,
		int length, char **value) {
	int i;
	int loc_length;
	char *loc_value;
	char *p;
	int ret;

	for (i = 0, loc_length = 0; i < length; i++) {
		loc_length += strlen(value[i]) + 1;
	}

	loc_value = (char *)malloc(loc_length * sizeof(char));
	if ( loc_value == NULL ) {
		return ARTIO_ERR_MEMORY_ALLOCATION;
	}

	for (i = 0, p = loc_value; i < length; i++) {
		strcpy(p, value[i]);
		p += strlen(value[i]) + 1;
	}

	ret = artio_parameter_list_insert(handle->parameters, key, 
				loc_length, loc_value, ARTIO_TYPE_STRING);
	free(loc_value);

	return ret;
}

int artio_parameter_get_string(artio_fileset *handle, const char *key, char *value, int max_length) {
	return artio_parameter_get_string_array(handle, key, 1, &value, max_length);
}

int artio_parameter_get_string_array(artio_fileset *handle, const char *key,
		int length, char **value, int max_length) {
	int i;
	char *p;
	int count;

	parameter *item = artio_parameter_list_search(handle->parameters, key);

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
