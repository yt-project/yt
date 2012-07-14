"""
libconfig (http://www.hyperrealm.com/libconfig/) wrapper in Cython.  Exposes a
minimum of the API.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2011 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from libc.stdio cimport fopen, fclose, FILE

cdef extern from "libconfig.h":
    cdef struct config_t:
        pass
    cdef struct config_setting_t:
        pass

    enum: CONFIG_TRUE
    enum: CONFIG_FALSE
    enum: CONFIG_TYPE_INT
    enum: CONFIG_TYPE_INT64
    enum: CONFIG_TYPE_FLOAT
    enum: CONFIG_TYPE_STRING
    enum: CONFIG_TYPE_BOOL
    enum: CONFIG_TYPE_ARRAY
    enum: CONFIG_TYPE_LIST
    enum: CONFIG_TYPE_GROUP

    void config_init (config_t *config)
    void config_destroy (config_t *config)
    int config_read (config_t * config, FILE * stream)
    int config_read_file (config_t * config, char * filename)
    int config_read_string (config_t * config, char * str)
    void config_write (config_t * config, FILE * stream)
    int config_write_file (config_t * config, char * filename)

    char * config_error_text (config_t * config)
    char * config_error_file (config_t * config)
    int config_error_line (config_t * config)

    config_setting_t * config_lookup (config_t * config, char * path)
    int config_lookup_int (config_t * config, char * path, int * value)
    int config_lookup_int64 (config_t * config, char * path, long long * value)
    int config_lookup_float (config_t * config, char * path, double * value)
    int config_lookup_bool (config_t * config, char * path, int * value)
    int config_lookup_string (config_t * config, char * path, char ** value)

    char * config_setting_name (config_setting_t * setting)

    int config_setting_get_int (config_setting_t * setting)
    long long config_setting_get_int64 (config_setting_t * setting)
    double config_setting_get_float (config_setting_t * setting)
    int config_setting_get_bool (config_setting_t * setting)
    char * config_setting_get_string (config_setting_t * setting)

    int config_setting_type (config_setting_t * setting)

    int config_setting_length (config_setting_t * setting)
    config_setting_t * config_setting_get_elem (
        config_setting_t * setting, unsigned int idx)

cdef class libconfigSetting:
    cdef config_setting_t *setting

    def keys(self):
        cdef int tt = config_setting_type(self.setting)
        cdef char *setting_name
        cdef config_setting_t *sub_setting
        if tt != CONFIG_TYPE_GROUP: return []
        nelem = config_setting_length(self.setting)
        tr = []
        for i in xrange(nelem):
            sub_setting = config_setting_get_elem(self.setting, i)
            setting_name = config_setting_name(sub_setting)
            tr.append(setting_name)
        return tr

cdef class libconfigConfiguration:
    cdef config_t cfg
    def __cinit__(self):
        config_init(&self.cfg)

    def __dealloc__(self):
        config_destroy(&self.cfg)

    def read_file(self, char *fn):
        cdef int status, error_line
        cdef char *error_text, *error_file
        status = config_read_file(&self.cfg, fn)
        if status == CONFIG_FALSE:
            error_text = config_error_text(&self.cfg)
            error_file = config_error_file(&self.cfg)
            error_line = config_error_line(&self.cfg)
            if error_file != NULL:
                print "Error encounters in %s at line %s" % (
                        error_file, error_line)
            print error_text

    def read_string(self, char *s):
        config_read_string(&self.cfg, s)

    def __getitem__(self, key):
        cdef config_setting_t *setting = config_lookup(&self.cfg, key)
        ms = libconfigSetting()
        ms.setting = setting
        rv = self.get_setting(ms)
        if rv == None: raise KeyError(key)
        return rv

    def get_setting(self, libconfigSetting ms):
        cdef config_setting_t *setting = ms.setting
        if setting == NULL:
            return None
        cdef int tt = config_setting_type(setting)
        cdef int kv_int
        cdef long long kv_int64
        cdef double kv_float
        cdef int kv_bool
        cdef char* kv_string
        cdef int nelem
        cdef config_setting_t *sub_setting
        if tt == CONFIG_TYPE_INT:
            kv_int = config_setting_get_int(setting)
            return kv_int
        elif tt == CONFIG_TYPE_INT64:
            kv_int64 = config_setting_get_int64(setting)
            return kv_int64
        elif tt == CONFIG_TYPE_FLOAT:
            kv_float = config_setting_get_float(setting)
            return kv_float
        elif tt == CONFIG_TYPE_BOOL:
            kv_bool = config_setting_get_bool(setting)
            return kv_bool
        elif tt == CONFIG_TYPE_STRING:
            kv_string = config_setting_get_string(setting)
            return kv_string
        elif tt == CONFIG_TYPE_ARRAY or tt == CONFIG_TYPE_LIST:
            nelem = config_setting_length(setting)
            tr = []
            for i in xrange(nelem):
                sub_ms = libconfigSetting()
                sub_setting = config_setting_get_elem(setting, i)
                sub_ms.setting = sub_setting
                tr.append(self.get_setting(sub_ms))
            return tr
        elif tt == CONFIG_TYPE_GROUP:
            group = libconfigSetting()
            group.setting = setting
            return group
        return None
