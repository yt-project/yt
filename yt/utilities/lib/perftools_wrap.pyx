"""
Turn on and off perftools profiling



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

# For more info:
# https://pygabriel.wordpress.com/2010/04/14/profiling-python-c-extensions/

# prof.pyx
cdef extern from "google/profiler.h":
    void ProfilerStart( char* fname )
    void ProfilerStop()

def profiler_start(fname):
    ProfilerStart(<char *>fname)

def profiler_stop():
    ProfilerStop()

