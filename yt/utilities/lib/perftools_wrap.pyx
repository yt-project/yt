"""
Turn on and off perftools profiling

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
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

