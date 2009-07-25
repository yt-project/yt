"""
Minimalist performance counting for yt

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2009 Matthew Turk.  All Rights Reserved.

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

from yt.config import ytcfg
from yt.funcs import *
import time

class PerformanceCounters(object):
    _shared_state = {}
    def __new__(cls, *args, **kwargs):
        self = object.__new__(cls, *args, **kwargs)
        self.__dict__ = cls._shared_state
        return self

    def __init__(self):
        self.counters = defaultdict(lambda: 0.0)
        self.counting = defaultdict(lambda: False)
        self._on = ytcfg.getboolean("yt", "time_functions")

    def __call__(self, name):
        if not self._on: return
        if self.counting[name]:
            self.counters[name] = time.time() - self.counters[name] 
            self.counting[name] = False
        else:
            self.counters[name] = time.time()
            self.counting[name] = True

    def call_func(self, func):
        if not self._on: return func
        @wraps(func)
        def func_wrapper(*args, **kwargs):
            self(func.func_name)
            func(*args, **kwargs)
            self(func.func_name)
        return func_wrapper

    def print_stats(self):
        print "Current counter status:\n"
        for i in sorted(self.counters):
            print "% 30s:" % (i),
            if self.counting[i]: print "still running"
            else: print "%0.3e" % (self.counters[i])

yt_counters = PerformanceCounters()
time_function = yt_counters.call_func
