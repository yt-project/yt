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
from datetime import datetime as dt
from bisect import insort
import atexit

class PerformanceCounters(object):
    _shared_state = {}
    def __new__(cls, *args, **kwargs):
        self = object.__new__(cls, *args, **kwargs)
        self.__dict__ = cls._shared_state
        return self

    def __init__(self):
        self.counters = defaultdict(lambda: 0.0)
        self.counting = defaultdict(lambda: False)
        self.starttime = defaultdict(lambda: 0)
        self.endtime = defaultdict(lambda: 0)
        self._on = ytcfg.getboolean("yt", "time_functions")
        self.exit()

    def __call__(self, name):
        if not self._on: return
        if self.counting[name]:
            self.counters[name] = time.time() - self.counters[name] 
            self.counting[name] = False
            self.endtime[name] = dt.now()
        else:
            self.counters[name] = time.time()
            self.counting[name] = True
            self.starttime[name] = dt.now()

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
        times = []
        for i in self.counters:
            insort(times, [self.starttime[i], i, 1]) # 1 for 'on'
            if not self.counting[i]:
                insort(times, [self.endtime[i], i, 0]) # 0 for 'off'
        #print times
        shift = 0
        multi = 5
        max = 20
        endline = ""
        for i in times:
            # A starting entry
            if i[2] == 1:
                shift += 1
            # An ending entry
            if i[2] == 0:
                # if shift > 1, this is a nested entry, so we want to record
                # this line to be printed later when the top level finish entry
                # is encountered.
                if shift > 1:
                    if self.counting[i[1]]:
                        endline = "%s%i : %s : still running\n%s" % (" "*shift*multi,shift, i[1],endline)
                    else:
                        endline = "%s%i : %s : %0.3e\n%s" % (" "*shift*multi,shift, i[1], self.counters[i[1]], endline)
                    shift -= 1
                # A top level entry.
                else:
                    if self.counting[i[1]]:
                        line = "%i : %s : still running\n%s" % (shift, i[1],endline)
                    else:
                        line = "%i : %s : %0.3e\n%s" % (shift, i[1], self.counters[i[1]],endline)
                    shift -= 1
                    endline = ""
                    print line

    def exit(self):
        if self._on:
            atexit.register(self.print_stats)

yt_counters = PerformanceCounters()
time_function = yt_counters.call_func
