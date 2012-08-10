"""
Default tests

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2010-2011 Matthew Turk.  All Rights Reserved.

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

from yt.mods import *
from output_tests import YTStaticOutputTest, create_test

class TestFieldStatistics(YTStaticOutputTest):

    tolerance = None

    def run(self):
        # We're going to calculate the field statistics for every single field.
        results = {}
        for field in self.pf.h.field_list:
            # Do it here so that it gets wiped each iteration
            dd = self.pf.h.all_data() 
            results[field] = (dd[field].std(),
                              dd[field].mean(),
                              dd[field].min(),
                              dd[field].max())
        self.result = results

    def compare(self, old_result):
        for field in sorted(self.result):
            for i in range(4):
                oi = old_result[field][i]
                ni = self.result[field][i]
                self.compare_value_delta(oi, ni, self.tolerance)

class TestAllProjections(YTStaticOutputTest):

    tolerance = None

    def run(self):
        results = {}
        for field in self.pf.h.field_list:
            if self.pf.field_info[field].particle_type: continue
            results[field] = []
            for ax in range(3):
                t = self.pf.h.proj(field, ax)
                results[field].append(t.field_data)
        self.result = results

    def compare(self, old_result):
        for field in sorted(self.result):
            for p1, p2 in zip(self.result[field], old_result[field]):
                self.compare_data_arrays(p1, p2, self.tolerance)
