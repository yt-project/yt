"""
Default tests



"""
from __future__ import absolute_import

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.mods import *
from .output_tests import YTDatasetTest, create_test

class TestFieldStatistics(YTDatasetTest):

    tolerance = None

    def run(self):
        # We're going to calculate the field statistics for every single field.
        results = {}
        for field in self.ds.field_list:
            # Do it here so that it gets wiped each iteration
            dd = self.ds.all_data() 
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

class TestAllProjections(YTDatasetTest):

    tolerance = None

    def run(self):
        results = {}
        for field in self.ds.field_list:
            if self.ds.field_info[field].particle_type: continue
            results[field] = []
            for ax in range(3):
                t = self.ds.proj(field, ax)
                results[field].append(t.field_data)
        self.result = results

    def compare(self, old_result):
        for field in sorted(self.result):
            for p1, p2 in zip(self.result[field], old_result[field]):
                self.compare_data_arrays(p1, p2, self.tolerance)

