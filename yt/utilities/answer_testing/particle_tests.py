import matplotlib
from yt.mods import *
import pylab
from output_tests import SingleOutputTest, YTStaticOutputTest, create_test

class TestParticleUniqueIDs(YTStaticOutputTest):

    def run(self):
        # Test to make sure that all the particles have unique IDs.
        all = self.pf.h.all_data()
        IDs = all["particle_index"]
        # Make sure the order is the same every time.
        IDs = IDs[IDs.argsort()]
        self.result = IDs
                    
    def compare(self, old_result):
        # Two things: there should be no repeats in either the new or
        # the old, and the two sets should be the same.
        if len(old_result) != len(set(old_result)): return False
        if len(self.result) != len(set(self.result)): return False
        if (self.result != old_result).all(): return False
        return True

    def plot(self):
        return []

create_test(TestParticleUniqueIDs, "particle_unique_ids_test")

class TestParticleExtrema(YTStaticOutputTest):

    def run(self):
        # Tests to make sure there are no particle positions aren't changing
        # drastically. This is very unlikely to be a problem.
        all = self.pf.h.all_data()
        min = np.empty(3,dtype='float64')
        max = min.copy()
        dims = ["particle_position_x","particle_position_y",
            "particle_position_z"]
        for i in xrange(3):
            min[i] = np.min(all[dims[i]])
            max[i] = np.max(all[dims[i]])
        self.result = (min,max)
    
    def compare(self, old_result):
        min,max = self.result
        old_min, old_max = old_result
        # The extrema should be very similar.
        self.compare_array_delta(min, old_min, 1e-7)
        self.compare_array_delta(max, old_max, 1e-7)
        # Also, the min/max shouldn't be outside the boundaries.
        if (min < self.pf.domain_left_edge).any(): return False
        if (max > self.pf.domain_right_edge).any(): return False
        return True
    
    def plot(self):
        return []

create_test(TestParticleExtrema, "particle_extrema_test")

