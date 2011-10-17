from yt.mods import *
import matplotlib; matplotlib.use("Agg")
import pylab
from output_tests import SingleOutputTest, YTStaticOutputTest, create_test

# Tests the number of halos returned by the HOP halo finder on a dataset
class TestHaloCountHOP(YTStaticOutputTest):
    threshold = 80.0

    def run(self):
        # Find the haloes using vanilla HOP.
        haloes = HaloFinder(self.pf, threshold=self.threshold, dm_only=False)
        # We only care about the number of haloes.
        self.result = len(haloes)
                    
    def compare(self, old_result):
        # The new value should be identical to the old one.
        self.compare_value_delta(self.result, old_result, 0)

    def plot(self):
        return []

create_test(TestHaloCountHOP, "halo_count_test_HOP", threshold=80.0)

# Tests the number of halos returned by the FOF halo finder on a dataset
class TestHaloCountFOF(YTStaticOutputTest):
    threshold = 80.0

    def run(self):
        # Find the haloes using FOF.
        haloes = FOFHaloFinder(self.pf, threshold=self.threshold, dm_only=False)
        # We only care about the number of haloes.
        self.result = len(haloes)
                    
    def compare(self, old_result):
        # The new value should be identical to the old one.
        self.compare_value_delta(self.result, old_result, 0)

    def plot(self):
        return []

create_test(TestHaloCountFOF, "halo_count_test_FOF", threshold=80.0)

# Tests the number of halos returned by the Parallel HOP halo finder on a 
# dataset
class TestHaloCountPHOP(YTStaticOutputTest):
    threshold = 80.0

    def run(self):
        # Find the haloes using parallel HOP.
        haloes = parallelHF(self.pf, threshold=self.threshold, dm_only=False)
        # We only care about the number of haloes.
        self.result = len(haloes)
                    
    def compare(self, old_result):
        # The new value should be identical to the old one.
        self.compare_value_delta(self.result, old_result, 0)

    def plot(self):
        return []

create_test(TestHaloCountPHOP, "halo_count_test_PHOP", threshold=80.0)

class TestHaloComposition(YTStaticOutputTest):
    threshold=80.0
    
    def run(self):
        # Find the haloes using vanilla HOP.
        haloes = HaloFinder(self.pf, threshold=self.threshold, dm_only=False)
        # The result is a list of the particle IDs, stored
        # as sets for easy comparison.
        IDs = []
        for halo in haloes:
            IDs.append(set(halo["particle_index"]))
        self.result = IDs
    
    def compare(self, old_result):
        # All the sets should be identical.
        pairs = zip(self.result, old_result)
        for pair in pairs:
            if len(pair[0] - pair[1]) != 0:
                return False
        return True
    
    def plot(self):
        return []

create_test(TestHaloComposition, "halo_composition_test", threshold=80.0)
