from __future__ import absolute_import
from yt.mods import *
import matplotlib
import pylab
from .output_tests import SingleOutputTest, YTDatasetTest, create_test
from yt.analysis_modules.halo_finding.api import *
import hashlib
import numpy as np

# Tests the number of halos returned by the HOP halo finder on a dataset
class TestHaloCountHOP(YTDatasetTest):
    threshold = 80.0

    def run(self):
        # Find the halos using vanilla HOP.
        halos = HaloFinder(self.ds, threshold=self.threshold, dm_only=False)
        # We only care about the number of halos.
        self.result = len(halos)
                    
    def compare(self, old_result):
        # The new value should be identical to the old one.
        self.compare_value_delta(self.result, old_result, 0)

    def plot(self):
        return []

# Tests the number of halos returned by the FOF halo finder on a dataset
class TestHaloCountFOF(YTDatasetTest):
    link = 0.2
    padding = 0.02

    def run(self):
        # Find the halos using FOF.
        halos = FOFHaloFinder(self.ds, link=self.link, dm_only=False, 
                               padding=self.padding)
        # We only care about the number of halos.
        self.result = len(halos)
                    
    def compare(self, old_result):
        # The new value should be identical to the old one.
        self.compare_value_delta(self.result, old_result, 0)

    def plot(self):
        return []

# Tests the number of halos returned by the Parallel HOP halo finder on a 
# dataset
class TestHaloCountPHOP(YTDatasetTest):
    threshold = 80.0

    def run(self):
        # Find the halos using parallel HOP.
        halos = parallelHF(self.ds, threshold=self.threshold, dm_only=False)
        # We only care about the number of halos.
        self.result = len(halos)
                    
    def compare(self, old_result):
        # The new value should be identical to the old one.
        self.compare_value_delta(self.result, old_result, 0)

    def plot(self):
        return []

class TestHaloComposition(YTDatasetTest):
    threshold=80.0
    
    def run(self):
        # Find the halos using vanilla HOP.
        halos = HaloFinder(self.ds, threshold=self.threshold, dm_only=False)
        # The result is a list of the particle IDs, stored
        # as sets for easy comparison.
        IDs = []
        for halo in halos:
            IDs.append(set(halo["particle_index"]))
        self.result = IDs
    
    def compare(self, old_result):
        # All the sets should be identical.
        pairs = zip(self.result, old_result)
        for pair in pairs:
            if len(pair[0] - pair[1]) != 0:
                return False
        return True
    
# Tests the content of the halos returned by the HOP halo finder on a dataset 
# by comparing the hash of the arrays of all the particles contained in each
# halo.  Evidently breaks on parallel runtime.  DO NOT USE.
class TestHaloCompositionHashHOP(YTDatasetTest):
    threshold=80.0
    
    def run(self):
        # Find the halos using vanilla HOP.
        halos = HaloFinder(self.ds, threshold=self.threshold, dm_only=False)
        # The result is a flattened array of the arrays of the particle IDs for
        # each halo
        IDs = []
        for halo in halos:
            IDs.append(halo["particle_index"])
        IDs = np.concatenate(IDs)
        self.result = IDs
    
    def compare(self, old_result):
        # All the lists of arrays should be identical.  To check this
        # faster, we take the 256-bit hash of these lists and compare them
        result_hash = hashlib.sha256(self.result.tostring()).hexdigest()
        old_result_hash = hashlib.sha256(old_result.tostring()).hexdigest()
        if result_hash == old_result_hash:
            return True
        else:
            return False

# Tests the content of the halos returned by the FOF halo finder on a dataset 
# by comparing the hash of the arrays of all the particles contained in each
# halo.  Evidently breaks on parallel runtime.  DO NOT USE.
class TestHaloCompositionHashFOF(YTDatasetTest):
    link = 0.2
    padding = 0.02
    
    def run(self):
        # Find the halos using vanilla FOF.
        halos = FOFHaloFinder(self.ds, link=self.link, dm_only=False, 
                               padding=self.padding)
        # The result is a flattened array of the arrays of the particle IDs for
        # each halo
        IDs = []
        for halo in halos:
            IDs.append(halo["particle_index"])
        IDs = np.concatenate(IDs)
        self.result = IDs
    
    def compare(self, old_result):
        # All the lists of arrays should be identical.  To check this
        # faster, we take the 256-bit hash of these lists and compare them
        result_hash = hashlib.sha256(self.result.tostring()).hexdigest()
        old_result_hash = hashlib.sha256(old_result.tostring()).hexdigest()
        if result_hash == old_result_hash:
            return True
        else:
            return False

# Tests the content of the halos returned by the Parallel HOP halo finder on a 
# dataset by comparing the hash of the arrays of all the particles contained 
# in each halo.  Evidently breaks on parallel runtime.  DO NOT USE.
class TestHaloCompositionHashPHOP(YTDatasetTest):
    threshold=80.0
    
    def run(self):
        # Find the halos using parallel HOP.
        halos = parallelHF(self.ds, threshold=self.threshold, dm_only=False)
        # The result is a flattened array of the arrays of the particle IDs for
        # each halo
        IDs = []
        for halo in halos:
            IDs.append(halo["particle_index"])
        IDs = np.concatenate(IDs)
        self.result = IDs
    
    def compare(self, old_result):
        # All the lists of arrays should be identical.  To check this
        # faster, we take the 256-bit hash of these lists and compare them
        result_hash = hashlib.sha256(self.result.tostring()).hexdigest()
        old_result_hash = hashlib.sha256(old_result.tostring()).hexdigest()
        if result_hash == old_result_hash:
            return True
        else:
            return False
