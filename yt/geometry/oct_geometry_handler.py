"""
Octree geometry handler




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import h5py
import numpy as np
import string, re, gc, time
from yt.extern.six.moves import cPickle
from yt.extern.six.moves import zip as izip
import weakref

from itertools import chain

from yt.funcs import *
from yt.utilities.logger import ytLogger as mylog
from yt.arraytypes import blankRecordArray
from yt.config import ytcfg
from yt.fields.field_info_container import NullFunc
from yt.geometry.geometry_handler import Index, YTDataChunk
from yt.utilities.definitions import MAXLEVEL
from yt.utilities.io_handler import io_registry
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface

from yt.data_objects.data_containers import data_object_registry

class OctreeIndex(Index):
    """The Index subclass for oct AMR datasets"""
    def _setup_geometry(self):
        mylog.debug("Initializing Octree Geometry Handler.")
        self._initialize_oct_handler()
        self._initialize_level_stats()

    def get_smallest_dx(self):
        """
        Returns (in code units) the smallest cell size in the simulation.
        """
        return (self.dataset.domain_width /
                (self.dataset.domain_dimensions * 2**(self.max_level))).min()

    def convert(self, unit):
        return self.dataset.conversion_factors[unit]

    def _initialize_level_stats(self):
        levels=sum([dom.level_count for dom in self.domains])
        desc = {'names': ['numcells','level'],
                'formats':['Int64']*2}
        max_level=self.dataset.max_level+1
        self.level_stats = blankRecordArray(desc, max_level)
        self.level_stats['level'] = [i for i in range(max_level)]
        self.level_stats['numcells'] = [0 for i in range(max_level)]
        for level in range(self.dataset.min_level+1):
            self.level_stats[level+1]['numcells']=2**(level*self.dataset.dimensionality)
        for level in range(self.max_level+1):
            self.level_stats[level+self.dataset.min_level+1]['numcells'] = levels[level]


    def print_stats(self):
        """
        Prints out (stdout) relevant information about the simulation
        """
        header = "%3s\t%14s\t%14s" % ("level", "# cells","# cells^3")
        print(header)
        print("%s" % (len(header.expandtabs())*"-"))
        for level in range(self.dataset.max_level+1):
            print("% 3i\t% 14i\t% 14i" % \
                  (level,
                   self.level_stats['numcells'][level],
                   np.ceil(self.level_stats['numcells'][level]**(1./3))))
        print("-" * 46)
        print("   \t% 14i" % (self.level_stats['numcells'].sum()))
        print("\n")

        dx = self.get_smallest_dx()
        try:
            print("z = %0.8f" % (self["CosmologyCurrentRedshift"]))
        except:
            pass
        print("t = %0.8e = %0.8e s = %0.8e years" % \
            (self.ds.current_time.in_units("code_time"),
             self.ds.current_time.in_units("s"),
             self.ds.current_time.in_units("yr")))
        print("\nSmallest Cell:")
        u=[]
        for item in ("Mpc", "pc", "AU", "cm"):
            print("\tWidth: %0.3e %s" % (dx.in_units(item), item))

