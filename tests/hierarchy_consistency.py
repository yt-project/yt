import numpy as na

from yt.utilities.answer_testing.output_tests import \
    YTStaticOutputTest, RegressionTestException
from yt.funcs import ensure_list


class HierarchyInconsistent(RegressionTestException):
    pass


class HierarchyConsistency(YTStaticOutputTest):
    name = "hierarchy_consistency"

    def run(self):
        self.result = \
            all(g in ensure_list(c.Parent) for g in self.pf.h.grids
                                            for c in g.Children)

    def compare(self, old_result):
        if not(old_result and self.result): raise HierarchyInconsistent()


class GridLocationsProperties(YTStaticOutputTest):
    name = "level_consistency"

    def run(self):
        self.result = dict(grid_left_edge=self.pf.h.grid_left_edge,
                           grid_right_edge=self.pf.h.grid_right_edge,
                           grid_levels=self.pf.h.grid_levels,
                           grid_particle_count=self.pf.h.grid_particle_count,
                           grid_dimensions=self.pf.h.grid_dimensions)

    def compare(self, old_result):
        # We allow now difference between these values
        self.compare_data_arrays(self.result, old_result, 0.0)


class GridRelationshipsChanged(RegressionTestException):
    pass


class GridRelationships(YTStaticOutputTest):

    name = "grid_relationships"

    def run(self):
        self.result = [[p.id for p in ensure_list(g.Parent) \
            if g.Parent is not None]
            for g in self.pf.h.grids]

    def compare(self, old_result):
        if len(old_result) != len(self.result):
            raise GridRelationshipsChanged()
        for plist1, plist2 in zip(old_result, self.result):
            if len(plist1) != len(plist2): raise GridRelationshipsChanged()
            if not all((p1 == p2 for p1, p2 in zip(plist1, plist2))):
                raise GridRelationshipsChanged()


class GridGlobalIndices(YTStaticOutputTest):
    name = "global_startindex"

    def run(self):
        self.result = na.array([g.get_global_startindex()
                                for g in self.pf.h.grids])

    def compare(self, old_result):
        self.compare_array_delta(old_result, self.result, 0.0)
