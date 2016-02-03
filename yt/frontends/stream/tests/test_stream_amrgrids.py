import numpy as np
from yt.utilities.exceptions import YTIntDomainOverflow

from yt import load_amr_grids, ProjectionPlot

from yt.testing import assert_raises

def test_qt_overflow():
    grid_data = []

    grid_dict = {}

    grid_dict['left_edge'] = [-1.0, -1.0, -1.0]
    grid_dict['right_edge'] = [1.0, 1.0, 1.0]
    grid_dict['dimensions'] = [8, 8, 8]
    grid_dict['level'] = 0

    grid_dict['density'] = np.ones((8,8,8))

    grid_data.append(grid_dict)

    domain_dimensions = np.array([8, 8, 8])

    spf = load_amr_grids(grid_data, domain_dimensions)

    def make_proj():
        p = ProjectionPlot(spf, 'x', ["density"], center='c', origin='native')
        return p
    yield assert_raises, YTIntDomainOverflow, make_proj

def test_refine_by():
    grid_data = []
    ref_by = 4
    lo = 0.0
    hi = 1.0
    fine_grid_width = (hi - lo) / ref_by
    for level in range(2):
        grid_dict = {}

        grid_dict['left_edge'] = [0.0 + 0.5*fine_grid_width*level]*3
        grid_dict['right_edge'] = [1.0 - 0.5*fine_grid_width*level]*3
        grid_dict['dimensions'] = [8, 8, 8]
        grid_dict['level'] = level

        grid_dict['density'] = np.ones((8,8,8))

        grid_data.append(grid_dict)

    domain_dimensions = np.array([8, 8, 8])

    load_amr_grids(grid_data, domain_dimensions, refine_by=ref_by)
