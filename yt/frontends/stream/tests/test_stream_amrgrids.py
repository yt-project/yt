from copy import deepcopy
import numpy as np
from yt.utilities.exceptions import \
    YTIllDefinedAMR, \
    YTIntDomainOverflow

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
    assert_raises(YTIntDomainOverflow, make_proj)

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

def test_validation():
    dims = np.array([4, 2, 4])
    grid_data = [
        dict(left_edge = [0.0, 0.0, 0.0],
             right_edge = [1.0, 1.0, 1.0],
             level = 0,
             dimensions = dims),
        dict(left_edge = [0.25, 0.25, 0.25],
             right_edge = [0.75, 0.75, 0.75],
             level = 1,
             dimensions = dims),
       ]
    bbox = np.array([[0, 1], [0, 1], [0, 1]])
    def load_grids():
        load_amr_grids(grid_data, dims, bbox=bbox, periodicity=(0, 0, 0),
                       length_unit=1.0, refine_by=2)
    assert_raises(YTIllDefinedAMR, load_grids)


def test_load_despite_rounding_errors_during_grid_construction():
    # from GH issue 2665
    grid_data = [
    dict(left_edge =[0.0, 0.0, 0.0],
         right_edge=[0.045,0.045,0.045],
         level=0,
         dimensions=[16,16,16]),
    dict(left_edge =[0.0225, 0.0225, 0.0225],
         right_edge=[0.045, 0.045, 0.045],
         level=1,
         dimensions=[16,16,16])
    ]

    for g in grid_data:
        g["density"] = np.random.random(g["dimensions"]) * 2 ** g["level"]

    ill_grid_data = deepcopy(grid_data)
    ill_grid_data[0].update({"left_edge": [1e-1, 0.0, 0.0]})

    bbox = np.array([[-0.18,0.18],[-0.18,0.18],[-0.18,0.18]])
    dimensions = [128,128,128]

    # this should not raise any exception
    load_amr_grids(grid_data, dimensions, bbox)

    # but this should
    assert_raises(YTIllDefinedAMR, load_amr_grids, ill_grid_data, dimensions, bbox)