from yt.utilities.exceptions import \
    YTBoundsDefinitionError

from yt.testing import \
    fake_random_ds, requires_file
from numpy.testing import \
    assert_raises
import yt
from numpy.testing import assert_array_less
import numpy as np


def test_cic_deposit():
    ds = fake_random_ds(64, nprocs = 8, particles=64**3)
    my_reg = ds.arbitrary_grid(ds.domain_left_edge, ds.domain_right_edge,
            dims=[1, 800, 800])
    f = ("deposit", "all_cic")
    assert_raises(YTBoundsDefinitionError, my_reg.__getitem__, f)


output_00080 = "output_00080/info_00080.txt"
@requires_file(output_00080)
def test_mesh_deposition():
    ds = yt.load(output_00080)
    ds.add_deposited_mesh_field(('index', 'x'), ptype='all')
    ds.add_deposited_mesh_field(('index', 'dx'), ptype='all')

    dx = ds.r['all', 'cell_index_dx']
    x = ds.r['all', 'cell_index_x']
    xpart = ds.r['all', 'particle_position_x']

    dist = (x - xpart) / dx

    upper = np.ones_like(dx) / 2

    # Some particles are out of their domain, skip them
    mask = np.isfinite(dist)

    assert_array_less(dist[mask], upper[mask])
    assert_array_less(-dist[mask], upper[mask])
