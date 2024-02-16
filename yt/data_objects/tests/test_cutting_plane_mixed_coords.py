import numpy as np

import yt
from yt.testing import fake_amr_ds


def test_cutting_plane_mixed_coords():
    ds = fake_amr_ds(geometry="spherical")
    normal = np.array([0.0, 0.0, 1.0])
    plane_center = np.array([0.0, 0.0, 0.5])
    slc = ds.cutting_mixed(normal, plane_center)
    frb = slc.to_frb(2.0, 800)
    bvals = frb[("index", "r")]
    mask = frb.get_mask(("index", "r"))

    assert np.allclose(bvals[mask].min().d, plane_center[2], atol=0.02)


def test_cutting_plane_mixed_coords_answer():
    # TODO: make this an answer test
    ds = yt.load_sample("KeplerianDisk")
    normal = np.array([0.0, 1.0, 0.0])
    plane_center = np.array([0.0, 0.0, 0.0])
    slc = ds.cutting_mixed(normal, plane_center)
    frb = slc.to_frb(8.0, 400)
    bvals = frb[("athena_pp", "dens")]
    mask = frb.get_mask(("athena_pp", "dens"))
    bvals[~mask] = np.nan

    # make a plot or something.
