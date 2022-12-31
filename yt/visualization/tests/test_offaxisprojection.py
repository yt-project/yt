import os
import shutil
import tempfile
import unittest

import numpy as np

from yt.testing import (
    assert_equal,
    assert_fname,
    assert_rel_equal,
    expand_keywords,
    fake_octree_ds,
    fake_random_ds,
)
from yt.visualization.api import OffAxisProjectionPlot, OffAxisSlicePlot
from yt.visualization.image_writer import write_projection
from yt.visualization.volume_rendering.api import off_axis_projection


class TestOffAxisProjection(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()
        cls.curdir = os.getcwd()
        os.chdir(cls.tmpdir)

    @classmethod
    def tearDownClass(cls):
        os.chdir(cls.curdir)
        shutil.rmtree(cls.tmpdir)

    def test_oap(self):
        """Tests functionality of off_axis_projection and write_projection."""

        # args for off_axis_projection
        test_ds = fake_random_ds(64)
        c = test_ds.domain_center
        norm = [0.5, 0.5, 0.5]
        W = test_ds.arr([0.5, 0.5, 1.0], "unitary")
        N = 256
        field = ("gas", "density")
        oap_args = [test_ds, c, norm, W, N, field]

        # kwargs for off_axis_projection
        oap_kwargs = {}
        oap_kwargs["weight"] = (None, "cell_mass")
        oap_kwargs["no_ghost"] = (True, False)
        oap_kwargs["interpolated"] = (False,)
        oap_kwargs["north_vector"] = ((1, 0, 0), (0, 0.5, 1.0))
        oap_kwargs_list = expand_keywords(oap_kwargs)

        # args or write_projection
        fn = "test_%d.png"

        # kwargs for write_projection
        wp_kwargs = {}
        wp_kwargs["colorbar"] = (True, False)
        wp_kwargs["colorbar_label"] = "test"
        wp_kwargs["title"] = "test"
        wp_kwargs["vmin"] = (None,)
        wp_kwargs["vmax"] = (1e3, 1e5)
        wp_kwargs["take_log"] = (True, False)
        wp_kwargs["figsize"] = ((8, 6), [1, 1])
        wp_kwargs["dpi"] = (100, 50)
        wp_kwargs["cmap_name"] = ("cmyt.arbre", "cmyt.kelp")
        wp_kwargs_list = expand_keywords(wp_kwargs)

        # test all off_axis_projection kwargs and write_projection kwargs
        # make sure they are able to be projected, then remove and try next
        # iteration
        for i, oap_kwargs in enumerate(oap_kwargs_list):
            image = off_axis_projection(*oap_args, **oap_kwargs)
            for wp_kwargs in wp_kwargs_list:
                write_projection(image, fn % i, **wp_kwargs)
                assert_fname(fn % i)

        # Test remaining parameters of write_projection
        write_projection(image, "test_2", xlabel="x-axis", ylabel="y-axis")
        assert_fname("test_2.png")

        write_projection(image, "test_3.pdf", xlabel="x-axis", ylabel="y-axis")
        assert_fname("test_3.pdf")

        write_projection(image, "test_4.eps", xlabel="x-axis", ylabel="y-axis")
        assert_fname("test_4.eps")


def test_field_cut_off_axis_octree():
    ds = fake_octree_ds()
    cut = ds.all_data().cut_region('obj[("gas", "density")]>0.5')
    p1 = OffAxisProjectionPlot(ds, [1, 0, 0], ("gas", "density"))
    p2 = OffAxisProjectionPlot(ds, [1, 0, 0], ("gas", "density"), data_source=cut)
    assert_equal(p2.frb[("gas", "density")].min() == 0.0, True)  # Lots of zeros
    assert_equal(
        (p1.frb[("gas", "density")] == p2.frb[("gas", "density")]).all(), False
    )
    p3 = OffAxisSlicePlot(ds, [1, 0, 0], ("gas", "density"))
    p4 = OffAxisSlicePlot(ds, [1, 0, 0], ("gas", "density"), data_source=cut)
    assert_equal(
        (p3.frb[("gas", "density")] == p4.frb[("gas", "density")]).all(), False
    )
    p4rho = p4.frb[("gas", "density")]
    assert_equal(np.nanmin(p4rho[p4rho > 0.0]) >= 0.5, True)


def test_offaxis_moment():
    ds = fake_random_ds(64)

    def _vlos_sq(field, data):
        return data["gas", "velocity_los"] ** 2

    ds.add_field(
        ("gas", "velocity_los_squared"),
        _vlos_sq,
        sampling_type="local",
        units="cm**2/s**2",
    )
    p1 = OffAxisProjectionPlot(
        ds,
        [1, 1, 1],
        [("gas", "velocity_los"), ("gas", "velocity_los_squared")],
        weight_field=("gas", "density"),
        moment=1,
        buff_size=(400, 400),
    )
    p2 = OffAxisProjectionPlot(
        ds,
        [1, 1, 1],
        ("gas", "velocity_los"),
        weight_field=("gas", "density"),
        moment=2,
        buff_size=(400, 400),
    )
    assert_rel_equal(
        np.sqrt(
            p1.frb["gas", "velocity_los_squared"] - p1.frb["gas", "velocity_los"] ** 2
        ),
        p2.frb["gas", "velocity_los"],
        10,
    )
