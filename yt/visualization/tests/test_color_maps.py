import os
import shutil
import tempfile
import unittest

import matplotlib.pyplot as plt
import numpy as np
from nose.tools import assert_raises

from yt import make_colormap, show_colormaps
from yt.testing import assert_almost_equal, assert_equal, requires_backend


class TestColorMaps(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.curdir = os.getcwd()
        os.chdir(self.tmpdir)

    def tearDown(self):
        os.chdir(self.curdir)
        shutil.rmtree(self.tmpdir)

    @requires_backend("Agg")
    def test_show_colormaps(self):
        show_colormaps()
        show_colormaps(subset=["jet", "cool"])
        show_colormaps(subset="yt_native", filename="yt_color_maps.png")

        # Test for non-existent color map
        with assert_raises(AttributeError) as ex:
            show_colormaps(subset="unknown", filename="yt_color_maps.png")
        desired = (
            "show_colormaps requires subset attribute to be 'all', "
            "'yt_native', or a list of valid colormap names."
        )
        assert_equal(str(ex.exception), desired)

    @requires_backend("Agg")
    def test_make_colormap(self):
        make_colormap(
            [([0, 0, 1], 10), ([1, 1, 1], 10), ([1, 0, 0], 10)],
            name="french_flag",
            interpolate=False,
        )
        show_colormaps("french_flag")

        cmap = make_colormap(
            [("dred", 5), ("blue", 2.0), ("orange", 0)], name="my_cmap"
        )
        assert_almost_equal(
            cmap["red"][1], np.array([0.00392157, 0.62400345, 0.62400345])
        )

        assert_almost_equal(
            cmap["blue"][2], np.array([0.00784314, 0.01098901, 0.01098901])
        )

        assert_almost_equal(cmap["green"][3], np.array([0.01176471, 0.0, 0.0]))


def test_cmyt_integration():
    for name in ["algae", "bds_highcontrast", "kelp", "arbre", "octarine", "kamae"]:
        cmap = plt.get_cmap(name)
        assert cmap.name == name
        name_r = name + "_r"
        cmap_r = plt.get_cmap(name_r)
        assert cmap_r.name == name_r

    for name in ["algae", "kelp", "arbre", "octarine", "pastel"]:
        cmap = plt.get_cmap("cmyt." + name)
        assert cmap.name == "cmyt." + name
