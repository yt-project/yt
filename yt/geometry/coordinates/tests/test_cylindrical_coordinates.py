# Some tests for the Cylindrical coordinates handler

import numpy as np
import matplotlib.pyplot as plt

from yt.testing import \
    fake_amr_ds, \
    assert_equal, \
    assert_almost_equal, \
    add_noise_fields

from yt.utilities.answer_testing.framework import GenericImageTest
from yt import SlicePlot

# Our canonical tests are that we can access all of our fields and we can
# compute our volume correctly.

def test_cylindrical_coordinates():
    # We're going to load up a simple AMR grid and check its volume
    # calculations and path length calculations.
    ds = fake_amr_ds(geometry="cylindrical")
    axes = ["r", "z", "theta"]
    for i, axis in enumerate(axes):
        dd = ds.all_data()
        fi = ("index", axis)
        fd = ("index", "d%s" % axis)
        ma = np.argmax(dd[fi])
        assert_equal(dd[fi][ma] + dd[fd][ma] / 2.0, ds.domain_right_edge[i].d)
        mi = np.argmin(dd[fi])
        assert_equal(dd[fi][mi] - dd[fd][mi] / 2.0, ds.domain_left_edge[i].d)
        assert_equal(dd[fd].max(), (ds.domain_width/ds.domain_dimensions)[i].d)
    assert_almost_equal(dd["cell_volume"].sum(dtype="float64"),
                        np.pi*ds.domain_width[0]**2 * ds.domain_width[1])
    assert_equal(dd["index", "path_element_r"], dd["index", "dr"])
    assert_equal(dd["index", "path_element_z"], dd["index", "dz"])
    assert_equal(dd["index", "path_element_theta"],
                 dd["index", "r"] * dd["index", "dtheta"])

def test_noise_plot_lin():
    ds = fake_amr_ds(geometry="cylindrical")
    add_noise_fields(ds)

    def create_image_lin(filename_prefix):
        fields = ['density'] + ['noise%d' % i for i in range(4)]
        p = SlicePlot(ds, 'z', fields)
        p.set_log('all', False) # <- the one diff line with test_noise_plot_log
        # lines bellow can later be replaced with the following snippet (see PR 2497)
        # 
        # fig = p.export_to_mpl_figure((2,2))
        # fig.savefig('noise_plot_lin.png')
        from mpl_toolkits.axes_grid1 import AxesGrid
        fig = plt.figure()

        grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols = (2, 2),
                axes_pad = 1.0,
                label_mode = "1",
                share_all = True,
                cbar_location="right",
                cbar_mode="each",
                cbar_size="3%",
                cbar_pad="0%")

        for i, field in enumerate(fields):
            plot = p.plots[field]
            plot.figure = fig
            plot.axes = grid[i].axes
            plot.cax = grid.cbar_axes[i]

        p._setup_plots()
        plt.savefig('noise_plot_lin.png')

    test = GenericImageTest(ds, create_image_lin, 12)
    test.prefix = "test_noise_plot_lin"
    test_noise_plot_lin.__name__ = test.description
    yield test

def test_noise_plot_log():
    ds = fake_amr_ds(geometry="cylindrical")
    add_noise_fields(ds)

    def create_image_log(filename_prefix):
        fields = ['density'] + ['noise%d' % i for i in range(4)]
        p = SlicePlot(ds, 'z', fields)
        # lines bellow can later be replaced with the following snippet (see PR 2497)
        # 
        # fig = p.export_to_mpl_figure((2,2))
        # fig.savefig('noise_plot_log.png')
        from mpl_toolkits.axes_grid1 import AxesGrid
        fig = plt.figure()

        grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols = (2, 2),
                axes_pad = 1.0,
                label_mode = "1",
                share_all = True,
                cbar_location="right",
                cbar_mode="each",
                cbar_size="3%",
                cbar_pad="0%")

        for i, field in enumerate(fields):
            plot = p.plots[field]
            plot.figure = fig
            plot.axes = grid[i].axes
            plot.cax = grid.cbar_axes[i]

        p._setup_plots()
        plt.savefig('noise_plot_log.png')

    test = GenericImageTest(ds, create_image_log, 12)
    test.prefix = "test_noise_plot_log"
    test_noise_plot_log.__name__ = test.description
    yield test