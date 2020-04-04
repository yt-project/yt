import matplotlib.pyplot as plt
from yt.testing import \
    fake_amr_ds, _geom_transforms, add_noise_fields
from yt.utilities.answer_testing.framework import \
    AxialPixelizationTest, \
    GenericImageTest
from yt import SlicePlot

def test_axial_pixelization():
    for geom in sorted(_geom_transforms):
        ds = fake_amr_ds(geometry=geom)
        yield AxialPixelizationTest(ds)


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