import numpy as np

from yt.loaders import load_uniform_grid
from yt.visualization.plot_window import ProjectionPlot


def test_ds_arr_invariance_under_projection_plot(tmp_path):
    data_array = np.random.random((10, 10, 10))
    bbox = np.array([[-100, 100], [-100, 100], [-100, 100]])
    data = {("gas", "density"): (data_array, "g*cm**(-3)")}
    ds = load_uniform_grid(data, data_array.shape, length_unit="kpc", bbox=bbox)

    start_source = np.array((0, 0, -0.5))
    end_source = np.array((0, 0, 0.5))
    start = ds.arr(start_source, "unitary")
    end = ds.arr(end_source, "unitary")

    start_i = start.copy()
    end_i = end.copy()

    p = ProjectionPlot(ds, 0, "number_density")
    p.annotate_line(start, end)
    p.save(tmp_path)

    # for lack of a unyt.testing.assert_unit_array_equal function
    np.testing.assert_array_equal(start_i, start)
    assert start_i.units == start.units
    np.testing.assert_array_equal(end_i, end)
    assert end_i.units == end.units
