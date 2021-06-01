from yt.testing import requires_file
from yt.utilities.answer_testing.framework import data_dir_load
from yt.visualization.plot_window import SlicePlot


@requires_file("bw_cartesian_3d/output0001.dat")
def test_code_units_xy_labels():
    ds = data_dir_load(
        "bw_cartesian_3d/output0001.dat", kwargs=dict(unit_system="code")
    )
    p = SlicePlot(ds, "x", ("gas", "density"))

    ax = p.plots[("gas", "density")].axes
    assert "code length" in ax.get_xlabel().replace("\\", "")
    assert "code length" in ax.get_ylabel().replace("\\", "")
