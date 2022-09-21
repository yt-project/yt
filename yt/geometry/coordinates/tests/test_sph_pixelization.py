import yt
from yt.testing import assert_rel_equal, requires_file
from yt.utilities.math_utils import compute_stddev_image

magneticum = "MagneticumCluster/snap_132"

mag_kwargs = dict(
    long_ids=True,
    field_spec="magneticum_box2_hr",
)


@requires_file(magneticum)
def test_sph_moment():
    ds = yt.load(magneticum, **mag_kwargs)

    def _vysq(field, data):
        return data["gas", "velocity_y"] ** 2

    ds.add_field(("gas", "vysq"), _vysq, sampling_type="local", units="cm**2/s**2")
    prj1 = yt.ProjectionPlot(
        ds,
        "y",
        [("gas", "velocity_y"), ("gas", "vysq")],
        weight_field=("gas", "density"),
        moment=1,
        buff_size=(400, 400),
    )
    prj2 = yt.ProjectionPlot(
        ds,
        "y",
        ("gas", "velocity_y"),
        moment=2,
        weight_field=("gas", "density"),
        buff_size=(400, 400),
    )
    sigy = compute_stddev_image(prj1.frb["gas", "vysq"], prj1.frb["gas", "velocity_y"])
    assert_rel_equal(sigy, prj2.frb["gas", "velocity_y"].d, 10)
