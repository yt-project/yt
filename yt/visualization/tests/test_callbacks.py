import contextlib
import shutil
import tempfile

from numpy.testing import assert_array_equal, assert_raises

import yt.units as u
from yt.config import ytcfg
from yt.loaders import load
from yt.testing import (
    assert_fname,
    fake_amr_ds,
    fake_hexahedral_ds,
    fake_random_ds,
    fake_tetrahedral_ds,
    requires_file,
)
from yt.utilities.exceptions import YTDataTypeUnsupported, YTPlotCallbackError
from yt.visualization.api import OffAxisSlicePlot, ProjectionPlot, SlicePlot
from yt.visualization.plot_container import accepts_all_fields

# These are a very simple set of tests that verify that each callback is or is
# not working.  They only check that it functions without an error; they do not
# check that it is providing correct results.

# These are the callbacks still to test:
#
#  X velocity
#  X magnetic_field
#  X quiver
#  X contour
#  X grids
#  X streamlines
#    units
#  X line
#    cquiver
#    clumps
#  X arrow
#  X marker
#  X sphere
#    hop_circles
#    hop_particles
#    coord_axes
#  X text
#  X particles
#    title
#    flash_ray_data
#  X timestamp
#  X scale
#    material_boundary
#  X ray
#  X line_integral_convolution

# cylindrical data for callback test
cyl_2d = "WDMerger_hdf5_chk_1000/WDMerger_hdf5_chk_1000.hdf5"
cyl_3d = "MHD_Cyl3d_hdf5_plt_cnt_0100/MHD_Cyl3d_hdf5_plt_cnt_0100.hdf5"


@contextlib.contextmanager
def _cleanup_fname():
    tmpdir = tempfile.mkdtemp()
    yield tmpdir
    shutil.rmtree(tmpdir)


def test_timestamp_callback():
    with _cleanup_fname() as prefix:
        ax = "z"
        vector = [1.0, 1.0, 1.0]
        ds = fake_amr_ds(fields=("density",), units=("g/cm**3",))
        p = ProjectionPlot(ds, ax, ("gas", "density"))
        p.annotate_timestamp()
        assert_fname(p.save(prefix)[0])
        p = SlicePlot(ds, ax, ("gas", "density"))
        p.annotate_timestamp()
        assert_fname(p.save(prefix)[0])
        p = OffAxisSlicePlot(ds, vector, ("gas", "density"))
        p.annotate_timestamp()
        assert_fname(p.save(prefix)[0])
        # Now we'll check a few additional minor things
        p = SlicePlot(ds, "x", ("gas", "density"))
        p.annotate_timestamp(corner="lower_right", redshift=True, draw_inset_box=True)
        p.save(prefix)

    with _cleanup_fname() as prefix:
        ds = fake_amr_ds(fields=("density",), units=("g/cm**3",), geometry="spherical")
        p = ProjectionPlot(ds, "r", ("gas", "density"))
        p.annotate_timestamp(coord_system="data")
        assert_raises(YTDataTypeUnsupported, p.save, prefix)
        p = ProjectionPlot(ds, "r", ("gas", "density"))
        p.annotate_timestamp(coord_system="axis")
        assert_fname(p.save(prefix)[0])


def test_timestamp_callback_code_units():
    # see https://github.com/yt-project/yt/issues/3869
    with _cleanup_fname() as prefix:
        ds = fake_random_ds(2, unit_system="code")
        p = SlicePlot(ds, "z", ("gas", "density"))
        p.annotate_timestamp()
        assert_fname(p.save(prefix)[0])


def test_scale_callback():
    with _cleanup_fname() as prefix:
        ax = "z"
        vector = [1.0, 1.0, 1.0]
        ds = fake_amr_ds(fields=("density",), units=("g/cm**3",))
        p = ProjectionPlot(ds, ax, ("gas", "density"))
        p.annotate_scale()
        assert_fname(p.save(prefix)[0])
        p = ProjectionPlot(ds, ax, ("gas", "density"), width=(0.5, 1.0))
        p.annotate_scale()
        assert_fname(p.save(prefix)[0])
        p = ProjectionPlot(ds, ax, ("gas", "density"), width=(1.0, 1.5))
        p.annotate_scale()
        assert_fname(p.save(prefix)[0])
        p = SlicePlot(ds, ax, ("gas", "density"))
        p.annotate_scale()
        assert_fname(p.save(prefix)[0])
        p = OffAxisSlicePlot(ds, vector, ("gas", "density"))
        p.annotate_scale()
        assert_fname(p.save(prefix)[0])
        # Now we'll check a few additional minor things
        p = SlicePlot(ds, "x", ("gas", "density"))
        p.annotate_scale(corner="upper_right", coeff=10.0, unit="kpc")
        assert_fname(p.save(prefix)[0])
        p = SlicePlot(ds, "x", ("gas", "density"))
        p.annotate_scale(text_args={"size": 24})
        assert_fname(p.save(prefix)[0])
        p = SlicePlot(ds, "x", ("gas", "density"))
        p.annotate_scale(text_args={"font": 24})
        assert_raises(YTPlotCallbackError)

    with _cleanup_fname() as prefix:
        ds = fake_amr_ds(fields=("density",), units=("g/cm**3",), geometry="spherical")
        p = ProjectionPlot(ds, "r", ("gas", "density"))
        p.annotate_scale()
        assert_raises(YTDataTypeUnsupported, p.save, prefix)
        p = ProjectionPlot(ds, "r", ("gas", "density"))
        p.annotate_scale(coord_system="axis")
        assert_raises(YTDataTypeUnsupported, p.save, prefix)


def test_line_callback():
    with _cleanup_fname() as prefix:
        ax = "z"
        vector = [1.0, 1.0, 1.0]
        ds = fake_amr_ds(fields=("density",), units=("g/cm**3",))
        p = ProjectionPlot(ds, ax, ("gas", "density"))
        p.annotate_line([0.1, 0.1, 0.1], [0.5, 0.5, 0.5])
        assert_fname(p.save(prefix)[0])
        p = SlicePlot(ds, ax, ("gas", "density"))
        p.annotate_line([0.1, 0.1, 0.1], [0.5, 0.5, 0.5])
        assert_fname(p.save(prefix)[0])
        p = OffAxisSlicePlot(ds, vector, ("gas", "density"))
        p.annotate_line([0.1, 0.1, 0.1], [0.5, 0.5, 0.5])
        assert_fname(p.save(prefix)[0])
        # Now we'll check a few additional minor things
        p = SlicePlot(ds, "x", ("gas", "density"))
        p.annotate_line(
            [0.1, 0.1], [0.5, 0.5], coord_system="axis", plot_args={"color": "red"}
        )
        p.save(prefix)

    with _cleanup_fname() as prefix:
        ds = fake_amr_ds(fields=("density",), units=("g/cm**3",), geometry="spherical")
        p = ProjectionPlot(ds, "r", ("gas", "density"))
        p.annotate_line([0.1, 0.1, 0.1], [0.5, 0.5, 0.5])
        assert_raises(YTDataTypeUnsupported, p.save, prefix)
        p = ProjectionPlot(ds, "r", ("gas", "density"))
        p.annotate_line([0.1, 0.1], [0.5, 0.5], coord_system="axis")
        assert_fname(p.save(prefix)[0])


def test_ray_callback():
    with _cleanup_fname() as prefix:
        ax = "z"
        vector = [1.0, 1.0, 1.0]
        ds = fake_amr_ds(fields=("density",), units=("g/cm**3",))
        ray = ds.ray((0.1, 0.2, 0.3), (0.6, 0.8, 0.5))
        oray = ds.ortho_ray(0, (0.3, 0.4))
        p = ProjectionPlot(ds, ax, ("gas", "density"))
        p.annotate_ray(oray)
        p.annotate_ray(ray)
        assert_fname(p.save(prefix)[0])
        p = SlicePlot(ds, ax, ("gas", "density"))
        p.annotate_ray(oray)
        p.annotate_ray(ray)
        assert_fname(p.save(prefix)[0])
        p = OffAxisSlicePlot(ds, vector, ("gas", "density"))
        p.annotate_ray(oray)
        p.annotate_ray(ray)
        assert_fname(p.save(prefix)[0])
        # Now we'll check a few additional minor things
        p = SlicePlot(ds, "x", ("gas", "density"))
        p.annotate_ray(oray)
        p.annotate_ray(ray, plot_args={"color": "red"})
        p.save(prefix)

    with _cleanup_fname() as prefix:
        ds = fake_amr_ds(fields=("density",), units=("g/cm**3",), geometry="spherical")
        ray = ds.ray((0.1, 0.2, 0.3), (0.6, 0.8, 0.5))
        oray = ds.ortho_ray(0, (0.3, 0.4))
        p = ProjectionPlot(ds, "r", ("gas", "density"))
        p.annotate_ray(oray)
        assert_raises(YTDataTypeUnsupported, p.save, prefix)
        p = ProjectionPlot(ds, "r", ("gas", "density"))
        p.annotate_ray(ray)
        assert_raises(YTDataTypeUnsupported, p.save, prefix)


def test_arrow_callback():
    with _cleanup_fname() as prefix:
        ax = "z"
        vector = [1.0, 1.0, 1.0]
        ds = fake_amr_ds(fields=("density",), units=("g/cm**3",))
        p = ProjectionPlot(ds, ax, ("gas", "density"))
        p.annotate_arrow([0.5, 0.5, 0.5])
        assert_fname(p.save(prefix)[0])
        p = SlicePlot(ds, ax, ("gas", "density"))
        p.annotate_arrow([0.5, 0.5, 0.5])
        assert_fname(p.save(prefix)[0])
        p = OffAxisSlicePlot(ds, vector, ("gas", "density"))
        p.annotate_arrow([0.5, 0.5, 0.5])
        assert_fname(p.save(prefix)[0])
        # Now we'll check a few additional minor things
        p = SlicePlot(ds, "x", ("gas", "density"))
        p.annotate_arrow([0.5, 0.5], coord_system="axis", length=0.05)
        p.annotate_arrow(
            [[0.5, 0.6], [0.5, 0.6], [0.5, 0.6]], coord_system="data", length=0.05
        )
        p.annotate_arrow(
            [[0.5, 0.6, 0.8], [0.5, 0.6, 0.8], [0.5, 0.6, 0.8]],
            coord_system="data",
            length=0.05,
        )
        p.annotate_arrow(
            [[0.5, 0.6, 0.8], [0.5, 0.6, 0.8]], coord_system="axis", length=0.05
        )
        p.annotate_arrow(
            [[0.5, 0.6, 0.8], [0.5, 0.6, 0.8]], coord_system="figure", length=0.05
        )
        p.annotate_arrow(
            [[0.5, 0.6, 0.8], [0.5, 0.6, 0.8]], coord_system="plot", length=0.05
        )
        p.save(prefix)

    with _cleanup_fname() as prefix:
        ds = fake_amr_ds(fields=("density",), units=("g/cm**3",), geometry="spherical")
        p = ProjectionPlot(ds, "r", ("gas", "density"))
        p.annotate_arrow([0.5, 0.5, 0.5])
        assert_raises(YTDataTypeUnsupported, p.save, prefix)
        p = ProjectionPlot(ds, "r", ("gas", "density"))
        p.annotate_arrow([0.5, 0.5], coord_system="axis")
        assert_fname(p.save(prefix)[0])


def test_marker_callback():
    with _cleanup_fname() as prefix:
        ax = "z"
        vector = [1.0, 1.0, 1.0]
        ds = fake_amr_ds(fields=("density",), units=("g/cm**3",))
        p = ProjectionPlot(ds, ax, ("gas", "density"))
        p.annotate_marker([0.5, 0.5, 0.5])
        assert_fname(p.save(prefix)[0])
        p = SlicePlot(ds, ax, ("gas", "density"))
        p.annotate_marker([0.5, 0.5, 0.5])
        assert_fname(p.save(prefix)[0])
        p = OffAxisSlicePlot(ds, vector, ("gas", "density"))
        p.annotate_marker([0.5, 0.5, 0.5])
        assert_fname(p.save(prefix)[0])
        # Now we'll check a few additional minor things
        p = SlicePlot(ds, "x", ("gas", "density"))
        coord = ds.arr([0.75, 0.75, 0.75], "unitary")
        coord.convert_to_units("kpc")
        p.annotate_marker(coord, coord_system="data")
        p.annotate_marker([0.5, 0.5], coord_system="axis", marker="*")
        p.annotate_marker([[0.5, 0.6], [0.5, 0.6], [0.5, 0.6]], coord_system="data")
        p.annotate_marker(
            [[0.5, 0.6, 0.8], [0.5, 0.6, 0.8], [0.5, 0.6, 0.8]], coord_system="data"
        )
        p.annotate_marker([[0.5, 0.6, 0.8], [0.5, 0.6, 0.8]], coord_system="axis")
        p.annotate_marker([[0.5, 0.6, 0.8], [0.5, 0.6, 0.8]], coord_system="figure")
        p.annotate_marker([[0.5, 0.6, 0.8], [0.5, 0.6, 0.8]], coord_system="plot")
        p.save(prefix)

    with _cleanup_fname() as prefix:
        ds = fake_amr_ds(fields=("density",), units=("g/cm**3",), geometry="spherical")
        p = ProjectionPlot(ds, "r", ("gas", "density"))
        p.annotate_marker([0.5, 0.5, 0.5])
        assert_raises(YTDataTypeUnsupported, p.save, prefix)
        p = ProjectionPlot(ds, "r", ("gas", "density"))
        p.annotate_marker([0.5, 0.5], coord_system="axis")
        assert_fname(p.save(prefix)[0])


def test_particles_callback():
    with _cleanup_fname() as prefix:
        ax = "z"
        ds = fake_amr_ds(fields=("density",), units=("g/cm**3",), particles=1)
        p = ProjectionPlot(ds, ax, ("gas", "density"))
        p.annotate_particles((10, "Mpc"))
        assert_fname(p.save(prefix)[0])
        p = SlicePlot(ds, ax, ("gas", "density"))
        p.annotate_particles((10, "Mpc"))
        assert_fname(p.save(prefix)[0])
        # Now we'll check a few additional minor things
        p = SlicePlot(ds, "x", ("gas", "density"))
        ad = ds.all_data()
        p.annotate_particles(
            (10, "Mpc"),
            p_size=1.0,
            col="k",
            marker="o",
            stride=1,
            ptype="all",
            alpha=1.0,
            data_source=ad,
        )
        p.save(prefix)

    with _cleanup_fname() as prefix:
        ds = fake_amr_ds(fields=("density",), units=("g/cm**3",), geometry="spherical")
        p = ProjectionPlot(ds, "r", ("gas", "density"))
        p.annotate_particles((10, "Mpc"))
        assert_raises(YTDataTypeUnsupported, p.save, prefix)


def test_sphere_callback():
    with _cleanup_fname() as prefix:
        ax = "z"
        vector = [1.0, 1.0, 1.0]
        ds = fake_amr_ds(fields=("density",), units=("g/cm**3",))
        p = ProjectionPlot(ds, ax, ("gas", "density"))
        p.annotate_sphere([0.5, 0.5, 0.5], 0.1)
        assert_fname(p.save(prefix)[0])
        p = SlicePlot(ds, ax, ("gas", "density"))
        p.annotate_sphere([0.5, 0.5, 0.5], 0.1)
        assert_fname(p.save(prefix)[0])
        p = OffAxisSlicePlot(ds, vector, ("gas", "density"))
        p.annotate_sphere([0.5, 0.5, 0.5], 0.1)
        assert_fname(p.save(prefix)[0])
        # Now we'll check a few additional minor things
        p = SlicePlot(ds, "x", ("gas", "density"))
        p.annotate_sphere([0.5, 0.5], 0.1, coord_system="axis", text="blah")
        p.save(prefix)

    with _cleanup_fname() as prefix:
        ds = fake_amr_ds(fields=("density",), units=("g/cm**3",), geometry="spherical")
        p = ProjectionPlot(ds, "r", ("gas", "density"))
        p.annotate_sphere([0.5, 0.5, 0.5], 0.1)
        assert_raises(YTDataTypeUnsupported, p.save, prefix)
        p = ProjectionPlot(ds, "r", ("gas", "density"))
        p.annotate_sphere([0.5, 0.5], 0.1, coord_system="axis", text="blah")
        assert_fname(p.save(prefix)[0])


def test_text_callback():
    with _cleanup_fname() as prefix:
        ax = "z"
        vector = [1.0, 1.0, 1.0]
        ds = fake_amr_ds(fields=("density",), units=("g/cm**3",))
        p = ProjectionPlot(ds, ax, ("gas", "density"))
        p.annotate_text([0.5, 0.5, 0.5], "dinosaurs!")
        assert_fname(p.save(prefix)[0])
        p = SlicePlot(ds, ax, ("gas", "density"))
        p.annotate_text([0.5, 0.5, 0.5], "dinosaurs!")
        assert_fname(p.save(prefix)[0])
        p = OffAxisSlicePlot(ds, vector, ("gas", "density"))
        p.annotate_text([0.5, 0.5, 0.5], "dinosaurs!")
        assert_fname(p.save(prefix)[0])
        # Now we'll check a few additional minor things
        p = SlicePlot(ds, "x", ("gas", "density"))
        p.annotate_text(
            [0.5, 0.5], "dinosaurs!", coord_system="axis", text_args={"color": "red"}
        )
        p.save(prefix)

    with _cleanup_fname() as prefix:
        ds = fake_amr_ds(fields=("density",), units=("g/cm**3",), geometry="spherical")
        p = ProjectionPlot(ds, "r", ("gas", "density"))
        p.annotate_text([0.5, 0.5, 0.5], "dinosaurs!")
        assert_raises(YTDataTypeUnsupported, p.save, prefix)
        p = ProjectionPlot(ds, "r", ("gas", "density"))
        p.annotate_text(
            [0.5, 0.5], "dinosaurs!", coord_system="axis", text_args={"color": "red"}
        )
        assert_fname(p.save(prefix)[0])


@requires_file(cyl_2d)
@requires_file(cyl_3d)
def test_velocity_callback():
    with _cleanup_fname() as prefix:
        ds = fake_amr_ds(
            fields=("density", "velocity_x", "velocity_y", "velocity_z"),
            units=("g/cm**3", "cm/s", "cm/s", "cm/s"),
        )
        for ax in "xyz":
            p = ProjectionPlot(
                ds, ax, ("gas", "density"), weight_field=("gas", "density")
            )
            p.annotate_velocity()
            assert_fname(p.save(prefix)[0])
            p = SlicePlot(ds, ax, ("gas", "density"))
            p.annotate_velocity()
            assert_fname(p.save(prefix)[0])
        # Test for OffAxis Slice
        p = SlicePlot(ds, [1, 1, 0], ("gas", "density"), north_vector=[0, 0, 1])
        p.annotate_velocity(factor=40, normalize=True)
        assert_fname(p.save(prefix)[0])
        # Now we'll check a few additional minor things
        p = SlicePlot(ds, "x", ("gas", "density"))
        p.annotate_velocity(factor=8, scale=0.5, scale_units="inches", normalize=True)
        assert_fname(p.save(prefix)[0])

    with _cleanup_fname() as prefix:
        ds = fake_hexahedral_ds(fields=[f"velocity_{ax}" for ax in "xyz"])
        sl = SlicePlot(ds, 1, ("connect1", "test"))
        sl.annotate_velocity()
        assert_fname(sl.save(prefix)[0])

    with _cleanup_fname() as prefix:
        ds = load(cyl_2d)
        slc = SlicePlot(ds, "theta", ("gas", "velocity_magnitude"))
        slc.annotate_velocity()
        assert_fname(slc.save(prefix)[0])

    with _cleanup_fname() as prefix:
        ds = load(cyl_3d)
        for ax in ["r", "z", "theta"]:
            slc = SlicePlot(ds, ax, ("gas", "velocity_magnitude"))
            slc.annotate_velocity()
            assert_fname(slc.save(prefix)[0])
            slc = ProjectionPlot(ds, ax, ("gas", "velocity_magnitude"))
            slc.annotate_velocity()
            assert_fname(slc.save(prefix)[0])

    with _cleanup_fname() as prefix:
        ds = fake_amr_ds(
            fields=("density", "velocity_r", "velocity_theta", "velocity_phi"),
            units=("g/cm**3", "cm/s", "cm/s", "cm/s"),
            geometry="spherical",
        )
        p = ProjectionPlot(ds, "r", ("gas", "density"))
        p.annotate_velocity(factor=40, normalize=True)
        assert_raises(YTDataTypeUnsupported, p.save, prefix)


@requires_file(cyl_2d)
@requires_file(cyl_3d)
def test_magnetic_callback():
    with _cleanup_fname() as prefix:
        ds = fake_amr_ds(
            fields=(
                "density",
                "magnetic_field_x",
                "magnetic_field_y",
                "magnetic_field_z",
            ),
            units=(
                "g/cm**3",
                "G",
                "G",
                "G",
            ),
        )
        for ax in "xyz":
            p = ProjectionPlot(
                ds, ax, ("gas", "density"), weight_field=("gas", "density")
            )
            p.annotate_magnetic_field()
            assert_fname(p.save(prefix)[0])
            p = SlicePlot(ds, ax, ("gas", "density"))
            p.annotate_magnetic_field()
            assert_fname(p.save(prefix)[0])
        # Test for OffAxis Slice
        p = SlicePlot(ds, [1, 1, 0], ("gas", "density"), north_vector=[0, 0, 1])
        p.annotate_magnetic_field(factor=40, normalize=True)
        assert_fname(p.save(prefix)[0])
        # Now we'll check a few additional minor things
        p = SlicePlot(ds, "x", ("gas", "density"))
        p.annotate_magnetic_field(
            factor=8, scale=0.5, scale_units="inches", normalize=True
        )
        assert_fname(p.save(prefix)[0])

    with _cleanup_fname() as prefix:
        ds = load(cyl_2d)
        slc = SlicePlot(ds, "theta", ("gas", "magnetic_field_strength"))
        slc.annotate_magnetic_field()
        assert_fname(slc.save(prefix)[0])

    with _cleanup_fname() as prefix:
        ds = load(cyl_3d)
        for ax in ["r", "z", "theta"]:
            slc = SlicePlot(ds, ax, ("gas", "magnetic_field_strength"))
            slc.annotate_magnetic_field()
            assert_fname(slc.save(prefix)[0])
            slc = ProjectionPlot(ds, ax, ("gas", "magnetic_field_strength"))
            slc.annotate_magnetic_field()
            assert_fname(slc.save(prefix)[0])

    with _cleanup_fname() as prefix:
        ds = fake_amr_ds(
            fields=(
                "density",
                "magnetic_field_r",
                "magnetic_field_theta",
                "magnetic_field_phi",
            ),
            units=(
                "g/cm**3",
                "G",
                "G",
                "G",
            ),
            geometry="spherical",
        )
        p = ProjectionPlot(ds, "r", ("gas", "density"))
        p.annotate_magnetic_field(
            factor=8, scale=0.5, scale_units="inches", normalize=True
        )
        assert_raises(YTDataTypeUnsupported, p.save, prefix)


@requires_file(cyl_2d)
@requires_file(cyl_3d)
def test_quiver_callback():
    with _cleanup_fname() as prefix:
        ds = fake_amr_ds(
            fields=("density", "velocity_x", "velocity_y", "velocity_z"),
            units=("g/cm**3", "cm/s", "cm/s", "cm/s"),
        )
        for ax in "xyz":
            p = ProjectionPlot(ds, ax, ("gas", "density"))
            p.annotate_quiver(("gas", "velocity_x"), ("gas", "velocity_y"))
            assert_fname(p.save(prefix)[0])
            p = ProjectionPlot(
                ds, ax, ("gas", "density"), weight_field=("gas", "density")
            )
            p.annotate_quiver(("gas", "velocity_x"), ("gas", "velocity_y"))
            assert_fname(p.save(prefix)[0])
            p = SlicePlot(ds, ax, ("gas", "density"))
            p.annotate_quiver(("gas", "velocity_x"), ("gas", "velocity_y"))
            assert_fname(p.save(prefix)[0])
        # Now we'll check a few additional minor things
        p = SlicePlot(ds, "x", ("gas", "density"))
        p.annotate_quiver(
            ("gas", "velocity_x"),
            ("gas", "velocity_y"),
            factor=8,
            scale=0.5,
            scale_units="inches",
            normalize=True,
            bv_x=0.5 * u.cm / u.s,
            bv_y=0.5 * u.cm / u.s,
        )
        assert_fname(p.save(prefix)[0])

    with _cleanup_fname() as prefix:
        ds = load(cyl_2d)
        slc = SlicePlot(ds, "theta", ("gas", "density"))
        slc.annotate_quiver(("gas", "velocity_r"), ("gas", "velocity_z"))
        assert_fname(slc.save(prefix)[0])

    with _cleanup_fname() as prefix:
        ds = load(cyl_3d)
        slc = SlicePlot(ds, "r", ("gas", "velocity_magnitude"))
        slc.annotate_quiver(("gas", "velocity_theta"), ("gas", "velocity_z"))
        assert_fname(slc.save(prefix)[0])
        slc = SlicePlot(ds, "z", ("gas", "velocity_magnitude"))
        slc.annotate_quiver(
            ("gas", "velocity_cartesian_x"), ("gas", "velocity_cartesian_y")
        )
        assert_fname(slc.save(prefix)[0])
        slc = SlicePlot(ds, "theta", ("gas", "velocity_magnitude"))
        slc.annotate_quiver(("gas", "velocity_r"), ("gas", "velocity_z"))
        assert_fname(slc.save(prefix)[0])

    with _cleanup_fname() as prefix:
        ds = fake_amr_ds(
            fields=("density", "velocity_x", "velocity_theta", "velocity_phi"),
            units=("g/cm**3", "cm/s", "cm/s", "cm/s"),
            geometry="spherical",
        )
        p = ProjectionPlot(ds, "r", ("gas", "density"))
        p.annotate_quiver(
            ("gas", "velocity_theta"),
            ("gas", "velocity_phi"),
            factor=8,
            scale=0.5,
            scale_units="inches",
            normalize=True,
            bv_x=0.5 * u.cm / u.s,
            bv_y=0.5 * u.cm / u.s,
        )
        assert_raises(YTDataTypeUnsupported, p.save, prefix)


@requires_file(cyl_2d)
def test_contour_callback():
    with _cleanup_fname() as prefix:
        ds = fake_amr_ds(fields=("density", "temperature"), units=("g/cm**3", "K"))
        for ax in "xyz":
            p = ProjectionPlot(ds, ax, ("gas", "density"))
            p.annotate_contour(("gas", "temperature"))
            assert_fname(p.save(prefix)[0])
            p = ProjectionPlot(
                ds, ax, ("gas", "density"), weight_field=("gas", "density")
            )
            p.annotate_contour(("gas", "temperature"))
            assert_fname(p.save(prefix)[0])
            p = SlicePlot(ds, ax, ("gas", "density"))
            p.annotate_contour(("gas", "temperature"))  # BREAKS WITH ndarray
            assert_fname(p.save(prefix)[0])
        # Now we'll check a few additional minor things
        p = SlicePlot(ds, "x", ("gas", "density"))
        p.annotate_contour(
            ("gas", "temperature"),
            ncont=10,
            factor=8,
            take_log=False,
            clim=(0.4, 0.6),
            plot_args={"linewidths": 2.0},
            label=True,
            text_args={"fontsize": "x-large"},
        )
        p.save(prefix)

        p = SlicePlot(ds, "x", ("gas", "density"))
        s2 = ds.slice(0, 0.2)
        p.annotate_contour(
            ("gas", "temperature"),
            ncont=10,
            factor=8,
            take_log=False,
            clim=(0.4, 0.6),
            plot_args={"linewidths": 2.0},
            label=True,
            text_args={"fontsize": "x-large"},
            data_source=s2,
        )
        p.save(prefix)

    with _cleanup_fname() as prefix:
        ds = load(cyl_2d)
        slc = SlicePlot(ds, "theta", ("gas", "plasma_beta"))
        slc.annotate_contour(
            ("gas", "plasma_beta"),
            ncont=2,
            factor=7.0,
            take_log=False,
            clim=(1.0e-1, 1.0e1),
            label=True,
            plot_args={"colors": ("c", "w"), "linewidths": 1},
            text_args={"fmt": "%1.1f"},
        )
        assert_fname(slc.save(prefix)[0])

    with _cleanup_fname() as prefix:
        ds = fake_amr_ds(
            fields=("density", "temperature"),
            units=("g/cm**3", "K"),
            geometry="spherical",
        )
        p = SlicePlot(ds, "r", ("gas", "density"))
        p.annotate_contour(
            ("gas", "temperature"),
            ncont=10,
            factor=8,
            take_log=False,
            clim=(0.4, 0.6),
            plot_args={"linewidths": 2.0},
            label=True,
            text_args={"fontsize": "x-large"},
        )
        assert_raises(YTDataTypeUnsupported, p.save, prefix)


@requires_file(cyl_2d)
def test_grids_callback():
    with _cleanup_fname() as prefix:
        ds = fake_amr_ds(fields=("density",), units=("g/cm**3",))
        for ax in "xyz":
            p = ProjectionPlot(ds, ax, ("gas", "density"))
            p.annotate_grids()
            assert_fname(p.save(prefix)[0])
            p = ProjectionPlot(
                ds, ax, ("gas", "density"), weight_field=("gas", "density")
            )
            p.annotate_grids()
            assert_fname(p.save(prefix)[0])
            p = SlicePlot(ds, ax, ("gas", "density"))
            p.annotate_grids()
            assert_fname(p.save(prefix)[0])
        # Now we'll check a few additional minor things
        p = SlicePlot(ds, "x", ("gas", "density"))
        p.annotate_grids(
            alpha=0.7,
            min_pix=10,
            min_pix_ids=30,
            draw_ids=True,
            id_loc="upper right",
            periodic=False,
            min_level=2,
            max_level=3,
            cmap="gist_stern",
        )
        p.save(prefix)

    with _cleanup_fname() as prefix:
        ds = load(cyl_2d)
        slc = SlicePlot(ds, "theta", ("gas", "density"))
        slc.annotate_grids()
        assert_fname(slc.save(prefix)[0])

    with _cleanup_fname() as prefix:
        ds = fake_amr_ds(fields=("density",), units=("g/cm**3",), geometry="spherical")
        p = SlicePlot(ds, "r", ("gas", "density"))
        p.annotate_grids(
            alpha=0.7,
            min_pix=10,
            min_pix_ids=30,
            draw_ids=True,
            id_loc="upper right",
            periodic=False,
            min_level=2,
            max_level=3,
            cmap="gist_stern",
        )
        assert_raises(YTDataTypeUnsupported, p.save, prefix)


@requires_file(cyl_2d)
def test_cell_edges_callback():
    with _cleanup_fname() as prefix:
        ds = fake_amr_ds(fields=("density",), units=("g/cm**3",))
        for ax in "xyz":
            p = ProjectionPlot(ds, ax, ("gas", "density"))
            p.annotate_cell_edges()
            assert_fname(p.save(prefix)[0])
            p = ProjectionPlot(
                ds, ax, ("gas", "density"), weight_field=("gas", "density")
            )
            p.annotate_cell_edges()
            assert_fname(p.save(prefix)[0])
            p = SlicePlot(ds, ax, ("gas", "density"))
            p.annotate_cell_edges()
            assert_fname(p.save(prefix)[0])
        # Now we'll check a few additional minor things
        p = SlicePlot(ds, "x", ("gas", "density"))
        p.annotate_cell_edges(alpha=0.7, line_width=0.9, color=(0.0, 1.0, 1.0))
        p.save(prefix)

    with _cleanup_fname() as prefix:
        ds = load(cyl_2d)
        slc = SlicePlot(ds, "theta", ("gas", "density"))
        slc.annotate_cell_edges()
        assert_fname(slc.save(prefix)[0])

    with _cleanup_fname() as prefix:
        ds = fake_amr_ds(fields=("density",), units=("g/cm**3",), geometry="spherical")
        p = SlicePlot(ds, "r", ("gas", "density"))
        p.annotate_cell_edges()
        assert_raises(YTDataTypeUnsupported, p.save, prefix)


def test_mesh_lines_callback():
    with _cleanup_fname() as prefix:

        ds = fake_hexahedral_ds()
        for field in ds.field_list:
            sl = SlicePlot(ds, 1, field)
            sl.annotate_mesh_lines(plot_args={"color": "black"})
            assert_fname(sl.save(prefix)[0])

        ds = fake_tetrahedral_ds()
        for field in ds.field_list:
            sl = SlicePlot(ds, 1, field)
            sl.annotate_mesh_lines(plot_args={"color": "black"})
            assert_fname(sl.save(prefix)[0])


@requires_file(cyl_2d)
@requires_file(cyl_3d)
def test_streamline_callback():

    with _cleanup_fname() as prefix:

        ds = fake_amr_ds(
            fields=("density", "velocity_x", "velocity_y", "magvel"),
            units=("g/cm**3", "cm/s", "cm/s", "cm/s"),
        )

        for ax in "xyz":
            # Projection plot tests
            p = ProjectionPlot(ds, ax, ("gas", "density"))
            p.annotate_streamlines(("gas", "velocity_x"), ("gas", "velocity_y"))
            assert_fname(p.save(prefix)[0])

            p = ProjectionPlot(
                ds, ax, ("gas", "density"), weight_field=("gas", "density")
            )
            p.annotate_streamlines(("gas", "velocity_x"), ("gas", "velocity_y"))
            assert_fname(p.save(prefix)[0])

            # Slice plot test
            p = SlicePlot(ds, ax, ("gas", "density"))
            p.annotate_streamlines(("gas", "velocity_x"), ("gas", "velocity_y"))
            assert_fname(p.save(prefix)[0])

            # Additional features
            p = SlicePlot(ds, ax, ("gas", "density"))
            p.annotate_streamlines(
                ("gas", "velocity_x"), ("gas", "velocity_y"), factor=32, density=4
            )
            assert_fname(p.save(prefix)[0])

            p = SlicePlot(ds, ax, ("gas", "density"))
            p.annotate_streamlines(
                ("gas", "velocity_x"),
                ("gas", "velocity_y"),
                field_color=("stream", "magvel"),
            )
            assert_fname(p.save(prefix)[0])

            p = SlicePlot(ds, ax, ("gas", "density"))
            p.annotate_streamlines(
                ("gas", "velocity_x"),
                ("gas", "velocity_y"),
                field_color=("stream", "magvel"),
                display_threshold=0.5,
                plot_args={
                    "cmap": ytcfg.get("yt", "default_colormap"),
                    "arrowstyle": "->",
                },
            )
            assert_fname(p.save(prefix)[0])

    # Axisymmetric dataset

    with _cleanup_fname() as prefix:

        ds = load(cyl_2d)
        slc = SlicePlot(ds, "theta", ("gas", "velocity_magnitude"))
        slc.annotate_streamlines(("gas", "velocity_r"), ("gas", "velocity_z"))
        assert_fname(slc.save(prefix)[0])

    with _cleanup_fname() as prefix:
        ds = load(cyl_3d)
        slc = SlicePlot(ds, "r", ("gas", "velocity_magnitude"))
        slc.annotate_streamlines(("gas", "velocity_theta"), ("gas", "velocity_z"))
        assert_fname(slc.save(prefix)[0])
        slc = SlicePlot(ds, "z", ("gas", "velocity_magnitude"))
        slc.annotate_streamlines(
            ("gas", "velocity_cartesian_x"), ("gas", "velocity_cartesian_y")
        )
        assert_fname(slc.save(prefix)[0])
        slc = SlicePlot(ds, "theta", ("gas", "velocity_magnitude"))
        slc.annotate_streamlines(("gas", "velocity_r"), ("gas", "velocity_z"))
        assert_fname(slc.save(prefix)[0])

    # Spherical dataset
    with _cleanup_fname() as prefix:

        ds = fake_amr_ds(
            fields=("density", "velocity_r", "velocity_theta", "velocity_phi"),
            units=("g/cm**3", "cm/s", "cm/s", "cm/s"),
            geometry="spherical",
        )
        p = SlicePlot(ds, "r", ("gas", "density"))
        p.annotate_streamlines(("gas", "velocity_theta"), ("gas", "velocity_phi"))
        assert_raises(YTDataTypeUnsupported, p.save, prefix)


@requires_file(cyl_2d)
@requires_file(cyl_3d)
def test_line_integral_convolution_callback():
    with _cleanup_fname() as prefix:
        ds = fake_amr_ds(
            fields=("density", "velocity_x", "velocity_y", "velocity_z"),
            units=("g/cm**3", "cm/s", "cm/s", "cm/s"),
        )
        for ax in "xyz":
            p = ProjectionPlot(ds, ax, ("gas", "density"))
            p.annotate_line_integral_convolution(
                ("gas", "velocity_x"), ("gas", "velocity_y")
            )
            assert_fname(p.save(prefix)[0])
            p = ProjectionPlot(
                ds, ax, ("gas", "density"), weight_field=("gas", "density")
            )
            p.annotate_line_integral_convolution(
                ("gas", "velocity_x"), ("gas", "velocity_y")
            )
            assert_fname(p.save(prefix)[0])
            p = SlicePlot(ds, ax, ("gas", "density"))
            p.annotate_line_integral_convolution(
                ("gas", "velocity_x"), ("gas", "velocity_y")
            )
            assert_fname(p.save(prefix)[0])
        # Now we'll check a few additional minor things
        p = SlicePlot(ds, "x", ("gas", "density"))
        p.annotate_line_integral_convolution(
            ("gas", "velocity_x"),
            ("gas", "velocity_y"),
            kernellen=100.0,
            lim=(0.4, 0.7),
            cmap=ytcfg.get("yt", "default_colormap"),
            alpha=0.9,
            const_alpha=True,
        )
        p.save(prefix)

    with _cleanup_fname() as prefix:
        ds = load(cyl_2d)
        slc = SlicePlot(ds, "theta", ("gas", "magnetic_field_strength"))
        slc.annotate_line_integral_convolution(
            ("gas", "magnetic_field_r"), ("gas", "magnetic_field_z")
        )
        assert_fname(slc.save(prefix)[0])

    with _cleanup_fname() as prefix:
        ds = load(cyl_3d)
        slc = SlicePlot(ds, "r", ("gas", "magnetic_field_strength"))
        slc.annotate_line_integral_convolution(
            ("gas", "magnetic_field_theta"), ("gas", "magnetic_field_z")
        )
        assert_fname(slc.save(prefix)[0])
        slc = SlicePlot(ds, "z", ("gas", "magnetic_field_strength"))
        slc.annotate_line_integral_convolution(
            ("gas", "magnetic_field_cartesian_x"), ("gas", "magnetic_field_cartesian_y")
        )
        assert_fname(slc.save(prefix)[0])
        slc = SlicePlot(ds, "theta", ("gas", "magnetic_field_strength"))
        slc.annotate_line_integral_convolution(
            ("gas", "magnetic_field_r"), ("gas", "magnetic_field_z")
        )
        assert_fname(slc.save(prefix)[0])

    with _cleanup_fname() as prefix:
        ds = fake_amr_ds(
            fields=("density", "velocity_r", "velocity_theta", "velocity_phi"),
            units=("g/cm**3", "cm/s", "cm/s", "cm/s"),
            geometry="spherical",
        )
        p = SlicePlot(ds, "r", ("gas", "density"))
        p.annotate_line_integral_convolution(
            ("gas", "velocity_theta"), ("gas", "velocity_phi")
        )
        assert_raises(YTDataTypeUnsupported, p.save, prefix)


def test_accepts_all_fields_decorator():
    fields = [
        ("gas", "density"),
        ("gas", "velocity_x"),
        ("gas", "pressure"),
        ("gas", "temperature"),
    ]
    units = ["g/cm**3", "cm/s", "dyn/cm**2", "K"]
    ds = fake_random_ds(16, fields=fields, units=units)
    plot = SlicePlot(ds, "z", fields=fields)

    # mocking a class method
    plot.fake_attr = {f: "not set" for f in fields}

    @accepts_all_fields
    def set_fake_field_attribute(self, field, value):
        self.fake_attr[field] = value
        return self

    # test on a single field
    plot = set_fake_field_attribute(plot, field=fields[0], value=1)
    assert_array_equal([plot.fake_attr[f] for f in fields], [1] + ["not set"] * 3)

    # test using "all" as a field
    plot = set_fake_field_attribute(plot, field="all", value=2)
    assert_array_equal(list(plot.fake_attr.values()), [2] * 4)
