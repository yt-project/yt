import matplotlib
import pytest
from matplotlib.colors import LogNorm, Normalize, SymLogNorm
from packaging.version import Version

import yt
from yt.testing import add_noise_fields, fake_particle_ds, fake_random_ds

MPL_VERSION = Version(matplotlib.__version__)


norms_to_test = [
    (Normalize(), "linear"),
    (LogNorm(), "log"),
]

if MPL_VERSION >= Version("3.4"):
    # don't test this in older versions of MPL because of a deprecation warning
    # that cannot be easily filtered
    norms_to_test.append((SymLogNorm(linthresh=0.01, vmin=-1, vmax=1), "symlog"))


@pytest.fixture(scope="session")
def myds():
    ds = fake_random_ds(16)
    add_noise_fields(ds)
    return ds


def test_lineplot_set_units(myds):
    ds = myds
    field = ("gas", "density")
    plot = yt.LinePlot(
        ds, field, start_point=[0, 0, 0], end_point=[1, 1, 1], npoints=512
    )
    plot.annotate_legend(field)
    plot.save("/tmp/tests/lineplot/defaults.png")

    plot.set_x_unit("cm")
    plot.save("/tmp/tests/lineplot/xunit.png")

    plot.set_unit(field, "kg/cm**3")
    plot.save("/tmp/tests/lineplot/xunit_zunit.png")

    plot.set_log(field, False)
    plot.save("/tmp/tests/lineplot/xunit_zunit_lin.png")


def test_profileplot_set_axis_properties(myds):
    ds = myds
    my_galaxy = ds.disk(ds.domain_center, [0.0, 0.0, 1.0], (10, "m"), (1, "m"))
    plot = yt.ProfilePlot(my_galaxy, ("gas", "density"), [("gas", "velocity_x")])
    plot.save("/tmp/tests/profileplot_defaults.png")

    plot.set_unit(("gas", "density"), "kg/cm**3")
    plot.save("/tmp/tests/profileplot_xunit.png")

    plot.set_log(("gas", "density"), False)
    plot.save("/tmp/tests/profileplot_xunit_xlin.png")

    plot.set_unit(("gas", "velocity_x"), "mile/hour")
    plot.save("/tmp/tests/profileplot_xunit_xlin_yunit.png")


@pytest.mark.parametrize("norm, name", norms_to_test)
def test_sliceplot_custom_norm(norm, name, myds):
    field = ("gas", "density")
    if name == "log" and field == "noise2":
        pytest.skip(reason="Can't use a log norm on a field with only negative values")
    p = yt.SlicePlot(myds, "z", field)
    p.set_norm(field, norm=norm)
    p.save(f"/tmp/tests/sliceplot/{name}")


@pytest.mark.skipif(
    MPL_VERSION < Version("3.2"), reason="TwoSlopeNorm is new in matplotlib 3.2"
)
@pytest.mark.parametrize("field", ("noise0", "noise1", "noise2", "noise3"))
def test_two_slope_norm(field, myds):
    from matplotlib.colors import TwoSlopeNorm

    p = yt.SlicePlot(myds, "z", field)
    p.set_norm(field, norm=TwoSlopeNorm(vcenter=0, vmin=-0.5, vmax=1))
    p.set_cmap(field, "RdBu")  # use a diverging colormap to help read the image
    p.save("/tmp/tests/sliceplot/norm/twoslope")


@pytest.mark.parametrize(
    "field", (("gas", "density"), "noise0", "noise1", "noise2", "noise3")
)
def test_sliceplot_set_log(field, myds):
    p = yt.SlicePlot(myds, "z", field)
    p.set_log(field, False)
    p.save("/tmp/tests/sliceplot/set_log")


def test_sliceplot_set_zlim_and_unit(myds):
    field = ("gas", "density")
    p = yt.SlicePlot(myds, "z", field)
    p.set_zlim(field, zmin=0)
    p.set_unit(field, "kg/m**3")
    p.save("/tmp/tests/sliceplot/set_zlim_and_unit.png")


def test_sliceplot_set_unit_and_zlim(myds):
    field = ("gas", "density")
    p = yt.SlicePlot(myds, "z", field)
    p.set_unit(field, "kg/m**3")
    p.set_zlim(field, zmin=0)
    p.save("/tmp/tests/sliceplot/set_unit_and_zlim.png")


@pytest.mark.parametrize("color", ("C0", None))
def test_sliceplot_set_background_color(myds, color):
    # see https://github.com/yt-project/yt/issues/3854
    normal = "z"
    p = yt.SlicePlot(myds, normal, ("gas", "density"), width=1.5)
    p.set_background_color(("gas", "density"), color="C0")
    p.save(f"/tmp/tests/sliceplot/background/{color}_log.png")
    p.set_log(("gas", "density"), False)
    p.save(f"/tmp/tests/sliceplot/background/{color}_lin.png")


@pytest.mark.parametrize(
    "field", (("gas", "density"), "noise0", "noise1", "noise2", "noise3")
)
def test_projectionplot_set_log(field, myds):
    p = yt.ProjectionPlot(myds, "z", field)
    p.set_log(field, False)
    p.save("/tmp/tests/projectionplot/set_log.png")


def test_projectionplot_set_zlim_and_unit(myds):
    field = ("gas", "density")
    p = yt.ProjectionPlot(myds, "z", field)
    p.set_zlim(field, zmin=0)
    p.set_unit(field, "kg/m**2")
    p.save("/tmp/tests/projectionplot/set_zlim_and_unit.png")


def test_projectionplot_set_unit_and_zlim(myds):
    field = ("gas", "density")
    p = yt.ProjectionPlot(myds, "z", field)
    p.set_unit(field, "kg/m**2")
    p.set_zlim(field, zmin=0)
    p.save("/tmp/tests/projectionplot/set_unit_and_zlim.png")


def test_phaseplot_set_colorbar_properties_implicit(myds):
    # see https://github.com/yt-project/yt/issues/2538
    my_sphere = myds.sphere("c", 1)
    p = yt.PhasePlot(
        my_sphere,
        ("gas", "noise1"),
        ("gas", "noise3"),
        [("gas", "density")],
        weight_field=None,
    )
    p.save("/tmp/tests/phaseplot/raw.png")  # TODO: drop this
    # using implicit units
    p.set_zlim(("gas", "density"), zmax=10)
    p.save("/tmp/tests/phaseplot/1_set_zlim.png")

    # changing units should affect the colorbar and not the image
    p.set_unit(("gas", "density"), "kg/AU**3")
    p.save("/tmp/tests/phaseplot/1_set_zlim_set_unit.png")


def test_phaseplot_set_colorbar_properties_explicit(myds):
    # see https://github.com/yt-project/yt/issues/2538
    my_sphere = myds.sphere("c", 1)
    p = yt.PhasePlot(
        my_sphere,
        ("gas", "noise1"),
        ("gas", "noise3"),
        [("gas", "density")],
        weight_field=None,
    )
    # using explit units, we expect the colorbar units to stay unchanged
    p.set_zlim(("gas", "density"), zmin=(1e36, "kg/AU**3"))
    p.save("/tmp/tests/phaseplot/2_set_zlim.png")

    # ... until we set them explicitly
    p.set_unit(("gas", "density"), "kg/AU**3")
    p.save("/tmp/tests/phaseplot/2_set_zlim_set_unit.png")


def test_particleprojectionplot_colorbar_properties():
    ds = fake_particle_ds(npart=100)
    field = ("all", "particle_mass")
    p = yt.ParticleProjectionPlot(ds, 2, field)
    p.set_buff_size(10)
    p.save("/tmp/tests/particleprojectionplot/raw.png")

    p.set_unit(field, "Msun")
    p.save("/tmp/tests/particleprojectionplot/set_unit.png")

    p.set_zlim(field, zmax=1e-35)
    p.save("/tmp/tests/particleprojectionplot/set_unit_set_zlim.png")

    p.set_log(field, False)
    p.save("/tmp/tests/particleprojectionplot/set_unit_set_zlim_set_log.png")


# @attr(ANSWER_TEST_TAG)
# def test_TEMPLATE():
#     def create_image(filename_prefix):
#         ...
#     test = GenericImageTest(ds, create_image, 12)
#     test.prefix = "test_TEMPLATE"
#     test_TEMPLATE.__name__ = test.description
#     yield test
