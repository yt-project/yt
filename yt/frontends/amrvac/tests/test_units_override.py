from yt.testing import assert_allclose_units, assert_raises, requires_file
from yt.units import YTQuantity
from yt.utilities.answer_testing.framework import data_dir_load

khi_cartesian_2D = "amrvac/kh_2d0000.dat"

# Tests for units: check that overriding certain units yields the correct derived units.
# The following are the correct normalisations
# based on length, numberdensity and temperature
length_unit = (1e9, "cm")
numberdensity_unit = (1e9, "cm**-3")
temperature_unit = (1e6, "K")
density_unit = (2.341670657200000e-15, "g*cm**-3")
mass_unit = (2.341670657200000e12, "g")
velocity_unit = (1.164508387441102e07, "cm*s**-1")
pressure_unit = (3.175492240000000e-01, "dyn*cm**-2")
time_unit = (8.587314705370271e01, "s")
magnetic_unit = (1.997608879907716, "gauss")


def _assert_normalisations_equal(ds):
    assert_allclose_units(ds.length_unit, YTQuantity(*length_unit))
    assert_allclose_units(ds.temperature_unit, YTQuantity(*temperature_unit))
    assert_allclose_units(ds.density_unit, YTQuantity(*density_unit))
    assert_allclose_units(ds.mass_unit, YTQuantity(*mass_unit))
    assert_allclose_units(ds.velocity_unit, YTQuantity(*velocity_unit))
    assert_allclose_units(ds.pressure_unit, YTQuantity(*pressure_unit))
    assert_allclose_units(ds.time_unit, YTQuantity(*time_unit))
    assert_allclose_units(ds.magnetic_unit, YTQuantity(*magnetic_unit))


@requires_file(khi_cartesian_2D)
def test_normalisations_length_temp_nb():
    # overriding length, temperature, numberdensity
    overrides = dict(
        length_unit=length_unit,
        temperature_unit=temperature_unit,
        numberdensity_unit=numberdensity_unit,
    )
    ds = data_dir_load(khi_cartesian_2D, kwargs={"units_override": overrides})
    _assert_normalisations_equal(ds)


@requires_file(khi_cartesian_2D)
def test_normalisations_length_temp_mass():
    # overriding length, temperature, mass
    overrides = dict(
        length_unit=length_unit, temperature_unit=temperature_unit, mass_unit=mass_unit
    )
    ds = data_dir_load(khi_cartesian_2D, kwargs={"units_override": overrides})
    _assert_normalisations_equal(ds)


@requires_file(khi_cartesian_2D)
def test_normalisations_length_time_mass():
    # overriding length, time, mass
    overrides = dict(length_unit=length_unit, time_unit=time_unit, mass_unit=mass_unit)
    ds = data_dir_load(khi_cartesian_2D, kwargs={"units_override": overrides})
    _assert_normalisations_equal(ds)


@requires_file(khi_cartesian_2D)
def test_normalisations_length_vel_nb():
    # overriding length, velocity, numberdensity
    overrides = dict(
        length_unit=length_unit,
        velocity_unit=velocity_unit,
        numberdensity_unit=numberdensity_unit,
    )
    ds = data_dir_load(khi_cartesian_2D, kwargs={"units_override": overrides})
    _assert_normalisations_equal(ds)


@requires_file(khi_cartesian_2D)
def test_normalisations_length_vel_mass():
    # overriding length, velocity, mass
    overrides = dict(
        length_unit=length_unit, velocity_unit=velocity_unit, mass_unit=mass_unit
    )
    ds = data_dir_load(khi_cartesian_2D, kwargs={"units_override": overrides})
    _assert_normalisations_equal(ds)


@requires_file(khi_cartesian_2D)
def test_normalisations_default():
    # test default normalisations, without overrides
    ds = data_dir_load(khi_cartesian_2D)
    assert_allclose_units(ds.length_unit, YTQuantity(1, "cm"))
    assert_allclose_units(ds.temperature_unit, YTQuantity(1, "K"))
    assert_allclose_units(
        ds.density_unit, YTQuantity(2.341670657200000e-24, "g*cm**-3")
    )
    assert_allclose_units(ds.mass_unit, YTQuantity(2.341670657200000e-24, "g"))
    assert_allclose_units(
        ds.velocity_unit, YTQuantity(1.164508387441102e04, "cm*s**-1")
    )
    assert_allclose_units(
        ds.pressure_unit, YTQuantity(3.175492240000000e-16, "dyn*cm**-2")
    )
    assert_allclose_units(ds.time_unit, YTQuantity(8.587314705370271e-05, "s"))
    assert_allclose_units(ds.magnetic_unit, YTQuantity(6.316993934686148e-08, "gauss"))


@requires_file(khi_cartesian_2D)
def test_normalisations_too_many_args():
    # test forbidden case: too many arguments (max 3 are allowed)
    overrides = dict(
        length_unit=length_unit,
        numberdensity_unit=numberdensity_unit,
        temperature_unit=temperature_unit,
        time_unit=time_unit,
    )
    with assert_raises(ValueError):
        data_dir_load(khi_cartesian_2D, kwargs={"units_override": overrides})


@requires_file(khi_cartesian_2D)
def test_normalisations_vel_and_length():
    # test forbidden case: both velocity and temperature are specified as overrides
    overrides = dict(
        length_unit=length_unit,
        velocity_unit=velocity_unit,
        temperature_unit=temperature_unit,
    )
    with assert_raises(ValueError):
        data_dir_load(khi_cartesian_2D, kwargs={"units_override": overrides})
