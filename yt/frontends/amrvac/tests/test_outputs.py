import numpy as np # NOQA
import pytest

import yt # NOQA
from yt.testing import assert_allclose_units, requires_file
from yt.frontends.amrvac.api import AMRVACDataset, AMRVACGrid
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils


blastwave_spherical_2D = "amrvac/bw_2d0000.dat"
khi_cartesian_2D = "amrvac/kh_2d0000.dat"
khi_cartesian_3D = "amrvac/kh_3D0000.dat"
jet_cylindrical_25D = "amrvac/Jet0003.dat"
riemann_cartesian_175D = "amrvac/R_1d0005.dat"
blastwave_cartesian_3D = "amrvac/bw_3d0000.dat"
blastwave_polar_2D = "amrvac/bw_polar_2D0000.dat"
blastwave_cylindrical_3D = "amrvac/bw_cylindrical_3D0000.dat"


@pytest.mark.answer_test
@pytest.mark.usefixtures('answer_file')
class TestAMRVac(fw.AnswerTest):
    @requires_file(khi_cartesian_2D)
    def test_AMRVACDataset(self):
        assert isinstance(utils.data_dir_load(khi_cartesian_2D), AMRVACDataset)

    @utils.requires_ds(blastwave_cartesian_3D)
    def test_domain_size(self):
        #"Check for correct box size, see bw_3d.par"
        ds = utils.data_dir_load(blastwave_cartesian_3D)
        for lb in ds.domain_left_edge:
            assert int(lb) == 0
        for rb in ds.domain_right_edge:
            assert int(rb) == 2
        for w in ds.domain_width:
            assert int(w) == 2

    @requires_file(blastwave_cartesian_3D)
    def test_grid_attributes(self):
        #"Check various grid attributes"
        ds = utils.data_dir_load(blastwave_cartesian_3D)
        grids = ds.index.grids
        assert ds.index.max_level == 2
        for g in grids:
            assert isinstance(g, AMRVACGrid)
            assert isinstance(g.LeftEdge, yt.units.yt_array.YTArray)
            assert isinstance(g.RightEdge, yt.units.yt_array.YTArray)
            assert isinstance(g.ActiveDimensions, np.ndarray)
            assert isinstance(g.Level, (np.int32, np.int64, int))

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(blastwave_polar_2D)
    def test_bw_polar_2d(self, a, d, w, f, ds_bw_polar_2D):
        self.hashes.update(self.small_patch_amr(ds_bw_polar_2D, f, w, a, d))

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(blastwave_cartesian_3D)
    def test_blastwave_cartesian_3D(self, a, d, w, f, ds_blastwave_cartesian_3D):
        self.hashes.update(self.small_patch_amr(ds_blastwave_cartesian_3D, f, w, a, d))

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(blastwave_spherical_2D)
    def test_blastwave_spherical_2D(self, a, d, w, f, ds_blastwave_spherical_2D):
        self.hashes.update(self.small_patch_amr(ds_blastwave_spherical_2D, f, w, a, d))

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(blastwave_cylindrical_3D)
    def test_blastwave_cylindrical_3D(self, a, d, w, f, ds_blastwave_cylindrical_3D):
        self.hashes.update(self.small_patch_amr(ds_blastwave_cylindrical_3D, f, w, a, d))

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(khi_cartesian_2D)
    def test_khi_cartesian_2D(self, a, d, w, f, ds_khi_cartesian_2D):
        self.hashes.update(self.small_patch_amr(ds_khi_cartesian_2D, f, w, a, d))

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(khi_cartesian_3D)
    def test_khi_cartesian_3D(self, a, d, w, f, ds_khi_cartesian_3D):
        self.hashes.update(self.small_patch_amr(ds_khi_cartesian_3D, f, w, a, d))

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(jet_cylindrical_25D)
    def test_jet_cylindrical_25D(self, a, d, w, f, ds_jet_cylindrical_25D):
        self.hashes.update(self.small_patch_amr(ds_jet_cylindrical_25D, f, w, a, d))

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(riemann_cartesian_175D)
    def test_riemann_cartesian_175D(self, a, d, w, f, ds_riemann_cartesian_175D):
        self.hashes.update(self.small_patch_amr(ds_riemann_cartesian_175D, f, w, a, d))


# Tests for units: verify that overriding certain units yields the correct derived units.
# The following are correct normalisations based on length, numberdensity and temperature
length_unit = (1e9, 'cm')
numberdensity_unit = (1e9, 'cm**-3')
temperature_unit = (1e6, 'K')
density_unit = (2.341670657200000e-15, 'g*cm**-3')
mass_unit = (2.341670657200000e+12, 'g')
velocity_unit = (1.164508387441102e+07, 'cm*s**-1')
pressure_unit = (3.175492240000000e-01, 'dyn*cm**-2')
time_unit = (8.587314705370271e+01, 's')
magnetic_unit = (1.997608879907716, 'gauss')

def _assert_normalisations_equal(ds):
    assert_allclose_units(ds.length_unit, YTQuantity(*length_unit))
    assert_allclose_units(ds.numberdensity_unit, YTQuantity(*numberdensity_unit))
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
    overrides = dict(length_unit=length_unit, temperature_unit=temperature_unit,
                     numberdensity_unit=numberdensity_unit)
    ds = data_dir_load(khi_cartesian_2D, kwargs={'units_override': overrides})
    _assert_normalisations_equal(ds)

@requires_file(khi_cartesian_2D)
def test_normalisations_length_temp_mass():
    # overriding length, temperature, mass
    overrides = dict(length_unit=length_unit, temperature_unit=temperature_unit,
                     mass_unit=mass_unit)
    ds = data_dir_load(khi_cartesian_2D, kwargs={'units_override': overrides})
    _assert_normalisations_equal(ds)

@requires_file(khi_cartesian_2D)
def test_normalisations_length_time_mass():
    # overriding length, time, mass
    overrides = dict(length_unit=length_unit, time_unit=time_unit, mass_unit=mass_unit)
    ds = data_dir_load(khi_cartesian_2D, kwargs={'units_override': overrides})
    _assert_normalisations_equal(ds)

@requires_file(khi_cartesian_2D)
def test_normalisations_length_vel_nb():
    # overriding length, velocity, numberdensity
    overrides = dict(length_unit=length_unit, velocity_unit=velocity_unit,
                     numberdensity_unit=numberdensity_unit)
    ds = data_dir_load(khi_cartesian_2D, kwargs={'units_override': overrides})
    _assert_normalisations_equal(ds)

@requires_file(khi_cartesian_2D)
def test_normalisations_length_vel_mass():
    # overriding length, velocity, mass
    overrides = dict(length_unit=length_unit, velocity_unit=velocity_unit,
                     mass_unit=mass_unit)
    ds = data_dir_load(khi_cartesian_2D, kwargs={'units_override': overrides})
    _assert_normalisations_equal(ds)

@requires_file(khi_cartesian_2D)
def test_normalisations_default():
    # test default normalisations, without overrides
    ds = data_dir_load(khi_cartesian_2D)
    assert_allclose_units(ds.length_unit, YTQuantity(1, 'cm'))
    assert_allclose_units(ds.numberdensity_unit, YTQuantity(1, 'cm**-3'))
    assert_allclose_units(ds.temperature_unit, YTQuantity(1, 'K'))
    assert_allclose_units(ds.density_unit, YTQuantity(2.341670657200000e-24, 'g*cm**-3'))
    assert_allclose_units(ds.mass_unit, YTQuantity(2.341670657200000e-24, 'g'))
    assert_allclose_units(ds.velocity_unit, YTQuantity(1.164508387441102e+04, 'cm*s**-1'))
    assert_allclose_units(ds.pressure_unit, YTQuantity(3.175492240000000e-16, 'dyn*cm**-2'))
    assert_allclose_units(ds.time_unit, YTQuantity(8.587314705370271e-05, 's'))
    assert_allclose_units(ds.magnetic_unit, YTQuantity(6.316993934686148e-08, 'gauss'))

@requires_file(khi_cartesian_2D)
def test_normalisations_too_many_args():
    # test forbidden case: too many arguments (max 3 are allowed)
    overrides = dict(length_unit=length_unit, numberdensity_unit=numberdensity_unit,
                     temperature_unit=temperature_unit, time_unit=time_unit)
    with assert_raises(ValueError):
        data_dir_load(khi_cartesian_2D, kwargs={'units_override': overrides})

@requires_file(khi_cartesian_2D)
def test_normalisations_vel_and_length():
    # test forbidden case: both velocity and temperature are specified as overrides
    overrides = dict(length_unit=length_unit, velocity_unit=velocity_unit,
                     temperature_unit=temperature_unit)
    with assert_raises(ValueError):
        data_dir_load(khi_cartesian_2D, kwargs={'units_override': overrides})
