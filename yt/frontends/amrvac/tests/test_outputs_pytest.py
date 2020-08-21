import numpy as np  # NOQA
import pytest

import yt  # NOQA
from yt.frontends.amrvac.api import AMRVACDataset, AMRVACGrid
from yt.testing import assert_allclose_units, assert_raises, requires_file
from yt.units import YTQuantity
from yt.utilities.answer_testing import utils
from yt.utilities.answer_testing.answer_tests import (
    field_values,
    grid_hierarchy,
    grid_values,
    parentage_relationships,
    projection_values,
)

# Test parameters
blastwave_spherical_2D = "amrvac/bw_2d0000.dat"
khi_cartesian_2D = "amrvac/kh_2d0000.dat"
khi_cartesian_3D = "amrvac/kh_3D0000.dat"
jet_cylindrical_25D = "amrvac/Jet0003.dat"
riemann_cartesian_175D = "amrvac/R_1d0005.dat"
blastwave_cartesian_3D = "amrvac/bw_3d0000.dat"
blastwave_polar_2D = "amrvac/bw_polar_2D0000.dat"
blastwave_cylindrical_3D = "amrvac/bw_cylindrical_3D0000.dat"
rmi_cartesian_dust_2D = "amrvac/Richtmyer_Meshkov_dust_2D/RM2D_dust_Kwok0000.dat"


ds_list = [
    blastwave_spherical_2D,
    khi_cartesian_2D,
    khi_cartesian_3D,
    jet_cylindrical_25D,
    riemann_cartesian_175D,
    blastwave_cartesian_3D,
    blastwave_polar_2D,
    blastwave_cylindrical_3D,
    rmi_cartesian_dust_2D,
]
a_list = [0, 1, 2]
w_list = [None, "density"]
dso_list = [None, ("sphere", ("max", (0.1, "unitary")))]


# Tests for units: verify that overriding certain units yields the correct derived units.
# The following are correct normalisations based on length, numberdensity and temperature
length_unit = (1e9, "cm")
numberdensity_unit = (1e9, "cm**-3")
temperature_unit = (1e6, "K")
density_unit = (2.341670657200000e-15, "g*cm**-3")
mass_unit = (2.341670657200000e12, "g")
velocity_unit = (1.164508387441102e07, "cm*s**-1")
pressure_unit = (3.175492240000000e-01, "dyn*cm**-2")
time_unit = (8.587314705370271e01, "s")
magnetic_unit = (1.997608879907716, "gauss")


def _get_fields_to_check(fname):
    # This function is called during test collection. If this frontend
    # is not being run, and therefore the data isn't present, this try
    # except block prevents pytest from failing needlessly
    try:
        ds = utils.data_dir_load(fname)
        fields = [("gas", "density"), ("gas", "velocity_magnitude")]
        field_ids = ["density", "velocity_magnitude"]
        raw_fields_labels = [fname for ftype, fname in ds.field_list]
        if "b1" in raw_fields_labels:
            fields.append(("gas", "magnetic_energy_density"))
            field_ids.append("magnetic_energy_density")
        if "e" in raw_fields_labels:
            fields.append(("gas", "energy_density"))
            field_ids.append("energy_density")
        if "rhod1" in raw_fields_labels:
            fields.append(("gas", "total_dust_density"))
            field_ids.append("total_dust_density")
            # note : not hitting dust velocity fields
        # return [fields, field_ids]
        return [ds, fields]
    except FileNotFoundError:
        return [None, [None,]]


def get_pairs():
    f_list = [_get_fields_to_check(fname) for fname in ds_list]
    pairs = [(i[0], f) for i in f_list for f in i[1]]
    return pairs


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


@pytest.mark.answer_test
class TestAMRVAC:
    answer_file = None

    @requires_file(khi_cartesian_2D)
    def test_AMRVACDataset(self):
        assert isinstance(utils.data_dir_load(khi_cartesian_2D), AMRVACDataset)

    @utils.requires_ds(blastwave_cartesian_3D)
    def test_domain_size(self):
        # "Check for correct box size, see bw_3d.par"
        ds = utils.data_dir_load(blastwave_cartesian_3D)
        for lb in ds.domain_left_edge:
            assert int(lb) == 0
        for rb in ds.domain_right_edge:
            assert int(rb) == 2
        for w in ds.domain_width:
            assert int(w) == 2

    @requires_file(blastwave_cartesian_3D)
    def test_grid_attributes(self):
        # "Check various grid attributes"
        ds = utils.data_dir_load(blastwave_cartesian_3D)
        grids = ds.index.grids
        assert ds.index.max_level == 2
        for g in grids:
            assert isinstance(g, AMRVACGrid)
            assert isinstance(g.LeftEdge, yt.units.yt_array.YTArray)
            assert isinstance(g.RightEdge, yt.units.yt_array.YTArray)
            assert isinstance(g.ActiveDimensions, np.ndarray)
            assert isinstance(g.Level, (np.int32, np.int64, int))

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", ds_list, indirect=True)
    def test_gh_pr(self, ds):
        self.hashes.update({"grid_hierarchy": grid_hierarchy(ds)})
        self.hashes.update({"parentage_relationships": parentage_relationships(ds)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f", get_pairs(), indirect=True)
    def test_gv(self, f, ds):
        if ds and f is None:
            pytest.skip("Data not found.")
        self.hashes.update({"grid_values": grid_values(ds, f)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f", get_pairs(), indirect=True)
    @pytest.mark.parametrize("d", dso_list, indirect=True)
    def test_fv(self, d, f, ds):
        if ds and f is None:
            pytest.skip("Data not found.")
        self.hashes.update({"field_values": field_values(ds, f, d)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f", get_pairs(), indirect=True)
    @pytest.mark.parametrize("d", dso_list, indirect=True)
    @pytest.mark.parametrize("a", a_list, indirect=True)
    @pytest.mark.parametrize("w", w_list, indirect=True)
    def test_pv(self, a, d, w, f, ds):
        if ds and f is None:
            pytest.skip("Data not found.")
        self.hashes.update({"projection_values": projection_values(ds, a, f, w, d)})

    @requires_file(khi_cartesian_2D)
    def test_normalisations_length_temp_nb(self):
        # overriding length, temperature, numberdensity
        overrides = dict(
            length_unit=length_unit,
            temperature_unit=temperature_unit,
            numberdensity_unit=numberdensity_unit,
        )
        ds = utils.data_dir_load(khi_cartesian_2D, kwargs={"units_override": overrides})
        _assert_normalisations_equal(ds)

    @requires_file(khi_cartesian_2D)
    def test_normalisations_length_temp_mass(self):
        # overriding length, temperature, mass
        overrides = dict(
            length_unit=length_unit,
            temperature_unit=temperature_unit,
            mass_unit=mass_unit,
        )
        ds = utils.data_dir_load(khi_cartesian_2D, kwargs={"units_override": overrides})
        _assert_normalisations_equal(ds)

    @requires_file(khi_cartesian_2D)
    def test_normalisations_length_time_mass(self):
        # overriding length, time, mass
        overrides = dict(
            length_unit=length_unit, time_unit=time_unit, mass_unit=mass_unit
        )
        ds = utils.data_dir_load(khi_cartesian_2D, kwargs={"units_override": overrides})
        _assert_normalisations_equal(ds)

    @requires_file(khi_cartesian_2D)
    def test_normalisations_length_vel_nb(self):
        # overriding length, velocity, numberdensity
        overrides = dict(
            length_unit=length_unit,
            velocity_unit=velocity_unit,
            numberdensity_unit=numberdensity_unit,
        )
        ds = utils.data_dir_load(khi_cartesian_2D, kwargs={"units_override": overrides})
        _assert_normalisations_equal(ds)

    @requires_file(khi_cartesian_2D)
    def test_normalisations_length_vel_mass(self):
        # overriding length, velocity, mass
        overrides = dict(
            length_unit=length_unit, velocity_unit=velocity_unit, mass_unit=mass_unit
        )
        ds = utils.data_dir_load(khi_cartesian_2D, kwargs={"units_override": overrides})
        _assert_normalisations_equal(ds)

    @requires_file(khi_cartesian_2D)
    def test_normalisations_default(self):
        # test default normalisations, without overrides
        ds = utils.data_dir_load(khi_cartesian_2D)
        assert_allclose_units(ds.length_unit, YTQuantity(1, "cm"))
        assert_allclose_units(ds.numberdensity_unit, YTQuantity(1, "cm**-3"))
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
        assert_allclose_units(
            ds.magnetic_unit, YTQuantity(6.316993934686148e-08, "gauss")
        )

    @requires_file(khi_cartesian_2D)
    def test_normalisations_too_many_args(self):
        # test forbidden case: too many arguments (max 3 are allowed)
        overrides = dict(
            length_unit=length_unit,
            numberdensity_unit=numberdensity_unit,
            temperature_unit=temperature_unit,
            time_unit=time_unit,
        )
        with assert_raises(ValueError):
            utils.data_dir_load(khi_cartesian_2D, kwargs={"units_override": overrides})

    @requires_file(khi_cartesian_2D)
    def test_normalisations_vel_and_length(self):
        # test forbidden case: both velocity and temperature are specified as overrides
        overrides = dict(
            length_unit=length_unit,
            velocity_unit=velocity_unit,
            temperature_unit=temperature_unit,
        )
        with assert_raises(ValueError):
            utils.data_dir_load(khi_cartesian_2D, kwargs={"units_override": overrides})
