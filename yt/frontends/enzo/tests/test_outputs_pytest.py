"""
Title:   test_enzo.py
Purpose: Contains Enzo frontend tests
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import numpy as np
import pytest

from yt.frontends.enzo.api import EnzoDataset
from yt.frontends.enzo.fields import NODAL_FLAGS
from yt.testing import \
    assert_allclose_units, \
    assert_almost_equal, \
    assert_equal, \
    requires_file, \
    units_override_check
from yt.visualization.plot_window import SlicePlot
from yt.utilities.answer_testing.answer_tests import \
    big_patch_amr, \
    small_patch_amr
from yt.utilities.answer_testing import utils


# Files containing data to be used in tests. Paths are relative to
# yt test_data_dir
toro1d = "ToroShockTube/DD0001/data0001"
kh2d = "EnzoKelvinHelmholtz/DD0011/DD0011"
m7 = "DD0010/moving7_0010"
g30 = "IsolatedGalaxy/galaxy0030/galaxy0030"
enzotiny = "enzo_tiny_cosmology/DD0046/DD0046"
ecp = "enzo_cosmology_plus/DD0046/DD0046"
two_sphere_test = 'ActiveParticleTwoSphere/DD0011/DD0011'
active_particle_cosmology = 'ActiveParticleCosmology/DD0046/DD0046'
mhdctot = "MHDCTOrszagTang/DD0004/data0004"
dnz = "DeeplyNestedZoom/DD0025/data0025"
p3mini = "PopIII_mini/DD0034/DD0034"


def color_conservation_test(ds):
    species_names = ds.field_info.species_names
    dd = ds.all_data()
    dens_yt = dd["density"].copy()
    # Enumerate our species here
    for s in sorted(species_names):
        if s == "El": continue
        dens_yt -= dd["%s_density" % s]
    dens_yt -= dd["metal_density"]
    delta_yt = np.abs(dens_yt / dd["density"])
    # Now we compare color conservation to Enzo's color conservation
    dd = ds.all_data()
    dens_enzo = dd["Density"].copy()
    for f in sorted(ds.field_list):
        ff = f[1]
        if not ff.endswith("_Density"):
            continue
        start_strings = ["Electron_", "SFR_", "Forming_Stellar_",
                         "Dark_Matter", "Star_Particle_"]
        if any([ff.startswith(ss) for ss in start_strings]):
            continue
        dens_enzo -= dd[f]
    delta_enzo = np.abs(dens_enzo / dd["Density"])
    np.testing.assert_almost_equal(delta_yt, delta_enzo)


@pytest.mark.answer_test
@pytest.mark.usefixtures('answer_file')
class TestEnzo:

    @pytest.mark.usefixtures('hashing')
    @pytest.mark.parametrize('ds', [toro1d], indirect=True)
    def test_toro1d(self, a, d, w, f, ds):
        self.hashes.update(small_patch_amr(ds, f, w, a, d))

    @pytest.mark.usefixtures('hashing')
    @pytest.mark.parametrize('ds', [kh2d], indirect=True)
    def test_kh2d(self, a, d, w, f, ds):
        self.hashes.update(small_patch_amr(ds, f, w, a, d))

    @pytest.mark.usefixtures('hashing')
    @pytest.mark.parametrize('ds', [m7], indirect=True)
    def test_moving7(self, a, d, w, f, ds):
        self.hashes.update(small_patch_amr(ds, f, w, a, d))

    @pytest.mark.big_data
    @pytest.mark.usefixtures('hashing')
    @pytest.mark.parametrize('ds', [g30], indirect=True)
    def test_galaxy0030(self, a, d, w, f, ds):
        color_conservation_test(ds)
        assert_equal(ds.particle_type_counts, {'io': 1124453})
        self.hashes.update(big_patch_amr(ds, f, w, a, d))

    @pytest.mark.big_data
    @pytest.mark.parametrize('ds', [ecp], indirect=True)
    def test_ecp(self, ds):
        color_conservation_test(ds)

    @requires_file(enzotiny)
    def test_units_override(self):
        units_override_check(enzotiny)

    @pytest.mark.big_data
    @pytest.mark.parametrize('ds', [ecp], indirect=True)
    def test_nuclei_density_fields(self, ds):
        ad = ds.all_data()
        assert_array_equal(ad["H_nuclei_density"],
                           (ad["H_p0_number_density"] + ad["H_p1_number_density"]))
        assert_array_equal(ad["He_nuclei_density"],
            (ad["He_p0_number_density"] +
             ad["He_p1_number_density"] + ad["He_p2_number_density"]))

    @pytest.mark.parametrize('ds', [enzotiny], indirect=True)
    def test_EnzoDataset(self, ds):
        assert isinstance(ds, EnzoDataset)

    @requires_file(two_sphere_test)
    @requires_file(active_particle_cosmology)
    def test_active_particle_datasets(self):
        ds_two_sphere_test = utils.data_dir_load(two_sphere_test)
        ds_active_particle_cosmology = utils.data_dir_load(active_particle_cosmology)
        # Set up lists for comparison
        pfields = ['GridID', 'creation_time', 'dynamical_time',
                   'identifier', 'level', 'metallicity', 'particle_mass']
        pfields += ['particle_position_%s' % d for d in 'xyz']
        pfields += ['particle_velocity_%s' % d for d in 'xyz']
        acc_part_fields = \
            [('AccretingParticle', pf) for pf in ['AccretionRate'] + pfields]
        real_acc_part_fields = sorted(
            [f for f in ds_two_sphere_test.field_list if f[0] == 'AccretingParticle'])
        # Set up lists for comparison
        apcos_fields = [('CenOstriker', pf) for pf in pfields]
        real_apcos_fields = sorted(
            [f for f in ds_active_particle_cosmology.field_list if f[0] == 'CenOstriker'])
        apcos_pcounts = {'CenOstriker': 899755, 'DarkMatter': 32768}
        assert 'AccretingParticle' in ds_two_sphere_test.particle_types_raw
        assert 'io' not in ds_two_sphere_test.particle_types_raw
        assert 'all' in ds_two_sphere_test.particle_types
        assert_equal(len(ds_two_sphere_test.particle_unions), 1)
        assert_equal(acc_part_fields, real_acc_part_fields)
        assert_equal(['CenOstriker', 'DarkMatter'], ds_active_particle_cosmology.particle_types_raw)
        assert 'all' in ds_active_particle_cosmology.particle_unions
        assert_equal(apcos_fields, real_apcos_fields)
        assert_equal(ds_active_particle_cosmology.particle_type_counts, apcos_pcounts)

    @pytest.mark.parametrize('ds', [mhdctot], indirect=True)
    def test_face_centered_mhdct_fields(self, ds):
        ad = ds.all_data()
        grid = ds.index.grids[0]
        dims = ds.domain_dimensions
        dims_prod = dims.prod()
        for field, flag in NODAL_FLAGS.items():
            assert_equal(ad[field].shape, (dims_prod, 2*sum(flag)))
            assert_equal(grid[field].shape, tuple(dims) + (2*sum(flag),))
        # Average of face-centered fields should be the same as
        # cell-centered field
        assert (ad['BxF'].sum(axis=-1)/2 == ad['Bx']).all()
        assert (ad['ByF'].sum(axis=-1)/2 == ad['By']).all()
        assert (ad['BzF'].sum(axis=-1)/2 == ad['Bz']).all()

    @pytest.mark.parametrize('ds', [dnz], indirect=True)
    def test_deeply_nested_zoom(self, ds):
        # Carefully chosen to just barely miss a grid in the middle of
        # the image
        center = [0.4915073260199302, 0.5052605316800006, 0.4905805557500548]
        plot = SlicePlot(ds, 'z', 'density', width=(0.001, 'pc'),
                         center=center)
        image = plot.frb['density']
        assert (image > 0).all()
        v, c = ds.find_max('density')
        assert_allclose_units(v, ds.quan(0.005878286377124154, 'g/cm**3'))
        c_actual = [0.49150732540021, 0.505260532936791, 0.49058055816398]
        c_actual = ds.arr(c_actual, 'code_length')
        assert_allclose_units(c, c_actual)
        assert_equal(max([g['density'].max() for g in ds.index.grids]), v)

    @pytest.mark.parametrize('ds', [kh2d], indirect=True)
    def test_2d_grid_shape(self, ds):
        r"""See issue #1601: we want to make sure that accessing data on
        a grid object returns a 3D array with a dummy dimension
        """
        g = ds.index.grids[1]
        assert g['density'].shape == (128, 100, 1)

    @pytest.mark.parametrize('ds', [p3mini], indirect=True)
    def test_nonzero_omega_radiation(self, ds):
        r"""Test support for non-zero omega_radiation cosmologies.
        """
        err_msg = "Simulation time not consistent with cosmology calculator."
        t_from_z = ds.cosmology.t_from_z(ds.current_redshift)
        tratio = ds.current_time / t_from_z
        assert_equal(ds.omega_radiation, ds.cosmology.omega_radiation)
        assert_almost_equal(tratio, 1, 4, err_msg=err_msg)
