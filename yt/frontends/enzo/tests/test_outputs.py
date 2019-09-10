"""
Title:   test_enzo.py
Purpose: Contains Enzo frontend tests
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
from collections import OrderedDict

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
import yt.utilities.answer_testing.framework as fw
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


# File to store answers (will eventually go in an ini file)
answer_file = 'enzo_answers.yaml'


#============================================
#                 TestEnzo
#============================================
@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
class TestEnzo(fw.AnswerTest):
    """
    Container for Enzo answer tests.

    Attributes:
    -----------
        pass

    Methods:
    --------
        pass
    """
    #-----
    # test_toro1d
    #-----
    @utils.requires_ds(toro1d)
    def test_toro1d(self, ds_toro1d):
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "max"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "density"]
        fields = ds_toro1d.field_list
        # Run the small_patch_amr test suite
        hashes = self.small_patch_amr(ds_toro1d, fields, weights, axes, ds_objs)
        # Add test name as a key
        hashes = {'toro1d' : hashes}
        # Save or compare answer
        utils.handle_hashes(self.save_dir, answer_file, hashes, self.answer_store)

    #-----
    # test_kh2d
    #-----
    @utils.requires_ds(kh2d)
    def test_kh2d(self, ds_kh2d):
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "max"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "density"]
        fields = ds_kh2d.field_list
        # Run the small_patch_amr test suite
        hashes = self.small_patch_amr(ds_kh2d, fields, weights, axes, ds_objs)
        hashes = {'kh2d' : hashes}
        # Save or compare answer
        utils.handle_hashes(self.save_dir, answer_file, hashes, self.answer_store)

    #-----
    # test_moving7
    #-----
    @utils.requires_ds(m7)
    def test_moving7(self, ds_m7):
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "max"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "density"]
        fields = ("temperature", "density", "velocity_magnitude",
            "velocity_divergence"
        )
        # Run the small_patch_amr test suite
        hashes = self.small_patch_amr(ds_m7, fields, weights, axes, ds_objs)
        hashes = {'moving7' : hashes}
        # Save or compare answer
        utils.handle_hashes(self.save_dir, answer_file, hashes, self.answer_store)

    #-----
    # test_galaxy0030
    #-----
    @pytest.mark.skipif(not pytest.config.getvalue('--answer-big-data'),
        reason="Skipping because --answer-big-data was not set."
    )
    @utils.requires_ds(g30)
    def test_galaxy0030(self, ds_g30):
        """
        Tests the galaxy_0030 big data simulation.

        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "max"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "density"]
        fields = ("temperature", "density", "velocity_magnitude",
                   "velocity_divergence")
        # Color conservation test
        self.color_conservation_test(ds_g30)
        # Run the big patch amr test suite
        hashes = self.big_patch_amr(ds_g30, fields, weights, axes, ds_objs)
        hashes = {'galaxy0030' : hashes}
        # Save or compare answer
        utils.handle_hashes(self.save_dir, answer_file, hashes, self.answer_store)

    #-----
    # test_simulated_halo_mass_function
    #-----
    @utils.requires_ds(enzotiny)
    def test_simulated_halo_mass_function(self, ds_enzotiny):
        """
        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
        # Set up hex digest
        hashes = OrderedDict()
        hashes['simulated_halo_mass_function'] = OrderedDict()
        # Loop over the different finders
        for finder in ["fof", "hop"]:
            shmf_hd = utils.generate_hash(
                self.simulated_halo_mass_function_test(ds_enzotiny, finder)
            )
            hashes['simulated_halo_mass_function'][finder] = shmf_hd
        hashes = {'simulated_halo_mass_function' : hashes}
        utils.handle_hashes(self.save_dir, answer_file, hashes, self.answer_store)

    #-----
    # test_analytic_halo_mass_function
    #-----
    @utils.requires_ds(enzotiny)
    def test_analytic_halo_mass_function(self, ds_enzotiny):
        """
        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
        # Set up hex digest
        hashes = OrderedDict()
        hashes['analytic_halo_mass_function'] = OrderedDict()
        # Loop over the different finders
        for fit in range(1,6):
            ahmf_hd = utils.generate_hash(
                self.analytic_halo_mass_function_test(ds_enzotiny, fit)
            )
            hashes['analytic_halo_mass_function'][str(fit)] = ahmf_hd
        hashes = {'analytic_halo_mass_function' : hashes}
        utils.handle_hashes(self.save_dir, answer_file, hashes, self.answer_store)

    #-----
    # test_ecp
    #-----
    @pytest.mark.skipif(not pytest.config.getvalue('--answer-big-data'),
        reason="Skipping test_jet because --answer-big-data was not set."
    )
    @utils.requires_ds(ecp)
    def test_ecp(self, ds_ecp):
        """
        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
        self.color_conservation_test(ds_ecp)

    #-----
    # test_units_override
    #-----
    @requires_file(enzotiny)
    def test_units_override(self, ds_enzotiny):
        """
        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
        units_override_check(ds_enzotiny, enzotiny)

    #-----
    # test_nuclei_density_fields
    #-----
    @pytest.mark.skipif(not pytest.config.getvalue('--answer-big-data'),
        reason="Skipping test_jet because --answer-big-data was not set."
    )
    @utils.requires_ds(ecp)
    def test_nuclei_density_fields(self, ds_ecp):
        ad = ds_ecp.all_data()
        # Hydrogen
        hd1 = utils.generate_hash(ad["H_nuclei_density"].tostring())
        hd2 = utils.generate_hash((ad["H_number_density"] +
            ad["H_p1_number_density"]).tostring())
        assert hd1 == hd2
        hd1 = utils.generate_hash(ad["He_nuclei_density"].tostring())
        hd2 = utils.generate_hash((ad["He_number_density"] +
            ad["He_p1_number_density"] +
            ad["He_p2_number_density"]).tostring()
        )
        assert hd1 == hd2

    #-----
    # test_EnzoDataset
    #-----
    @requires_file(enzotiny)
    def test_EnzoDataset(self, ds_enzotiny):
        """
        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
        assert isinstance(ds_enzotiny, EnzoDataset)

    #-----
    # test_active_particle_dataset
    #-----
    @requires_file(two_sphere_test)
    @requires_file(active_particle_cosmology)
    def test_active_particle_datasets(self, ds_two_sphere_test, ds_active_particle_cosmology):
        """
        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
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

    #-----
    # test_face_centered_mhdct_fields
    #-----
    @requires_file(mhdctot)
    def test_face_centered_mhdct_fields(self, ds_mhdctot):
        """
        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
        ad = ds_mhdctot.all_data()
        grid = ds_mhdctot.index.grids[0]
        dims = ds_mhdctot.domain_dimensions
        dims_prod = dims.prod()
        for field, flag in NODAL_FLAGS.items():
            assert_equal(ad[field].shape, (dims_prod, 2*sum(flag)))
            assert_equal(grid[field].shape, tuple(dims) + (2*sum(flag),))
        # Average of face-centered fields should be the same as
        # cell-centered field
        assert (ad['BxF'].sum(axis=-1)/2 == ad['Bx']).all()
        assert (ad['ByF'].sum(axis=-1)/2 == ad['By']).all()
        assert (ad['BzF'].sum(axis=-1)/2 == ad['Bz']).all()

    #-----
    # test_deeply_nested_zoom
    #-----
    @utils.requires_ds(dnz)
    def test_deeply_nested_zoom(self, ds_dnz):
        """
        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
        # Carefully chosen to just barely miss a grid in the middle of
        # the image
        center = [0.4915073260199302, 0.5052605316800006, 0.4905805557500548]
        plot = SlicePlot(ds_dnz, 'z', 'density', width=(0.001, 'pc'),
                         center=center)
        image = plot.frb['density']
        assert (image > 0).all()
        v, c = ds_dnz.find_max('density')
        assert_allclose_units(v, ds_dnz.quan(0.005878286377124154, 'g/cm**3'))
        c_actual = [0.49150732540021, 0.505260532936791, 0.49058055816398]
        c_actual = ds_dnz.arr(c_actual, 'code_length')
        assert_allclose_units(c, c_actual)
        assert_equal(max([g['density'].max() for g in ds_dnz.index.grids]), v)

    #-----
    # test_2d_grid_shape
    #-----
    @requires_file(kh2d)
    def test_2d_grid_shape(self, ds_kh2d):
        """
        See issue #1601: we want to make sure that accessing data on
        a grid object returns a 3D array with a dummy dimension

        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
        g = ds_kh2d.index.grids[1]
        assert g['density'].shape == (128, 100, 1)

    #-----
    # test_nonzero_omega_radiation
    #-----
    @requires_file(p3mini)
    def test_nonzero_omega_radiation(self, ds_p3mini):
        """
        Test support for non-zero omega_radiation cosmologies.

        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
        err_msg = "Simulation time not consistent with cosmology calculator."
        t_from_z = ds_p3mini.cosmology.t_from_z(ds_p3mini.current_redshift)
        tratio = ds_p3mini.current_time / t_from_z
        assert_equal(ds_p3mini.omega_radiation, ds_p3mini.cosmology.omega_radiation)
        assert_almost_equal(tratio, 1, 4, err_msg=err_msg)
