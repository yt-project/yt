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

import yt
from yt.frontends.enzo.api import EnzoDataset
from yt.frontends.enzo.fields import NODAL_FLAGS
from yt.testing import assert_allclose_units
from yt.visualization.plot_window import SlicePlot

import framework as fw
import utils


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


#============================================
#                 TestEnzo
#============================================
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
    def test_toro1d(self):
        # Load data
        ds = utils.data_dir_load(toro1d)
        # Make sure we have the right data file
        assert str(ds) == 'data0001'
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "max"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "density"]
        fields = ds.field_list
        # Run the small_patch_amr test suite
        hashes = self.small_patch_amr(ds, fields, weights, axes, ds_objs)
        # Save or compare answer
        utils.handle_hashes(self.save_dir, 'toro1d', hashes, self.answer_store)

    #-----
    # test_kh2d
    #-----
    def test_kh2d(self):
        # Load data
        ds = utils.data_dir_load(kh2d)
        # Make sure we're dealing with the right file
        assert str(ds) == 'DD0011'
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "max"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "density"]
        fields = ds.field_list
        # Run the small_patch_amr test suite
        hashes = self.small_patch_amr(ds, fields, weights, axes, ds_objs)
        # Save or compare answer
        utils.handle_hashes(self.save_dir, 'kh2d', hashes, self.answer_store)

    #-----
    # test_moving7
    #-----
    def test_moving7(self):
        # Load data
        ds = utils.data_dir_load(m7)
        # Make sure we have the right data file
        assert str(ds) == 'moving7_0010'
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "max"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "density"]
        fields = ("temperature", "density", "velocity_magnitude",
            "velocity_divergence"
        )
        # Run the small_patch_amr test suite
        hashes = self.small_patch_amr(ds, fields, weights, axes, ds_objs)
        # Save or compare answer
        utils.handle_hashes(self.save_dir, 'moving7', hashes, self.answer_store)

    #-----
    # test_galaxy0030
    #-----
    def test_galaxy0030(self):
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
        ds = utils.data_dir_load(g30)
        # Make sure we have the right data file
        assert str(ds) == 'galaxy0030'
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "max"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "density"]
        fields = ("temperature", "density", "velocity_magnitude",
                   "velocity_divergence")
        # Color conservation test
        self.color_conservation_test(ds)
        # Run the big patch amr test suite
        hashes = self.big_patch_amr(ds, fields, weights, axes, ds_objs)
        # Save or compare answer
        utils.handle_hashes(self.save_dir, 'galaxy0030', hashes, self.answer_store)

    #-----
    # test_simulated_halo_mass_function
    #-----
    def test_simulated_halo_mass_function(self):
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
        # Load in the data
        ds = utils.data_dir_load(enzotiny)
        # Set up hex digest
        shmf_hd = b''
        # Loop over the different finders
        for finder in ["fof", "hop"]:
            shmf_hd += self.simulated_halo_mass_function_test(ds, finder)
        hashes = {}
        hashes['simulated_halo_mass_function'] = utils.generate_hash(shmf_hd)
        utils.handle_hashes(self.save_dir, 'shmf-enzotiny', hashes, self.answer_store)

    #-----
    # test_analytic_halo_mass_function
    #-----
    def test_analytic_halo_mass_function(self):
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
        # Load in the data
        ds = utils.data_dir_load(enzotiny)
        # Set up hex digest
        ahmf_hd = b''
        # Loop over the different finders
        for fit in range(1,6):
            ahmf_hd += self.analytic_halo_mass_function_test(ds, fit)
        hashes = {}
        hashes['analytic_halo_mass_function'] = utils.generate_hash(ahmf_hd)
        utils.handle_hashes(self.save_dir, 'ahmf-enzotiny', hashes, self.answer_store)

    #-----
    # test_ecp
    #-----
    def test_ecp(self):
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
        ds = utils.data_dir_load(ecp)
        self.color_conservation_test(ds)

    #-----
    # test_units_override
    #-----
    def test_units_override(self):
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
        yt.testing.units_override_check(enzotiny)

    #-----
    # test_nuclei_density_fields
    #-----
    def test_nuclei_density_fields(self):
        ds = utils.data_dir_load(ecp)
        ad = ds.all_data()
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
    # test_enzo_dataset
    #-----
    def test_enzo_dataset(self):
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
        assert isinstance(utils.data_dir_load(enzotiny), EnzoDataset)

    #-----
    # test_active_particle_dataset
    #-----
    def test_active_particle_datasets(self):
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
        # Read data for two sphere test
        two_sph = utils.data_dir_load(two_sphere_test)
        # Read data for active particle cosmology
        apcos = utils.data_dir_load(active_particle_cosmology)
        # Set up lists for comparison
        pfields = ['GridID', 'creation_time', 'dynamical_time',
                   'identifier', 'level', 'metallicity', 'particle_mass']
        pfields += ['particle_position_%s' % d for d in 'xyz']
        pfields += ['particle_velocity_%s' % d for d in 'xyz']
        acc_part_fields = \
            [('AccretingParticle', pf) for pf in ['AccretionRate'] + pfields]
        real_acc_part_fields = sorted(
            [f for f in two_sph.field_list if f[0] == 'AccretingParticle'])
        # Set up lists for comparison
        apcos_fields = [('CenOstriker', pf) for pf in pfields]
        real_apcos_fields = sorted(
            [f for f in apcos.field_list if f[0] == 'CenOstriker'])
        apcos_pcounts = {'CenOstriker': 899755, 'DarkMatter': 32768}
        assert 'AccretingParticle' in two_sph.particle_types_raw
        assert 'io' not in two_sph.particle_types_raw
        assert 'all' in two_sph.particle_types
        assert len(two_sph.particle_unions) ==  1
        assert acc_part_fields == real_acc_part_fields
        assert ['CenOstriker', 'DarkMatter'] == apcos.particle_types_raw
        assert 'all' in apcos.particle_unions
        assert apcos_fields == real_apcos_fields
        assert apcos.particle_type_counts == apcos_pcounts

    #-----
    # test_face_centered_mhdct_fields
    #-----
    def test_face_centered_mhdct_fields(self):
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
        ds = utils.data_dir_load(mhdctot)
        ad = ds.all_data()
        grid = ds.index.grids[0]
        dims = ds.domain_dimensions
        dims_prod = dims.prod()
        for field, flag in NODAL_FLAGS.items():
            assert ad[field].shape == (dims_prod, 2*sum(flag))
            assert grid[field].shape == tuple(dims) + (2*sum(flag),)
        # Average of face-centered fields should be the same as
        # cell-centered field
        assert (ad['BxF'].sum(axis=-1)/2 == ad['Bx']).all()
        assert (ad['ByF'].sum(axis=-1)/2 == ad['By']).all()
        assert (ad['BzF'].sum(axis=-1)/2 == ad['Bz']).all()

    #-----
    # test_deeply_nested_zoom
    #-----
    def test_deeply_nested_zoom(self):
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
        # Load data
        ds = utils.data_dir_load(dnz)
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
        assert max([g['density'].max() for g in ds.index.grids]) == v

    #-----
    # test_2d_grid_shape
    #-----
    def test_2d_grid_shape(self):
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
        ds = utils.data_dir_load(kh2d)
        g = ds.index.grids[1]
        assert g['density'].shape == (128, 100, 1)

    #-----
    # test_nonzero_omega_radiation
    #-----
    def test_nonzero_omega_radiation(self):
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
        ds = utils.data_dir_load(p3mini)
        err_msg = "Simulation time not consistent with cosmology calculator."
        t_from_z = ds.cosmology.t_from_z(ds.current_redshift)
        tratio = ds.current_time / t_from_z
        assert ds.omega_radiation == ds.cosmology.omega_radiation
        np.testing.assert_almost_equal(tratio, 1, 4, err_msg=err_msg)
