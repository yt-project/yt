"""
Title: test_boxlib.py
Purpose: Boxlib frontend tests
Notes:
    Copyright (c) 2017, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import numpy as np

from yt.testing import \
    assert_equal, \
    requires_file, \
    units_override_check
from yt.frontends.boxlib.api import \
    OrionDataset, \
    NyxDataset, \
    WarpXDataset, \
    CastroDataset, \
    MaestroDataset

import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils


# Test data
radadvect = "RadAdvect/plt00000"
rt = "RadTube/plt00500"
star = "StarParticles/plrd01000"
LyA = "Nyx_LyA/plt00000"
RT_particles = "RT_particles/plt00050"
langmuir = "LangmuirWave/plt00020_v2"
plasma = "PlasmaAcceleration/plt00030_v2"
beam = "GaussianBeam/plt03008"
raw_fields = "Laser/plt00015"
nyx_no_particles = "nyx_sedov_plt00086"
msubch = 'maestro_subCh_plt00248'


# We don't do anything needing ghost zone generation right now, because these
# are non-periodic datasets.
_orion_fields = ("temperature", "density", "velocity_magnitude")
_nyx_fields = ("Ne", "Temp", "particle_mass_density")
_warpx_fields = ("Ex", "By", "jz")
_castro_fields = ("Temp", "density", "particle_count")
_raw_fields = [('raw', 'Bx'), ('raw', 'Ey'), ('raw', 'jz')]


#============================================
#                TestBoxLib
#============================================
class TestBoxLib(fw.AnswerTest):
    """
    Container for boxlib frontend answer tests.

    Attributes:
    -----------
        pass

    Methods:
    --------
        pass
    """
    #-----
    # test_radavect
    #-----
    @utils.requires_ds(radadvect)
    def test_radadvect(self):
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
        ds = utils.data_dir_load(radadvect)
        assert_equal(str(ds), "plt00000")
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "max"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "density"]
        # Run the small_patch_amr test suite
        hashes = self.small_patch_amr(ds, _orion_fields, weights, axes, ds_objs)
        # Save or compare answer
        utils.handle_hashes(self.save_dir, 'boxlib-test-radadvect', hashes, self.answer_store)

    #-----
    # test_radtube
    #-----
    @utils.requires_ds(rt)
    def test_radtube(self, ds_rt):
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
        ds = ds_rt
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "max"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "density"]
        # Run the small_patch_amr test suite
        hashes = self.small_patch_amr(ds, _orion_fields, weights, axes, ds_objs)
        # Save or compare answer
        utils.handle_hashes(self.save_dir, 'boxlib-test-radtube', hashes, self.answer_store)

    #-----
    # test_star
    #-----
    @utils.requires_ds(star)
    def test_star(self, ds_star):
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
        ds = ds_star
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "max"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "density"]
        # Run the small_patch_amr test suite
        hashes = self.small_patch_amr(ds, _orion_fields, weights, axes, ds_objs)
        # Save or compare answer
        utils.handle_hashes(self.save_dir, 'boxlib-test-star', hashes, self.answer_store)

    #-----
    # test_LyA
    #-----
    @utils.requires_ds(LyA)
    def test_LyA(self, ds_LyA):
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
        ds = ds_LyA
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "c"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "Ne"]
        # Run the small_patch_amr test suite
        hashes = self.small_patch_amr(ds, _nyx_fields, weights, axes, ds_objs)
        # Save or compare answer
        utils.handle_hashes(self.save_dir, 'boxlib-test-LyA', hashes, self.answer_store)

    #-----
    # test_nyx_particle_io
    #-----
    @requires_file(LyA)
    def test_nyx_particle_io(self, ds_LyA):
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
        ds = ds_LyA
        grid = ds.index.grids[0]
        npart_grid_0 = 7908  # read directly from the header
        assert_equal(grid['particle_position_x'].size, npart_grid_0)
        assert_equal(grid['DM', 'particle_position_y'].size, npart_grid_0)
        assert_equal(grid['all', 'particle_position_z'].size, npart_grid_0)
        ad = ds.all_data()
        npart = 32768  # read directly from the header
        assert_equal(ad['particle_velocity_x'].size, npart)
        assert_equal(ad['DM', 'particle_velocity_y'].size, npart)
        assert_equal(ad['all', 'particle_velocity_z'].size, npart)
        assert(np.all(ad['particle_mass'] == ad['particle_mass'][0]))
        left_edge = ds.arr([0.0, 0.0, 0.0], 'code_length')
        right_edge = ds.arr([4.0, 4.0, 4.0], 'code_length')
        center = 0.5*(left_edge + right_edge)
        reg = ds.region(center, left_edge, right_edge)
        assert(np.all(np.logical_and(reg['particle_position_x'] <= right_edge[0], 
                                     reg['particle_position_x'] >= left_edge[0])))
        assert(np.all(np.logical_and(reg['particle_position_y'] <= right_edge[1], 
                                     reg['particle_position_y'] >= left_edge[1])))
        assert(np.all(np.logical_and(reg['particle_position_z'] <= right_edge[2], 
                                     reg['particle_position_z'] >= left_edge[2])))

    #-----
    # test_RT_particles
    #-----
    @utils.requires_ds(RT_particles)
    def test_RT_particles(self, ds_RT_particles):
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
        ds = ds_RT_particles
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "max"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "density"]
        # Run the small_patch_amr test suite
        hashes = self.small_patch_amr(ds, _castro_fields, weights, axes, ds_objs)
        # Save or compare answer
        utils.handle_hashes(self.save_dir, 'boxlib-test-RT-particles', hashes, self.answer_store)

    #-----
    # test_castro_particle_io
    #-----
    @requires_file(RT_particles)
    def test_castro_particle_io(self, ds_RT_particles):
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
        ds = ds_RT_particles
        grid = ds.index.grids[2]
        npart_grid_2 = 49  # read directly from the header
        assert_equal(grid['particle_position_x'].size, npart_grid_2)
        assert_equal(grid['Tracer', 'particle_position_y'].size, npart_grid_2)
        assert_equal(grid['all', 'particle_position_y'].size, npart_grid_2)
        ad = ds.all_data()
        npart = 49  # read directly from the header
        assert_equal(ad['particle_velocity_x'].size, npart)
        assert_equal(ad['Tracer', 'particle_velocity_y'].size, npart)
        assert_equal(ad['all', 'particle_velocity_y'].size, npart)
        left_edge = ds.arr([0.0, 0.0, 0.0], 'code_length')
        right_edge = ds.arr([0.25, 1.0, 1.0], 'code_length')
        center = 0.5*(left_edge + right_edge)
        reg = ds.region(center, left_edge, right_edge)
        assert(np.all(np.logical_and(reg['particle_position_x'] <= right_edge[0], 
                                     reg['particle_position_x'] >= left_edge[0])))
        assert(np.all(np.logical_and(reg['particle_position_y'] <= right_edge[1], 
                                     reg['particle_position_y'] >= left_edge[1])))

    #-----
    # test_langmuir
    #-----
    @utils.requires_ds(langmuir)
    def test_langmuir(self):
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
        ds = utils.data_dir_load(langmuir)
        assert_equal(str(ds), "plt00020_v2")
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "c"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "Ex"]
        # Run the small_patch_amr test suite
        hashes = self.small_patch_amr(ds, _warpx_fields, weights, axes, ds_objs)
        # Save or compare answer
        utils.handle_hashes(self.save_dir, 'boxlib-test-langmuir', hashes, self.answer_store)

    #-----
    # test_plasma
    #-----
    @utils.requires_ds(plasma)
    def test_plasma(self, ds_plasma):
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
        ds = ds_plasma
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "c"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "Ex"]
        # Run the small_patch_amr test suite
        hashes = self.small_patch_amr(ds, _warpx_fields, weights, axes, ds_objs)
        # Save or compare answer
        utils.handle_hashes(self.save_dir, 'boxlib-test-plasma', hashes, self.answer_store)

    #-----
    # test_beam
    #-----
    @utils.requires_ds(beam)
    def test_beam(self):
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
        ds = utils.data_dir_load(beam)
        assert_equal(str(ds), "plt03008")
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "c"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "Ex"]
        # Run the small_patch_amr test suite
        hashes = self.small_patch_amr(ds, _warpx_fields, weights, axes, ds_objs)
        # Save or compare answer
        utils.handle_hashes(self.save_dir, 'boxlib-test-beam', hashes, self.answer_store)

    #-----
    # test_warpx_particle_io
    #-----
    @requires_file(plasma)
    def test_warpx_particle_io(self, ds_plasma):
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
        ds = ds_plasma
        grid = ds.index.grids[0]
        # read directly from the header
        npart0_grid_0 = 344  
        npart1_grid_0 = 69632
        assert_equal(grid['particle0', 'particle_position_x'].size, npart0_grid_0)
        assert_equal(grid['particle1', 'particle_position_y'].size, npart1_grid_0)
        assert_equal(grid['all', 'particle_position_z'].size, npart0_grid_0 + npart1_grid_0)
        # read directly from the header
        npart0 = 1360  
        npart1 = 802816  
        ad = ds.all_data()
        assert_equal(ad['particle0', 'particle_velocity_x'].size, npart0)
        assert_equal(ad['particle1', 'particle_velocity_y'].size, npart1)
        assert_equal(ad['all', 'particle_velocity_z'].size, npart0 + npart1)
        np.all(ad['particle1', 'particle_mass'] == ad['particle1', 'particle_mass'][0])
        np.all(ad['particle0', 'particle_mass'] == ad['particle0', 'particle_mass'][0])
        left_edge = ds.arr([-7.5e-5, -7.5e-5, -7.5e-5], 'code_length')
        right_edge = ds.arr([2.5e-5, 2.5e-5, 2.5e-5], 'code_length')
        center = 0.5*(left_edge + right_edge)
        reg = ds.region(center, left_edge, right_edge)
        assert(np.all(np.logical_and(reg['particle_position_x'] <= right_edge[0], 
                                     reg['particle_position_x'] >= left_edge[0])))
        assert(np.all(np.logical_and(reg['particle_position_y'] <= right_edge[1], 
                                     reg['particle_position_y'] >= left_edge[1])))
        assert(np.all(np.logical_and(reg['particle_position_z'] <= right_edge[2], 
                                     reg['particle_position_z'] >= left_edge[2])))

    #-----
    # test_raw_fields
    #-----
    @utils.requires_ds(raw_fields)
    def test_raw_fields(self):
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
        gv_hd = b''
        ds_fn = utils.data_dir_load(raw_fields)
        for field in _raw_fields:
            gv_hd += self.grid_values_test(ds_fn, field)
        hd = {'grid_values' : utils.generate_hash(gv_hd)}
        utils.handle_hashes(self.save_dir, 'boxlib-test-raw-fields', hd, self.answer_store)

    #-----
    # test_OrionDataset
    #-----
    @requires_file(rt)
    def test_OrionDataset(self, ds_rt):
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
        assert isinstance(ds_rt, OrionDataset)

    #-----
    # teset_NyxDataset
    #-----
    @requires_file(LyA)
    def test_NyxDataset(self, ds_LyA):
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
        assert isinstance(ds_LyA, NyxDataset)

    #-----
    # test_CastroDataset
    #-----
    @requires_file(RT_particles)
    def test_CastroDataset(self, ds_RT_particles):
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
        assert isinstance(ds_RT_particles, CastroDataset)

    #-----
    # test_WarpXDataset
    #-----
    @requires_file(plasma)
    def test_WarpXDataset(self, ds_plasma):
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
        assert isinstance(ds_plasma, WarpXDataset)

    #-----
    # test_units_override
    #-----
    @requires_file(rt)
    def test_units_override(self, ds_rt):
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
        units_override_check(ds_rt, rt)

    #-----
    # test_nyx_no_part
    #-----
    @requires_file(nyx_no_particles)
    def test_nyx_no_part(self):
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
        ds = utils.data_dir_load(nyx_no_particles)
        assert isinstance(ds, NyxDataset)
        fields = sorted(
            [('boxlib', 'H'), ('boxlib', 'He'), ('boxlib', 'MachNumber'),
             ('boxlib', 'Ne'), ('boxlib', 'Rank'), ('boxlib', 'StateErr'),
             ('boxlib', 'Temp'), ('boxlib', 'X(H)'), ('boxlib', 'X(He)'),
             ('boxlib', 'density'), ('boxlib', 'divu'), ('boxlib', 'eint_E'),
             ('boxlib', 'eint_e'), ('boxlib', 'entropy'), ('boxlib', 'forcex'),
             ('boxlib', 'forcey'), ('boxlib', 'forcez'), ('boxlib', 'kineng'),
             ('boxlib', 'logden'), ('boxlib', 'magmom'), ('boxlib', 'magvel'),
             ('boxlib', 'magvort'), ('boxlib', 'pressure'), ('boxlib', 'rho_E'),
             ('boxlib', 'rho_H'), ('boxlib', 'rho_He'), ('boxlib', 'rho_e'),
             ('boxlib', 'soundspeed'), ('boxlib', 'x_velocity'), ('boxlib', 'xmom'),
             ('boxlib', 'y_velocity'), ('boxlib', 'ymom'), ('boxlib', 'z_velocity'),
             ('boxlib', 'zmom')])
        assert_equal(sorted(ds.field_list), fields)

    #-----
    # test_maestro_parameters
    #-----
    @requires_file(msubch)
    def test_maestro_parameters(self):
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
        ds = utils.data_dir_load(msubch)
        assert isinstance(ds, MaestroDataset)
        # Check a string parameter
        assert(ds.parameters['plot_base_name']=="subCh_hot_baserun_plt")
        assert(type(ds.parameters['plot_base_name']) is str)
        # Check boolean parameters: T or F
        assert(ds.parameters['use_thermal_diffusion'] is False)
        assert(type(ds.parameters['use_thermal_diffusion']) is bool)
        assert(ds.parameters['do_burning'] is True)
        assert(type(ds.parameters['do_burning']) is bool)
        # Check a float parameter with a decimal point
        assert(ds.parameters['sponge_kappa']==float('10.00000000'))
        assert(type(ds.parameters['sponge_kappa']) is float)
        # Check a float parameter with E exponent notation
        assert(ds.parameters['small_dt']==float('0.1000000000E-09'))
        # Check an int parameter
        assert(ds.parameters['s0_interp_type']==3)
        assert(type(ds.parameters['s0_interp_type']) is int)
