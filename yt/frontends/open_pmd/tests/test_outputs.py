"""
Title: test_open_pmd.py
Purpose: openPMD frontend tests
Notes:
    Copyright (c) 2016, Fabian Koller (HZDR).
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import numpy as np
import pytest

from yt.convenience import load
from yt.frontends.open_pmd.data_structures import OpenPMDDataset
from yt.testing import \
    assert_almost_equal, \
    assert_array_equal, \
    assert_equal, \
    requires_file
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils

# Test data
twoD = "example-2d/hdf5/data00000100.h5"
threeD = "example-3d/hdf5/data00000100.h5"
noFields = "no_fields/data00000400.h5"
noParticles = "no_particles/data00000400.h5"
groupBased = "singleParticle/simData.h5"


#============================================
#                TestOpenPMD
#============================================
@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
class TestOpenPMD(fw.AnswerTest):
    """
    Container for open_pmd answer tests.

    Attributes:
    -----------
        pass

    Methods:
    --------
        pass
    """
    #-----
    # test_3d_out
    #-----
    @requires_file(threeD)
    def test_3d_out(self, ds_threeD):
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
        field_list = [('all', 'particle_charge'),
                      ('all', 'particle_mass'),
                      ('all', 'particle_momentum_x'),
                      ('all', 'particle_momentum_y'),
                      ('all', 'particle_momentum_z'),
                      ('all', 'particle_positionCoarse_x'),
                      ('all', 'particle_positionCoarse_y'),
                      ('all', 'particle_positionCoarse_z'),
                      ('all', 'particle_positionOffset_x'),
                      ('all', 'particle_positionOffset_y'),
                      ('all', 'particle_positionOffset_z'),
                      ('all', 'particle_weighting'),
                      ('io', 'particle_charge'),
                      ('io', 'particle_mass'),
                      ('io', 'particle_momentum_x'),
                      ('io', 'particle_momentum_y'),
                      ('io', 'particle_momentum_z'),
                      ('io', 'particle_positionCoarse_x'),
                      ('io', 'particle_positionCoarse_y'),
                      ('io', 'particle_positionCoarse_z'),
                      ('io', 'particle_positionOffset_x'),
                      ('io', 'particle_positionOffset_y'),
                      ('io', 'particle_positionOffset_z'),
                      ('io', 'particle_weighting'),
                      ('openPMD', 'E_x'),
                      ('openPMD', 'E_y'),
                      ('openPMD', 'E_z'),
                      ('openPMD', 'rho')]
        domain_dimensions = [26, 26, 201] * np.ones_like(ds_threeD.domain_dimensions)
        domain_width = [2.08e-05, 2.08e-05, 2.01e-05] * \
            np.ones_like(ds_threeD.domain_left_edge)
        assert isinstance(ds_threeD, OpenPMDDataset)
        assert_equal(ds_threeD.dimensionality, 3)
        assert_equal(ds_threeD.particle_types_raw, ('io',))
        assert "all" in ds_threeD.particle_unions
        assert_array_equal(ds_threeD.field_list, field_list)
        assert_array_equal(ds_threeD.domain_dimensions, domain_dimensions)
        assert_almost_equal(ds_threeD.current_time,
            3.28471214521e-14 * np.ones_like(ds_threeD.current_time))
        assert_almost_equal(ds_threeD.domain_right_edge - \
            ds_threeD.domain_left_edge, domain_width)

    #-----
    # test_2d_out
    #-----
    @requires_file(twoD)
    def test_2d_out(self, ds_twoD):
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
        field_list = [('Hydrogen1+', 'particle_charge'),
                      ('Hydrogen1+', 'particle_mass'),
                      ('Hydrogen1+', 'particle_momentum_x'),
                      ('Hydrogen1+', 'particle_momentum_y'),
                      ('Hydrogen1+', 'particle_momentum_z'),
                      ('Hydrogen1+', 'particle_positionCoarse_x'),
                      ('Hydrogen1+', 'particle_positionCoarse_y'),
                      ('Hydrogen1+', 'particle_positionCoarse_z'),
                      ('Hydrogen1+', 'particle_positionOffset_x'),
                      ('Hydrogen1+', 'particle_positionOffset_y'),
                      ('Hydrogen1+', 'particle_positionOffset_z'),
                      ('Hydrogen1+', 'particle_weighting'),
                      ('all', 'particle_charge'),
                      ('all', 'particle_mass'),
                      ('all', 'particle_momentum_x'),
                      ('all', 'particle_momentum_y'),
                      ('all', 'particle_momentum_z'),
                      ('all', 'particle_positionCoarse_x'),
                      ('all', 'particle_positionCoarse_y'),
                      ('all', 'particle_positionCoarse_z'),
                      ('all', 'particle_positionOffset_x'),
                      ('all', 'particle_positionOffset_y'),
                      ('all', 'particle_positionOffset_z'),
                      ('all', 'particle_weighting'),
                      ('electrons', 'particle_charge'),
                      ('electrons', 'particle_mass'),
                      ('electrons', 'particle_momentum_x'),
                      ('electrons', 'particle_momentum_y'),
                      ('electrons', 'particle_momentum_z'),
                      ('electrons', 'particle_positionCoarse_x'),
                      ('electrons', 'particle_positionCoarse_y'),
                      ('electrons', 'particle_positionCoarse_z'),
                      ('electrons', 'particle_positionOffset_x'),
                      ('electrons', 'particle_positionOffset_y'),
                      ('electrons', 'particle_positionOffset_z'),
                      ('electrons', 'particle_weighting'),
                      ('openPMD', 'B_x'),
                      ('openPMD', 'B_y'),
                      ('openPMD', 'B_z'),
                      ('openPMD', 'E_x'),
                      ('openPMD', 'E_y'),
                      ('openPMD', 'E_z'),
                      ('openPMD', 'J_x'),
                      ('openPMD', 'J_y'),
                      ('openPMD', 'J_z'),
                      ('openPMD', 'rho')]
        domain_dimensions = [51, 201, 1] * np.ones_like(ds_twoD.domain_dimensions)
        domain_width = [3.06e-05, 2.01e-05, 1e+0] * \
            np.ones_like(ds_twoD.domain_left_edge)
        assert isinstance(ds_twoD, OpenPMDDataset)
        assert_equal(ds_twoD.dimensionality, 2)
        assert_equal(ds_twoD.particle_types_raw, ('Hydrogen1+', 'electrons'))
        assert "all" in ds_twoD.particle_unions
        assert_array_equal(ds_twoD.field_list, field_list)
        assert_array_equal(ds_twoD.domain_dimensions, domain_dimensions)
        assert_almost_equal(ds_twoD.current_time,
            3.29025596712e-14 * np.ones_like(ds_twoD.current_time))
        assert_almost_equal(ds_twoD.domain_right_edge - \
            ds_twoD.domain_left_edge, domain_width)

    #-----
    # test_no_fields_out
    #-----
    @requires_file(noFields)
    def test_no_fields_out(self, ds_noFields):
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
        field_list = [('all', 'particle_charge'),
                      ('all', 'particle_id'),
                      ('all', 'particle_mass'),
                      ('all', 'particle_momentum_x'),
                      ('all', 'particle_momentum_y'),
                      ('all', 'particle_momentum_z'),
                      ('all', 'particle_positionCoarse_x'),
                      ('all', 'particle_positionCoarse_y'),
                      ('all', 'particle_positionCoarse_z'),
                      ('all', 'particle_positionOffset_x'),
                      ('all', 'particle_positionOffset_y'),
                      ('all', 'particle_positionOffset_z'),
                      ('all', 'particle_weighting'),
                      ('io', 'particle_charge'),
                      ('io', 'particle_id'),
                      ('io', 'particle_mass'),
                      ('io', 'particle_momentum_x'),
                      ('io', 'particle_momentum_y'),
                      ('io', 'particle_momentum_z'),
                      ('io', 'particle_positionCoarse_x'),
                      ('io', 'particle_positionCoarse_y'),
                      ('io', 'particle_positionCoarse_z'),
                      ('io', 'particle_positionOffset_x'),
                      ('io', 'particle_positionOffset_y'),
                      ('io', 'particle_positionOffset_z'),
                      ('io', 'particle_weighting')]
        domain_dimensions = [1, 1, 1] * np.ones_like(ds_noFields.domain_dimensions)
        domain_width = [1, 1, 1] * np.ones_like(ds_noFields.domain_left_edge)
        assert isinstance(ds_noFields, OpenPMDDataset)
        assert_equal(ds_noFields.dimensionality, 3)
        assert_equal(ds_noFields.particle_types_raw, ('io', ))
        assert "all" in ds_noFields.particle_unions
        assert_array_equal(ds_noFields.field_list, field_list)
        assert_array_equal(ds_noFields.domain_dimensions, domain_dimensions)
        assert_almost_equal(ds_noFields.current_time,
            1.3161023868481013e-13 * np.ones_like(ds_noFields.current_time))
        assert_almost_equal(ds_noFields.domain_right_edge - \
            ds_noFields.domain_left_edge, domain_width)

    #-----
    # test_no_particles_out
    #-----
    @requires_file(noParticles)
    def test_no_particles_out(self, ds_noParticles):
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
        field_list = [('openPMD', 'E_x'),
                      ('openPMD', 'E_y'),
                      ('openPMD', 'E_z'),
                      ('openPMD', 'rho')]
        domain_dimensions = [51, 201, 1] * np.ones_like(ds_noParticles.domain_dimensions)
        domain_width = [3.06e-05, 2.01e-05, 1e+0] * \
            np.ones_like(ds_noParticles.domain_left_edge)
        assert isinstance(ds_noParticles, OpenPMDDataset)
        assert_equal(ds_noParticles.dimensionality, 2)
        assert_equal(ds_noParticles.particle_types_raw, ('io', ))
        assert "all" not in ds_noParticles.particle_unions
        assert_array_equal(ds_noParticles.field_list, field_list)
        assert_array_equal(ds_noParticles.domain_dimensions, domain_dimensions)
        assert_almost_equal(ds_noParticles.current_time,
            1.3161023868481013e-13 * np.ones_like(ds_noParticles.current_time))
        assert_almost_equal(ds_noParticles.domain_right_edge - \
            ds_noParticles.domain_left_edge, domain_width)

    #-----
    # test_groupBased_out
    #-----
    @requires_file(groupBased)
    def test_groupBased_out(self):
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
        dss = load(groupBased)
        field_list = [('all', 'particle_charge'),
                      ('all', 'particle_mass'),
                      ('all', 'particle_momentum_x'),
                      ('all', 'particle_momentum_y'),
                      ('all', 'particle_momentum_z'),
                      ('all', 'particle_positionCoarse_x'),
                      ('all', 'particle_positionCoarse_y'),
                      ('all', 'particle_positionCoarse_z'),
                      ('all', 'particle_positionOffset_x'),
                      ('all', 'particle_positionOffset_y'),
                      ('all', 'particle_positionOffset_z'),
                      ('all', 'particle_weighting'),
                      ('io', 'particle_charge'),
                      ('io', 'particle_mass'),
                      ('io', 'particle_momentum_x'),
                      ('io', 'particle_momentum_y'),
                      ('io', 'particle_momentum_z'),
                      ('io', 'particle_positionCoarse_x'),
                      ('io', 'particle_positionCoarse_y'),
                      ('io', 'particle_positionCoarse_z'),
                      ('io', 'particle_positionOffset_x'),
                      ('io', 'particle_positionOffset_y'),
                      ('io', 'particle_positionOffset_z'),
                      ('io', 'particle_weighting'),
                      ('openPMD', 'J_x'),
                      ('openPMD', 'J_y'),
                      ('openPMD', 'J_z'),
                      ('openPMD', 'e-chargeDensity')]
        domain_dimensions = [32, 64, 64] *\
            np.ones_like(dss[0].domain_dimensions)
        domain_width = [0.0002752, 0.0005504, 0.0005504] * \
            np.ones_like(dss[0].domain_left_edge)
        assert_equal(len(dss), 101)
        for ds in dss:
            assert_equal(str(ds), "simData.h5")
            assert_equal(ds.dimensionality, 3)
            assert_equal(ds.particle_types_raw, ('io', ))
            assert_array_equal(ds.field_list, field_list)
            assert_array_equal(ds.domain_dimensions,
                domain_dimensions)
            assert ds.current_time >= np.zeros_like(ds.current_time)
            assert ds.current_time <= 1.6499999999999998e-12 * \
                np.ones_like(ds.current_time)
            assert_almost_equal(ds.domain_right_edge - \
                ds.domain_left_edge, domain_width)
