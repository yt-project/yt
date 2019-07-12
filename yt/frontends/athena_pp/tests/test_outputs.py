"""
title: test_athena_pp.py
Purpose: Athena++ frontend tests
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import numpy as np

from yt.convenience import load
from yt.frontends.athena_pp.api import AthenaPPDataset
from yt.testing import \
    units_override_check, \
    assert_allclose
from yt.utilities.answer_testing import \
    framework as fw, \
    utils


# Test data
disk = "KeplerianDisk/disk.out1.00000.athdf"
AM06 = "AM06/AM06.out1.00400.athdf"


#============================================
#                TestAthenaPP
#============================================
class TestAthenaPP(fw.AnswerTest):
    """
    Container for athena++ frontent tests.

    Attributes:
    -----------
        pass

    Methods:
    --------
        pass
    """
    #-----
    # test_disk
    #-----
    def test_disk(self):
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
        fields = ("density", "velocity_r")
        ds = utils.data_dir_load(disk)
        assert str(ds) == "disk.out1.00000"
        dd = ds.all_data()
        vol = (ds.domain_right_edge[0]**3-ds.domain_left_edge[0]**3)/3.0
        vol *= np.cos(ds.domain_left_edge[1])-np.cos(ds.domain_right_edge[1])
        vol *= ds.domain_right_edge[2].v-ds.domain_left_edge[2].v
        assert_allclose(dd.quantities.total_quantity("cell_volume"), vol)
        ga_hd = b''
        for field in fields:
            def field_func(name):
                return dd[field]
            ga_hd += self.generic_array_test(ds, field_func, args=[field])
        hashes = {}
        hashes['generic_array'] = utils.generate_hash(ga_hd)
        utils.handle_hashes('athenapp_disk', hashes, self.answer_store)

    #-----
    # test_AM06
    #-----
    def test_AM06(self):
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
        # Arrays for testing
        axes = [0, 1, 2]
        center = "max"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "density"]
        fields = ("temperature",
            "density",
            "velocity_magnitude",
            "magnetic_field_x"
        )
        ds = utils.data_dir_load(AM06)
        assert str(ds) == "AM06.out1.00400"
        # Run the small_patch_amr test suite
        hashes = self.small_patch_amr(ds, fields, weights, axes, ds_objs)
        # Save or compare answer
        utils.handle_hashes('athenapp_AM06', hashes, self.answer_store)

    #-----
    # test_AM06_override
    #-----
    def test_AM06_override(self):
        """
        Verify that overriding units causes derived unit values to be
        updated. See issue #1259.

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
        uo_AM06 = {
            'length_unit': (1.0, 'kpc'),
            'mass_unit': (1.0, 'Msun'),
            'time_unit': (1.0, 'Myr'),
        }
        ds = load(AM06, units_override=uo_AM06)
        assert float(ds.magnetic_unit.in_units('gauss')) == 9.01735778342523e-08

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
        units_override_check(AM06)

    #-----
    # test_AthenaPPDataset
    #-----
    def test_AthenaPPDataset(self):
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
        assert isinstance(utils.data_dir_load(AM06), AthenaPPDataset)
