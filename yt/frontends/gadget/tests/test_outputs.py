"""
Title: test_gadget.py
Purpose: Gadget frontend tests
Notes:
Copyright (c) 2015, yt Development Team.
Distributed under the terms of the Modified BSD License.
The full license is in the file COPYING.txt, distributed with this
software.
"""
from collections import OrderedDict
from itertools import product
import os
import shutil
import tempfile

import pytest

import yt
from yt.testing import requires_file
from yt.frontends.gadget.api import GadgetHDF5Dataset, GadgetDataset
from yt.frontends.gadget.testing import fake_gadget_binary

import framework as fw
import utils


# Test data
isothermal_h5 = "IsothermalCollapse/snap_505.hdf5"
isothermal_bin = "IsothermalCollapse/snap_505"
BE_Gadget = "BigEndianGadgetBinary/BigEndianGadgetBinary"
LE_SnapFormat2 = "Gadget3-snap-format2/Gadget3-snap-format2"
keplerian_ring = "KeplerianRing/keplerian_ring_0020.hdf5"


# Globals
# This maps from field names to weight field names to use for projections
iso_fields = OrderedDict(
    [
        (("gas", "density"), None),
        (("gas", "temperature"), None),
        (("gas", "temperature"), ('gas', 'density')),
        (('gas', 'velocity_magnitude'), None),
        (("deposit", "all_density"), None),
        (("deposit", "all_count"), None),
        (("deposit", "all_cic"), None),
        (("deposit", "PartType0_density"), None),
    ]
)
iso_kwargs = dict(bounding_box=[[-3, 3], [-3, 3], [-3, 3]])


#============================================
#                 TestGadget
#============================================
class TestGadget(fw.AnswerTest):
    """
    Container for the gadget frontend answer tests.

    Attributes:
    -----------
        pass

    Methods:
    --------
        pass
    """
    #-----
    # test_gadget_binary
    #-----
    @pytest.mark.usefixtures('temp_dir')
    def test_gadget_binary(self):
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
        header_specs = ['default', 'default+pad32', ['default', 'pad32']]
        for header_spec, endian, fmt in product(header_specs, '<>', [1, 2]):
            fake_snap = fake_gadget_binary(
                header_spec=header_spec,
                endian=endian,
                fmt=fmt
            )
            ds = yt.load(fake_snap, header_spec=header_spec)
            assert isinstance(ds, GadgetDataset)
            ds.field_list

    #-----
    # test_gadget_hdf5
    #-----
    @requires_file(isothermal_h5)
    def test_gadget_hdf5(self):
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
        assert isinstance(utils.data_dir_load(isothermal_h5, kwargs=iso_kwargs),
                          GadgetHDF5Dataset)

    #-----
    # test_non_cosmo_dataset
    #-----
    @requires_file(keplerian_ring)
    def test_non_cosmo_dataset(self):
        """
        Non-cosmological datasets may not have the cosmological parametrs in the
        Header. The code should fall back gracefully when they are not present,
        with the Redshift set to 0.

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
        data = utils.data_dir_load(keplerian_ring)
        assert data.current_redshift == 0.0
        assert data.cosmological_simulation == 0

    #-----
    # test_iso_collapse
    #-----
    @utils.requires_ds(isothermal_h5)
    def test_iso_collapse(self):
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
        ds = utils.data_dir_load(isothermal_h5, kwargs=iso_kwargs)
        hashes = self.sph_answer(ds, 'snap_505', 2**17, iso_fields)
        utils.handle_hashes(self.save_dir, 'gadget-test-iso-collapse', hashes, self.answer_store)

    #-----
    # test_pid_uniqueness
    #-----
    @utils.requires_ds(LE_SnapFormat2)
    def test_pid_uniqueness(self):
        """
        ParticleIDs should be unique.

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
        ds = utils.data_dir_load(LE_SnapFormat2)
        ad = ds.all_data()
        pid = ad['ParticleIDs']
        assert len(pid) == len(set(pid.v))

    #-----
    # test_bigendian_field_access
    #-----
    @utils.requires_ds(BE_Gadget)
    def test_bigendian_field_access(self):
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
        ds = utils.data_dir_load(BE_Gadget)
        data = ds.all_data()
        data['Halo', 'Velocities']
