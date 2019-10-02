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

import pytest

import yt
from yt.testing import requires_file
from yt.frontends.gadget.api import GadgetHDF5Dataset, GadgetDataset
from yt.frontends.gadget.testing import fake_gadget_binary
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils


# Test data
isothermal_h5 = "IsothermalCollapse/snap_505.hdf5"
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
@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason=("--with-answer-testing not set."))
@pytest.mark.usefixtures('answer_file')
class TestGadget(fw.AnswerTest):
    #-----
    # test_gadget_binary
    #-----
    @pytest.mark.usefixtures('temp_dir')
    def test_gadget_binary(self):
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
    def test_gadget_hdf5(self, ds_isothermal_h5):
        assert isinstance(ds_isothermal_h5, GadgetHDF5Dataset)

    #-----
    # test_non_cosmo_dataset
    #-----
    @requires_file(keplerian_ring)
    def test_non_cosmo_dataset(self):
        r"""Non-cosmological datasets may not have the cosmological parametrs in the
        Header. The code should fall back gracefully when they are not present,
        with the Redshift set to 0.
        """
        data = utils.data_dir_load(keplerian_ring)
        assert data.current_redshift == 0.0
        assert data.cosmological_simulation == 0

    #-----
    # test_iso_collapse
    #-----
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(isothermal_h5)
    def test_iso_collapse(self, ds_isothermal_h5):
        self.hashes = self.sph_answer(ds_isothermal_h5, 'snap_505', 2**17, iso_fields)

    #-----
    # test_pid_uniqueness
    #-----
    @utils.requires_ds(LE_SnapFormat2)
    def test_pid_uniqueness(self):
        ds = utils.data_dir_load(LE_SnapFormat2)
        ad = ds.all_data()
        pid = ad['ParticleIDs']
        assert len(pid) == len(set(pid.v))

    #-----
    # test_bigendian_field_access
    #-----
    @utils.requires_ds(BE_Gadget)
    def test_bigendian_field_access(self):
        ds = utils.data_dir_load(BE_Gadget)
        data = ds.all_data()
        data['Halo', 'Velocities']
