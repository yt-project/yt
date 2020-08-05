"""
Title: test_gadget.py
Purpose: Gadget frontend tests
Notes:
    Copyright (c) 2015, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
from itertools import product

import pytest

import yt
from yt.frontends.gadget.api import GadgetDataset, GadgetHDF5Dataset
from yt.frontends.gadget.testing import fake_gadget_binary
from yt.testing import ParticleSelectionComparison, requires_file
from yt.utilities.answer_testing.answer_tests import sph_answer
from yt.utilities.answer_testing.utils import data_dir_load, requires_ds

# Test data
isothermal_h5 = "IsothermalCollapse/snap_505.hdf5"
isothermal_bin = "IsothermalCollapse/snap_505"
BE_Gadget = "BigEndianGadgetBinary/BigEndianGadgetBinary"
LE_SnapFormat2 = "Gadget3-snap-format2/Gadget3-snap-format2"
keplerian_ring = "KeplerianRing/keplerian_ring_0020.hdf5"
snap_33 = "snapshot_033/snap_033.0.hdf5"
snap_33_dir = "snapshot_033/"


iso_kwargs = dict(bounding_box=[[-3, 3], [-3, 3], [-3, 3]])


@pytest.mark.answer_test
@pytest.mark.usefixtures("answer_file")
class TestGadget:
    @pytest.mark.usefixtures("temp_dir")
    def test_gadget_binary(self):
        header_specs = ["default", "default+pad32", ["default", "pad32"]]
        for header_spec, endian, fmt in product(header_specs, "<>", [1, 2]):
            fake_snap = fake_gadget_binary(
                header_spec=header_spec, endian=endian, fmt=fmt
            )
            ds = yt.load(fake_snap, header_spec=header_spec)
            assert isinstance(ds, GadgetDataset)
            ds.field_list

    @requires_file(isothermal_h5)
    def test_gadget_hdf5(self):
        assert isinstance(
            data_dir_load(isothermal_h5, kwargs=iso_kwargs), GadgetHDF5Dataset
        )

    @requires_file(keplerian_ring)
    def test_non_cosmo_dataset(self):
        r"""Non-cosmological datasets may not have the cosmological parameters in the
        Header. The code should fall back gracefully when they are not present,
        with the Redshift set to 0.
        """
        data = data_dir_load(keplerian_ring)
        assert data.current_redshift == 0.0
        assert data.cosmological_simulation == 0

    @pytest.mark.usefixtures("hashing")
    @requires_ds(isothermal_h5)
    def test_iso_collapse(self, f, w, d, a):
        ds = data_dir_load(isothermal_h5, kwargs=iso_kwargs)
        self.hashes.update(sph_answer(ds, "snap_505", 2 ** 17, f, w, d, a))

    @requires_ds(LE_SnapFormat2)
    def test_pid_uniqueness(self):
        ds = data_dir_load(LE_SnapFormat2)
        ad = ds.all_data()
        pid = ad["ParticleIDs"]
        assert len(pid) == len(set(pid.v))

    @requires_file(snap_33)
    @requires_file(snap_33_dir)
    def test_multifile_read(self):
        """
        Tests to make sure multi-file gadget snapshot can be loaded by passing '.0' file
        or by passing the directory containing the multi-file snapshot.
        """
        assert isinstance(data_dir_load(snap_33), GadgetDataset)
        assert isinstance(data_dir_load(snap_33_dir), GadgetDataset)

    @requires_file(snap_33)
    def test_particle_subselection(self):
        # This checks that we correctly subselect from a dataset, first by making
        # sure we get all the particles, then by comparing manual selections against
        # them.
        ds = data_dir_load(snap_33)
        psc = ParticleSelectionComparison(ds)
        psc.run_defaults()

    @requires_ds(BE_Gadget)
    def test_bigendian_field_access(self):
        ds = data_dir_load(BE_Gadget)
        data = ds.all_data()
        data["Halo", "Velocities"]
