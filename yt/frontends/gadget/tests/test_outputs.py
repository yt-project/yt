"""
Gadget frontend tests




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from collections import OrderedDict
from itertools import product
import os

import yt
from yt.testing import requires_file
from yt.utilities.answer_testing.framework import \
    data_dir_load, \
    requires_ds, \
    sph_answer
from yt.frontends.gadget.api import GadgetHDF5Dataset, GadgetDataset
from yt.frontends.gadget.testing import fake_gadget_binary

isothermal_h5 = "IsothermalCollapse/snap_505.hdf5"
isothermal_bin = "IsothermalCollapse/snap_505"
BE_Gadget = "BigEndianGadgetBinary/BigEndianGadgetBinary"
LE_SnapFormat2 = "Gadget3-snap-format2/Gadget3-snap-format2"
keplerian_ring = "KeplerianRing/keplerian_ring_0020.hdf5"

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


def test_gadget_binary():
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
        os.remove(fake_snap)


@requires_file(isothermal_h5)
def test_gadget_hdf5():
    assert isinstance(data_dir_load(isothermal_h5, kwargs=iso_kwargs),
                      GadgetHDF5Dataset)


@requires_file(keplerian_ring)
def test_non_cosmo_dataset():
    """
    Non-cosmological datasets may not have the cosmological parametrs in the
    Header. The code should fall back gracefully when they are not present,
    with the Redshift set to 0.
    """
    data = data_dir_load(keplerian_ring)
    assert data.current_redshift == 0.0
    assert data.cosmological_simulation == 0


@requires_ds(isothermal_h5)
def test_iso_collapse():
    ds = data_dir_load(isothermal_h5, kwargs=iso_kwargs)
    for test in sph_answer(ds, 'snap_505', 2**17, iso_fields):
        test_iso_collapse.__name__ = test.description
        yield test


@requires_ds(LE_SnapFormat2)
def test_pid_uniqueness():
    """
    ParticleIDs should be unique.
    """
    ds = data_dir_load(LE_SnapFormat2)
    ad = ds.all_data()
    pid = ad['ParticleIDs']
    assert len(pid) == len(set(pid.v))
