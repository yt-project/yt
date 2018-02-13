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

from yt.testing import requires_file
from yt.utilities.answer_testing.framework import \
    data_dir_load, \
    requires_ds, \
    sph_answer
from yt.frontends.gadget.api import GadgetHDF5Dataset, GadgetDataset

isothermal_h5 = "IsothermalCollapse/snap_505.hdf5"
isothermal_bin = "IsothermalCollapse/snap_505"
BE_Gadget = "BigEndianGadgetBinary/BigEndianGadgetBinary"
LE_SnapFormat2 = "Gadget3-snap-format2/Gadget3-snap-format2"
keplerian_ring = "KeplerianRing/keplerian_ring_test_data.hdf5"

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


@requires_file(isothermal_h5)
@requires_file(isothermal_bin)
@requires_file(BE_Gadget)
@requires_file(LE_SnapFormat2)
def test_gadget_dataset():
    assert isinstance(data_dir_load(isothermal_h5, kwargs=iso_kwargs),
                      GadgetHDF5Dataset)
    assert isinstance(data_dir_load(isothermal_bin, kwargs=iso_kwargs),
                      GadgetDataset)
    assert isinstance(data_dir_load(BE_Gadget), GadgetDataset)
    assert isinstance(data_dir_load(LE_SnapFormat2), GadgetDataset)


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
