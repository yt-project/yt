"""
Gadget frontend tests using the IsothermalCollapse dataset




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
    requires_ds, \
    data_dir_load, \
    sph_answer
from yt.frontends.gadget.api import GadgetHDF5Dataset


isothermal_h5 = "IsothermalCollapse/snap_505.hdf5"
isothermal_bin = "IsothermalCollapse/snap_505"
gdg = "GadgetDiskGalaxy/snapshot_200.hdf5"

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

gdg_fields = iso_fields.copy()
gdg_fields["deposit", "PartType4_density"] = None
gdg_kwargs = dict(bounding_box=[[-1e5, 1e5], [-1e5, 1e5], [-1e5, 1e5]])


@requires_file(isothermal_h5)
@requires_file(isothermal_bin)
def test_GadgetDataset():
    kwargs = dict(bounding_box=[[-3,3], [-3,3], [-3,3]])
    assert isinstance(data_dir_load(isothermal, kwargs=kwargs),
                      GadgetHDF5Dataset)
    assert isinstance(data_dir_load(isothermal_bin, kwargs=iso_kwargs),
                      GadgetDataset)


@requires_ds(isothermal_h5)
def test_iso_collapse():
    for test in sph_answer(isothermal_h5, 'snap_505', 2**17,
                           iso_fields, ds_kwargs=iso_kwargs):
        yield test

@requires_ds(gdg, big_data=True)
def test_gadget_disk_galaxy():
    for test in sph_answer(gdg, 'snapshot_200', 11907080, gdg_fields,
                           ds_kwargs=gdg_kwargs):
        yield test
