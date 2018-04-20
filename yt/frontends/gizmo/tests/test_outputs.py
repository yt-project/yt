"""
Gizmo frontend tests




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from collections import OrderedDict

from yt.utilities.answer_testing.framework import \
    data_dir_load, \
    requires_ds, \
    sph_answer
from yt.frontends.gizmo.api import GizmoDataset

FIRE_m12i = 'FIRE_M12i_ref11/snapshot_600.hdf5'

# This maps from field names to weight field names to use for projections
fields = OrderedDict(
    [
        (("gas", "density"), None),
        (("gas", "temperature"), ('gas', 'density')),
        (("gas", "metallicity"), ('gas', 'density')),
        (("gas", "O_metallicity"), ('gas', 'density')),
        (('gas', 'velocity_magnitude'), None),
        (("deposit", "all_count"), None),
        (("deposit", "all_cic"), None),
        (("deposit", "PartType0_density"), None),
    ]
)

g64 = "gizmo_64/output/snap_N64L16_135.hdf5"
gmhd = "gizmo_mhd_mwdisk/gizmo_mhd_mwdisk.hdf5"


@requires_ds(g64, big_data=True)
def test_gizmo_64():
    ds = data_dir_load(g64)
    assert isinstance(ds, GizmoDataset)
    for test in sph_answer(ds, 'snap_N64L16_135', 524288, fields):
        test_gizmo_64.__name__ = test.description
        yield test


@requires_ds(gmhd, big_data=True)
def test_gizmo_mhd():
    """
    Magnetic fields should be loaded correctly when they are present.
    """
    kwargs = {
        'bounding_box': [[-400, 400]] * 3,
        'unit_system': 'code'
    }
    ds = data_dir_load(gmhd, kwargs=kwargs)
    ad = ds.all_data()
    ptype = 'PartType0'

    # Test vector magnetic field
    fmag = 'particle_magnetic_field'
    f = ad[ptype, fmag]
    assert str(f.units) == 'code_magnetic'
    assert f.shape == (409013, 3)

    # Test component magnetic fields
    for axis in 'xyz':
        f = ad[ptype, fmag + '_' + axis]
        assert str(f.units) == 'code_magnetic'
        assert f.shape == (409013,)
