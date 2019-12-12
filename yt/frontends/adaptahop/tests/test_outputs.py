"""
AdaptaHOP frontend tests



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.testing import requires_file
from yt.utilities.answer_testing.framework import \
    data_dir_load

from yt.frontends.adaptahop.data_structures import AdaptaHOPDataset

r0 = "output_00080/info_00080.txt"
r1 = "output_00080_halos/tree_bricks080"

@requires_file(r0)
@requires_file(r1)
def test_opening():
    ds_parent = data_dir_load(r0)
    ds = data_dir_load(r1, kwargs=dict(parent_ds=ds_parent))

    ds.index

    assert isinstance(ds, AdaptaHOPDataset)

@requires_file(r0)
@requires_file(r1)
def test_field_access():
    ds_parent = data_dir_load(r0)
    ds = data_dir_load(r1, kwargs=dict(parent_ds=ds_parent))

    fields = [(ptype, field) for (ptype, field) in ds.derived_field_list
              if (ptype == 'halos') and (field != 'particle_identities')]

    ad = ds.all_data()
    def test_access(ptype, field):
        ad[ptype, field]
    for (ptype, field) in fields:
        yield test_access, ptype, field
