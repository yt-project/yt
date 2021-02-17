"""
AdaptaHOP frontend tests



"""

import numpy as np

from yt.frontends.adaptahop.data_structures import AdaptaHOPDataset
from yt.testing import requires_file
from yt.utilities.answer_testing.framework import data_dir_load

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

    skip_list = ("particle_identities", "mesh_id")
    fields = [
        (ptype, field)
        for (ptype, field) in ds.derived_field_list
        if (ptype == "halos") and (field not in skip_list)
    ]

    ad = ds.all_data()

    def test_access(ptype, field):
        ad[ptype, field]

    for (ptype, field) in fields:
        test_access(ptype, field)


@requires_file(r0)
@requires_file(r1)
def test_get_halo():
    ds_parent = data_dir_load(r0)
    ds = data_dir_load(r1, kwargs=dict(parent_ds=ds_parent))

    halo = ds.halo(1, ptype="io")

    # Check halo objet has position, velocity, mass and members attributes
    for attr_name in ("mass", "position", "velocity", "member_ids"):
        getattr(halo, attr_name)

    members = np.sort(halo.member_ids)

    # Check sphere contains all the members
    id_sphere = halo.sphere["io", "particle_identity"].astype(int)
    assert len(np.lib.arraysetops.setdiff1d(members, id_sphere)) == 0

    # Check region contains *only* halo particles
    id_reg = np.sort(halo["io", "particle_identity"].astype(int))

    np.testing.assert_equal(members, id_reg)
