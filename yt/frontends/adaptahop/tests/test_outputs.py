"""
AdaptaHOP frontend tests



"""

import numpy as np
import pytest

from yt.frontends.adaptahop.data_structures import AdaptaHOPDataset
from yt.testing import requires_file
from yt.utilities.answer_testing.framework import data_dir_load

r0 = "output_00080/info_00080.txt"
brick_old = "output_00080_halos/tree_bricks080"
brick_new = "output_00080_new_halos/tree_bricks080"


@requires_file(r0)
@requires_file(brick_old)
@requires_file(brick_new)
@pytest.mark.parametrize("brick", [brick_old, brick_new])
def test_opening(brick):
    ds_parent = data_dir_load(r0)
    ds = data_dir_load(brick, kwargs=dict(parent_ds=ds_parent))

    ds.index

    assert isinstance(ds, AdaptaHOPDataset)


@requires_file(r0)
@requires_file(brick_old)
@requires_file(brick_new)
@pytest.mark.parametrize("brick", [brick_old, brick_new])
def test_field_access(brick):
    ds_parent = data_dir_load(r0)
    ds = data_dir_load(brick, kwargs=dict(parent_ds=ds_parent))

    skip_list = ("particle_identities", "mesh_id", "radius_profile", "rho_profile")
    fields = [
        (ptype, field)
        for (ptype, field) in ds.derived_field_list
        if (ptype == "halos") and (field not in skip_list)
    ]

    ad = ds.all_data()

    for (ptype, field) in fields:
        ad[(ptype, field)]


@requires_file(r0)
@requires_file(brick_old)
@requires_file(brick_new)
@pytest.mark.parametrize("brick", [brick_old, brick_new])
def test_get_halo(brick):
    ds_parent = data_dir_load(r0)
    ds = data_dir_load(brick, kwargs=dict(parent_ds=ds_parent))

    halo = ds.halo(1, ptype="io")

    # Check halo object has position, velocity, mass and members attributes
    for attr_name in ("mass", "position", "velocity", "member_ids"):
        getattr(halo, attr_name)

    members = np.sort(halo.member_ids)

    # Check sphere contains all the members
    id_sphere = halo.sphere["io", "particle_identity"].astype(int)
    assert len(np.lib.arraysetops.setdiff1d(members, id_sphere)) == 0

    # Check region contains *only* halo particles
    id_reg = np.sort(halo["io", "particle_identity"].astype(int))

    np.testing.assert_equal(members, id_reg)
