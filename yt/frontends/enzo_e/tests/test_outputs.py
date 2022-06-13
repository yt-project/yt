import numpy as np

from yt.frontends.enzo_e.api import EnzoEDataset
from yt.frontends.enzo_e.fields import NODAL_FLAGS
from yt.testing import assert_array_equal, assert_equal, requires_file
from yt.utilities.answer_testing.framework import (
    FieldValuesTest,
    PixelizedProjectionValuesTest,
    create_obj,
    data_dir_load,
    requires_ds,
)
from yt.utilities.on_demand_imports import _h5py as h5py

_fields = (
    ("gas", "density"),
    ("gas", "specific_total_energy"),
    ("gas", "velocity_x"),
    ("gas", "velocity_y"),
)

_pfields = (
    ("all", "particle_position_x"),
    ("all", "particle_position_y"),
    ("all", "particle_position_z"),
    ("all", "particle_velocity_x"),
    ("all", "particle_velocity_y"),
    ("all", "particle_velocity_z"),
)

hello_world = "hello-0210/hello-0210.block_list"
ep_cosmo = "ENZOP_DD0140/ENZOP_DD0140.block_list"
orszag_tang = "ENZOE_orszag-tang_0.5/ENZOE_orszag-tang_0.5.block_list"


@requires_file(hello_world)
def test_EnzoEDataset():
    assert isinstance(data_dir_load(hello_world), EnzoEDataset)


@requires_ds(hello_world)
def test_hello_world():
    ds = data_dir_load(hello_world)

    dso = [None, ("sphere", ("max", (0.25, "unitary")))]
    for dobj_name in dso:
        for field in _fields:
            for axis in [0, 1, 2]:
                for weight_field in [None, ("gas", "density")]:
                    yield PixelizedProjectionValuesTest(
                        hello_world, axis, field, weight_field, dobj_name
                    )
            yield FieldValuesTest(hello_world, field, dobj_name)
        dobj = create_obj(ds, dobj_name)
        s1 = dobj[("index", "ones")].sum()
        s2 = sum(mask.sum() for block, mask in dobj.blocks)
        assert_equal(s1, s2)


@requires_ds(ep_cosmo)
def test_particle_fields():
    ds = data_dir_load(ep_cosmo)

    dso = [None, ("sphere", ("max", (0.1, "unitary")))]
    for dobj_name in dso:
        for field in _pfields:
            yield FieldValuesTest(ep_cosmo, field, dobj_name, particle_type=True)
        dobj = create_obj(ds, dobj_name)
        s1 = dobj[("index", "ones")].sum()
        s2 = sum(mask.sum() for block, mask in dobj.blocks)
        assert_equal(s1, s2)


@requires_file(hello_world)
def test_hierarchy():
    ds = data_dir_load(hello_world)

    fh = h5py.File(ds.index.grids[0].filename, mode="r")
    for grid in ds.index.grids:
        assert_array_equal(
            grid.LeftEdge.d, fh[grid.block_name].attrs["enzo_GridLeftEdge"]
        )
        assert_array_equal(ds.index.grid_left_edge[grid.id], grid.LeftEdge)
        assert_array_equal(ds.index.grid_right_edge[grid.id], grid.RightEdge)
        for child in grid.Children:
            assert (child.LeftEdge >= grid.LeftEdge).all()
            assert (child.RightEdge <= grid.RightEdge).all()
            assert_equal(child.Parent.id, grid.id)
    fh.close()


@requires_file(ep_cosmo)
def test_critical_density():
    ds = data_dir_load(ep_cosmo)

    c1 = (
        (ds.r["dark", "particle_mass"].sum() + ds.r["gas", "cell_mass"].sum())
        / ds.domain_width.prod()
        / ds.critical_density
    )
    c2 = (
        ds.omega_matter
        * (1 + ds.current_redshift) ** 3
        / (ds.omega_matter * (1 + ds.current_redshift) ** 3 + ds.omega_lambda)
    )

    assert np.abs(c1 - c2) / max(c1, c2) < 1e-3


@requires_file(orszag_tang)
def test_face_centered_bfields():
    # this is based on the enzo frontend test, test_face_centered_mhdct_fields
    ds = data_dir_load(orszag_tang)

    ad = ds.all_data()
    assert len(ds.index.grids) == 4

    for field, flag in NODAL_FLAGS.items():
        dims = ds.domain_dimensions
        assert_equal(ad[field].shape, (dims.prod(), 2 * sum(flag)))

        # the x and y domains are each split across 2 blocks. The z domain
        # is all located on a single block
        block_dims = (dims[0] // 2, dims[1] // 2, dims[2])
        for grid in ds.index.grids:
            assert_equal(grid[field].shape, block_dims + (2 * sum(flag),))

    # Average of face-centered fields should be the same as cell-centered field
    assert_array_equal(
        ad["enzoe", "bfieldi_x"].sum(axis=-1) / 2, ad["enzoe", "bfield_x"]
    )
    assert_array_equal(
        ad["enzoe", "bfieldi_y"].sum(axis=-1) / 2, ad["enzoe", "bfield_y"]
    )
    assert_array_equal(
        ad["enzoe", "bfieldi_z"].sum(axis=-1) / 2, ad["enzoe", "bfield_z"]
    )
