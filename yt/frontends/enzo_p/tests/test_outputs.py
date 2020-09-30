import numpy as np
import pytest

from yt.frontends.enzo_p.api import EnzoPDataset
from yt.testing import assert_array_equal, assert_equal
from yt.utilities.answer_testing.testing_utilities import create_obj
from yt.utilities.answer_testing.answer_tests import (
    field_values,
    pixelized_projection_values,
)
from yt.utilities.on_demand_imports import _h5py as h5py

# Test data
hello_world = "hello-0210/hello-0210.block_list"
ep_cosmo = "ENZOP_DD0140/ENZOP_DD0140.block_list"


# Global field info
_fields = ("density", "total_energy", "velocity_x", "velocity_y")
_pfields = (
    "particle_position_x",
    "particle_position_y",
    "particle_position_z",
    "particle_velocity_x",
    "particle_velocity_y",
    "particle_velocity_z",
)


ds_list = [
    hello_world,
    ep_cosmo,
]
fields = [
    _fields,
    _pfields,
]
objs = [
    [None, ("sphere", ("max", (0.25, "unitary")))],
    [None, ("sphere", ("max", (0.1, "unitary")))],
]
axes = [0, 1, 2]
weights = [None, "density"]


fv_pairs = []
ppv_pairs = []
sum_pairs = []


for i, ds in enumerate(ds_list):
    for field in fields[i]:
        for obj in objs[i]:
            fv_pairs.append((ds, field, obj))


for field in fields[0]:
    for obj in objs[0]:
        for weight in weights:
            for a in axes:
                ppv_pairs.append((ds_list[0], field, obj, weight, a))


for i, ds in enumerate(ds_list):
    for obj in objs[i]:
        sum_pairs.append((ds, obj))


@pytest.mark.answer_test
class TestEnzoP:
    answer_file = None
    saved_hashes = None

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f, d", fv_pairs, indirect=True)
    def test_field_values(self, ds, f, d):
        particle_type = f[0] in ds.particle_types
        if str(ds) == "ENZOP_DD0140":
            particle_type = True
        fv = field_values(ds, f, d, particle_type=particle_type)
        self.hashes.update({"field_values": fv})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f, d, w, a", ppv_pairs, indirect=True)
    def test_pixelized_projection_values(self, ds, f, d, w, a):
        ppv = pixelized_projection_values(ds, a, f, w, d)
        self.hashes.update({"pixelized_projection_values": ppv})

    @pytest.mark.parametrize("ds, d", sum_pairs, indirect=True)
    def test_sum(self, ds, d):
        dobj = create_obj(ds, d)
        s1 = dobj["ones"].sum()
        s2 = sum(mask.sum() for block, mask in dobj.blocks)
        assert_equal(s1, s2)

    @pytest.mark.parametrize("ds", [hello_world], indirect=True)
    def test_EnzoPDataset(self, ds):
        assert isinstance(ds, EnzoPDataset)

    @pytest.mark.parametrize("ds", [hello_world], indirect=True)
    def test_hierarchy(self, ds):
        fh = h5py.File(ds.index.grids[0].filename, "r")
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

    @pytest.mark.parametrize("ds", [ep_cosmo], indirect=True)
    def test_critical_density(self, ds):
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
