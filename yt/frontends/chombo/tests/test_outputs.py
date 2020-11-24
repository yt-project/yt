import pytest

from yt.frontends.chombo.api import ChomboDataset, Orion2Dataset, PlutoDataset
from yt.testing import requires_file, units_override_check
from yt.utilities.answer_testing.answer_tests import (
    field_values,
    grid_hierarchy,
    grid_values,
    parentage_relationships,
    projection_values,
)

# Test data
gc = "GaussianCloud/data.0077.3d.hdf5"
tb = "TurbBoxLowRes/data.0005.3d.hdf5"
iso = "IsothermalSphere/data.0000.3d.hdf5"
zp = "ZeldovichPancake/plt32.2d.hdf5"
kho = "KelvinHelmholtz/data.0004.hdf5"


ds_list = [
    gc,
    tb,
    iso,
    zp,
    kho,
]
axes = [0, 1, 2]
weights = [None, "density"]
weights_zp = [None, "rhs"]
objs = [None, ("sphere", ("max", (0.1, "unitary")))]
obj_zp = [None, ("sphere", ("c", (0.1, "unitary")))]
fields = ["density", "velocity_magnitude", "magnetic_field_x"]
fields_zp = ["rhs", "phi"]

pair_list = [
    [gc, fields, objs, weights],
    [tb, fields, objs, weights],
    [iso, fields, objs, weights],
    [zp, fields_zp, obj_zp, weights_zp],
    [kho, fields, objs, weights],
]

gv_pairs = [(i[0], f) for i in ds_list for f in i[1]]
fv_pairs = [(i[0], f, d) for i in pair_list for f in i[1] for d in i[2]]
pv_pairs = [
    (i[0], f, d, w) for i in pair_list for f in i[1] for d in i[2] for w in i[3]
]


@pytest.mark.answer_test
class TestChombo:
    answer_file = None
    saved_hashes = None

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", ds_list, indirect=True)
    def test_grid_hierarchy_parentage_relationships(self, ds):
        self.hashes.update({"grid_hierarchy": grid_hierarchy(ds)})
        self.hashes.update({"parentage_relationships": parentage_relationships(ds)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f", gv_pairs, indirect=True)
    def test_grid_values(self, f, ds):
        self.hashes.update({"grid_values": grid_values(ds, f)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f, d", fv_pairs, indirect=True)
    def test_field_values(self, d, f, ds):
        self.hashes.update({"field_values": field_values(ds, f, d)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f, d, w", pv_pairs, indirect=True)
    @pytest.mark.parametrize("a", axes, indirect=True)
    def test_projection_values(self, a, d, w, f, ds):
        self.hashes.update({"projection_values": projection_values(ds, a, f, w, d)})

    @pytest.mark.parametrize("ds", [zp], indirect=True)
    def test_ChomboDataset(self, ds):
        assert isinstance(ds, ChomboDataset)

    @pytest.mark.parametrize("ds", [gc], indirect=True)
    def test_Orion2Dataset(self, ds):
        assert isinstance(ds, Orion2Dataset)

    @pytest.mark.parametrize("ds", [kho], indirect=True)
    def test_PlutoDataset(self, ds):
        assert isinstance(ds, PlutoDataset)

    @requires_file(zp)
    def test_units_override_zp(self):
        units_override_check(zp)

    @requires_file(gc)
    def test_units_override_gc(self):
        units_override_check(gc)

    @requires_file(kho)
    def test_units_override_kho(self):
        units_override_check(kho)
