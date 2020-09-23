import numpy as np  # NOQA
import pytest

import yt  # NOQA
from yt.frontends.amrvac.api import AMRVACDataset, AMRVACGrid
from yt.testing import requires_file
from yt.units import YTArray
from yt.utilities.answer_testing import utils
from yt.utilities.answer_testing.answer_tests import (
    field_values,
    grid_hierarchy,
    grid_values,
    parentage_relationships,
    projection_values,
)

# Test parameters
blastwave_spherical_2D = "amrvac/bw_2d0000.dat"
khi_cartesian_2D = "amrvac/kh_2d0000.dat"
khi_cartesian_3D = "amrvac/kh_3D0000.dat"
jet_cylindrical_25D = "amrvac/Jet0003.dat"
riemann_cartesian_175D = "amrvac/R_1d0005.dat"
blastwave_cartesian_3D = "amrvac/bw_3d0000.dat"
blastwave_polar_2D = "amrvac/bw_polar_2D0000.dat"
blastwave_cylindrical_3D = "amrvac/bw_cylindrical_3D0000.dat"
rmi_cartesian_dust_2D = "amrvac/Richtmyer_Meshkov_dust_2D/RM2D_dust_Kwok0000.dat"


ds_list = [
    blastwave_spherical_2D,
    khi_cartesian_2D,
    khi_cartesian_3D,
    jet_cylindrical_25D,
    riemann_cartesian_175D,
    blastwave_cartesian_3D,
    blastwave_polar_2D,
    blastwave_cylindrical_3D,
    rmi_cartesian_dust_2D,
]
a_list = [0, 1, 2]
w_list = [None, "density"]
dso_list = [None, ("sphere", ("max", (0.1, "unitary")))]


def _get_fields_to_check(fname):
    # This function is called during test collection. If this frontend
    # is not being run, and therefore the data isn't present, this try
    # except block prevents pytest from failing needlessly
    try:
        ds = utils.data_dir_load(fname)
    except FileNotFoundError:
        return [None, [None,]]
    fields = [("gas", "density"), ("gas", "velocity_magnitude")]
    field_ids = ["density", "velocity_magnitude"]
    raw_fields_labels = [fname for ftype, fname in ds.field_list]
    if "b1" in raw_fields_labels:
        fields.append(("gas", "magnetic_energy_density"))
        field_ids.append("magnetic_energy_density")
    if "e" in raw_fields_labels:
        fields.append(("gas", "energy_density"))
        field_ids.append("energy_density")
    if "rhod1" in raw_fields_labels:
        fields.append(("gas", "total_dust_density"))
        field_ids.append("total_dust_density")
        # note : not hitting dust velocity fields
    # return [fields, field_ids]
    return [ds, fields]


def get_pairs():
    f_list = [_get_fields_to_check(fname) for fname in ds_list]
    pairs = [(i[0], f) for i in f_list for f in i[1]]
    return pairs


@pytest.mark.answer_test
class TestAMRVAC:
    answer_file = None
    saved_hashes = None

    @requires_file(khi_cartesian_2D)
    def test_AMRVACDataset(self):
        assert isinstance(utils.data_dir_load(khi_cartesian_2D), AMRVACDataset)

    @utils.requires_ds(blastwave_cartesian_3D)
    def test_domain_size(self):
        # "Check for correct box size, see bw_3d.par"
        ds = utils.data_dir_load(blastwave_cartesian_3D)
        for lb in ds.domain_left_edge:
            assert int(lb) == 0
        for rb in ds.domain_right_edge:
            assert int(rb) == 2
        for w in ds.domain_width:
            assert int(w) == 2

    @requires_file(blastwave_cartesian_3D)
    def test_grid_attributes(self):
        # "Check various grid attributes"
        ds = utils.data_dir_load(blastwave_cartesian_3D)
        grids = ds.index.grids
        assert ds.index.max_level == 2
        for g in grids:
            assert isinstance(g, AMRVACGrid)
            assert isinstance(g.LeftEdge, YTArray)
            assert isinstance(g.RightEdge, YTArray)
            assert isinstance(g.ActiveDimensions, np.ndarray)
            assert isinstance(g.Level, (np.int32, np.int64, int))

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", ds_list, indirect=True)
    def test_gh_pr(self, ds):
        self.hashes.update({"grid_hierarchy": grid_hierarchy(ds)})
        self.hashes.update({"parentage_relationships": parentage_relationships(ds)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f", get_pairs(), indirect=True)
    def test_gv(self, f, ds):
        if ds and f is None:
            pytest.skip("Data not found.")
        self.hashes.update({"grid_values": grid_values(ds, f)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f", get_pairs(), indirect=True)
    @pytest.mark.parametrize("d", dso_list, indirect=True)
    def test_fv(self, d, f, ds):
        if ds and f is None:
            pytest.skip("Data not found.")
        self.hashes.update({"field_values": field_values(ds, f, d)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f", get_pairs(), indirect=True)
    @pytest.mark.parametrize("d", dso_list, indirect=True)
    @pytest.mark.parametrize("a", a_list, indirect=True)
    @pytest.mark.parametrize("w", w_list, indirect=True)
    def test_pv(self, a, d, w, f, ds):
        if ds and f is None:
            pytest.skip("Data not found.")
        self.hashes.update({"projection_values": projection_values(ds, a, f, w, d)})
