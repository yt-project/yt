import numpy as np  # NOQA
import pytest

import yt  # NOQA
from yt.frontends.amrvac.api import AMRVACDataset, AMRVACGrid
from yt.testing import requires_file
from yt.units import YTArray
from yt.utilities.answer_testing.testing_utilities import data_dir_load
from yt.utilities.answer_testing.testing_utilities import requires_ds
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
axes = [0, 1, 2]
weights = [None, "density"]
objs = [None, ("sphere", ("max", (0.1, "unitary")))]


def _get_fields_to_check(fname):
    # This function is called during test collection. If this frontend
    # is not being run, and therefore the data isn't present, this try
    # except block prevents pytest from failing needlessly. Since this
    # isn't a fixture, we can't use pytest.skip. As such, in order to
    # indicate to the test function that the data wasn't found, we
    # return None for the ds and None for the field
    try:
        ds = data_dir_load(fname)
    except FileNotFoundError:
        return [None, [None,]]
    fields = [("gas", "density"), ("gas", "velocity_magnitude")]
    raw_fields_labels = [fname for ftype, fname in ds.field_list]
    if "b1" in raw_fields_labels:
        fields.append(("gas", "magnetic_energy_density"))
    if "e" in raw_fields_labels:
        fields.append(("gas", "energy_density"))
    if "rhod1" in raw_fields_labels:
        fields.append(("gas", "total_dust_density"))
        # note : not hitting dust velocity fields
    return [ds, fields]


def get_pairs():
    """
    In order to properly parametrize the test, we need to group each ds
    with each field that is to be checked for that ds. So we want:

    [(ds1, ds1_field1), (ds1, ds1_field2), ..., (ds2, ds2_field1), ...]

    This function creates and returns the above list.
    """
    # Start by grouping each ds with all of the fields to be checked
    # for that ds
    ds_fieldlist_pairs = []
    for ds_name in ds_list:
        pair = _get_fields_to_check(ds_name)
        ds_fieldlist_pairs.append(pair)
    # Now we can "distribute" the dataset name to each field in that dataset's
    # list of fields in order to get the form required by pytest:
    # [(ds1, ds1_f1), (ds1, ds1_f2), ..., (ds2, ds2_f1), ... (dsN, dsN_fN)]
    pairs = []
    for ds, field_list in ds_fieldlist_pairs:
        for field in field_list:
            pairs.append((ds, field))
    return pairs


@pytest.mark.answer_test
class TestAMRVAC:
    answer_file = None
    saved_hashes = None

    @requires_file(khi_cartesian_2D)
    def test_AMRVACDataset(self):
        assert isinstance(data_dir_load(khi_cartesian_2D), AMRVACDataset)

    @requires_ds(blastwave_cartesian_3D)
    def test_domain_size(self):
        # Check for correct box size, see bw_3d.par
        ds = data_dir_load(blastwave_cartesian_3D)
        for lb in ds.domain_left_edge:
            assert int(lb) == 0
        for rb in ds.domain_right_edge:
            assert int(rb) == 2
        for w in ds.domain_width:
            assert int(w) == 2

    @requires_file(blastwave_cartesian_3D)
    def test_grid_attributes(self):
        # Check various grid attributes
        ds = data_dir_load(blastwave_cartesian_3D)
        grids = ds.index.grids
        assert ds.index.max_level == 2
        for g in grids:
            assert isinstance(g, AMRVACGrid)
            assert isinstance(g.LeftEdge, YTArray)
            assert isinstance(g.RightEdge, YTArray)
            assert isinstance(g.ActiveDimensions, np.ndarray)
            assert isinstance(g.Level, (int, np.integer))

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", ds_list, indirect=True)
    def test_grid_hierarchy_parentage_relationships(self, ds):
        self.hashes.update({"grid_hierarchy": grid_hierarchy(ds)})
        self.hashes.update({"parentage_relationships": parentage_relationships(ds)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f", get_pairs(), indirect=True)
    def test_grid_values(self, f, ds):
        if ds and f is None:
            pytest.skip(f"Data `{str(ds)}` not found; skipping.")
        self.hashes.update({"grid_values": grid_values(ds, f)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f", get_pairs(), indirect=True)
    @pytest.mark.parametrize("d", objs, indirect=True)
    def test_field_values(self, d, f, ds):
        if ds and f is None:
            pytest.skip(f"Data `{str(ds)}` not found; skipping.")
        self.hashes.update({"field_values": field_values(ds, f, d)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f", get_pairs(), indirect=True)
    @pytest.mark.parametrize("d", objs, indirect=True)
    @pytest.mark.parametrize("a", axes, indirect=True)
    @pytest.mark.parametrize("w", weights, indirect=True)
    def test_projection_values(self, a, d, w, f, ds):
        if ds and f is None:
            pytest.skip(f"Data `{str(ds)}` not found; skipping.")
        self.hashes.update({"projection_values": projection_values(ds, a, f, w, d)})
