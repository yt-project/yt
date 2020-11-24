import pytest

from yt.frontends.fits.data_structures import (
    EventsFITSDataset,
    FITSDataset,
    SkyDataFITSDataset,
    SpectralCubeFITSDataset,
)
from yt.frontends.fits.misc import setup_counts_fields
from yt.testing import requires_file, units_override_check
from yt.utilities.answer_testing.answer_tests import (
    field_values,
    grid_hierarchy,
    grid_values,
    parentage_relationships,
    projection_values,
)
from yt.utilities.answer_testing.testing_utilities import data_dir_load

# Test data
grs = "radio_fits/grs-50-cube.fits"
vf = "UnigridData/velocity_field_20.fits"
acis = "xray_fits/acisf05356N003_evt2.fits"
A2052 = "xray_fits/A2052_merged_0.3-2_match-core_tmap_bgecorr.fits"


def get_acis():
    try:
        ds = data_dir_load(acis, cls=EventsFITSDataset)
        ebounds = [(0.1, 2.0), (2.0, 5.0)]
        setup_counts_fields(ds, ebounds)
        return ds
    except FileNotFoundError:
        return "/does/not/exist"


grs_kwargs = {"kwargs": {"nan_mask": 0.0}, "cls": SpectralCubeFITSDataset}
vf_kwargs = {"cls": FITSDataset}
A2052_kwargs = {"cls": SkyDataFITSDataset}

_fields_grs = ("temperature",)
_fields_vels = ("velocity_x", "velocity_y", "velocity_z")
_fields_acis = ("counts_0.1-2.0", "counts_2.0-5.0")
_fields_A2052 = ("flux",)

axes = [0, 1, 2]
objs = [None, ("sphere", ("c", (0.1, "unitary")))]
weights = [None, "ones"]
fields = [
    _fields_grs,
    _fields_vels,
    _fields_acis,
    _fields_A2052,
]
ds_list = [
    [grs, grs_kwargs],
    [vf, vf_kwargs],
    get_acis(),
    [A2052, A2052_kwargs],
]


def get_pairs():
    pairs = []
    for i, ds in enumerate(ds_list):
        for f in fields[i]:
            pairs.append((ds, f))
    return pairs


@pytest.mark.answer_test
class TestFits:
    answer_file = None
    saved_hashes = None

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", ds_list, indirect=True)
    def test_grid_hierarchy_parentage_relationships(self, ds):
        self.hashes.update({"grid_hierarchy": grid_hierarchy(ds)})
        self.hashes.update({"parentage_relationships": parentage_relationships(ds)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f", get_pairs(), indirect=True)
    def test_grid_values(self, f, ds):
        self.hashes.update({"grid_values": grid_values(ds, f)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f", get_pairs(), indirect=True)
    @pytest.mark.parametrize("d", objs, indirect=True)
    def test_field_values(self, d, f, ds):
        self.hashes.update({"field_values": field_values(ds, f, d)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f", get_pairs(), indirect=True)
    @pytest.mark.parametrize("d", objs, indirect=True)
    @pytest.mark.parametrize("a", axes, indirect=True)
    @pytest.mark.parametrize("w", weights, indirect=True)
    def test_projection_values(self, a, d, w, f, ds):
        self.hashes.update({"projection_values": projection_values(ds, a, f, w, d)})

    @requires_file(vf)
    def test_units_override(self):
        units_override_check(vf)

    @pytest.mark.parametrize("ds", [[vf, vf_kwargs]], indirect=True)
    def test_FITSDataset(self, ds):
        assert isinstance(ds, FITSDataset)

    @pytest.mark.parametrize("ds", [[grs, grs_kwargs]], indirect=True)
    def test_SpectralCubeFITSDataset(self, ds):
        assert isinstance(ds, SpectralCubeFITSDataset)

    @pytest.mark.parametrize("ds", [get_acis(),], indirect=True)
    def test_EventsFITSDataset(self, ds):
        assert isinstance(ds, EventsFITSDataset)

    @pytest.mark.parametrize("ds", [[A2052, A2052_kwargs]], indirect=True)
    def test_SkyDataFITSDataset(self, ds):
        assert isinstance(ds, SkyDataFITSDataset)
