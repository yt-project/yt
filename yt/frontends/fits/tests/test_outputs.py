from yt.testing import (
    assert_equal,
    requires_file,
    requires_module,
    units_override_check,
)
from yt.utilities.answer_testing.framework import (
    data_dir_load,
    requires_ds,
    small_patch_amr,
)

from ..data_structures import (
    EventsFITSDataset,
    FITSDataset,
    SkyDataFITSDataset,
    SpectralCubeFITSDataset,
)

_fields_grs = (("fits", "temperature"),)

grs = "radio_fits/grs-50-cube.fits"


@requires_ds(grs)
def test_grs():
    ds = data_dir_load(grs, cls=SpectralCubeFITSDataset, kwargs={"nan_mask": 0.0})
    assert_equal(str(ds), "grs-50-cube.fits")
    for test in small_patch_amr(ds, _fields_grs, input_center="c", input_weight="ones"):
        test_grs.__name__ = test.description
        yield test


_fields_vels = (("fits", "velocity_x"), ("fits", "velocity_y"), ("fits", "velocity_z"))

vf = "UnigridData/velocity_field_20.fits"


@requires_module("astropy")
@requires_ds(vf)
def test_velocity_field():
    ds = data_dir_load(vf, cls=FITSDataset)
    assert_equal(str(ds), "velocity_field_20.fits")
    for test in small_patch_amr(
        ds, _fields_vels, input_center="c", input_weight="ones"
    ):
        test_velocity_field.__name__ = test.description
        yield test


acis = "xray_fits/acisf05356N003_evt2.fits.gz"

_fields_acis = (("gas", "counts_0.1-2.0"), ("gas", "counts_2.0-5.0"))


@requires_ds(acis)
def test_acis():
    from yt.frontends.fits.misc import setup_counts_fields

    ds = data_dir_load(acis, cls=EventsFITSDataset)
    ebounds = [(0.1, 2.0), (2.0, 5.0)]
    setup_counts_fields(ds, ebounds)
    assert_equal(str(ds), "acisf05356N003_evt2.fits.gz")
    for test in small_patch_amr(
        ds, _fields_acis, input_center="c", input_weight="ones"
    ):
        test_acis.__name__ = test.description
        yield test


A2052 = "xray_fits/A2052_merged_0.3-2_match-core_tmap_bgecorr.fits"

_fields_A2052 = (("fits", "flux"),)


@requires_ds(A2052)
def test_A2052():
    ds = data_dir_load(A2052, cls=SkyDataFITSDataset)
    assert_equal(str(ds), "A2052_merged_0.3-2_match-core_tmap_bgecorr.fits")
    for test in small_patch_amr(
        ds, _fields_A2052, input_center="c", input_weight="ones"
    ):
        test_A2052.__name__ = test.description
        yield test


@requires_file(vf)
def test_units_override():
    units_override_check(vf)


@requires_file(vf)
def test_FITSDataset():
    assert isinstance(data_dir_load(vf), FITSDataset)


@requires_file(grs)
def test_SpectralCubeFITSDataset():
    assert isinstance(data_dir_load(grs), SpectralCubeFITSDataset)


@requires_file(acis)
def test_EventsFITSDataset():
    assert isinstance(data_dir_load(acis), EventsFITSDataset)


@requires_file(A2052)
def test_SkyDataFITSDataset():
    assert isinstance(data_dir_load(A2052), SkyDataFITSDataset)
