"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest

from yt.frontends.fits.data_structures import FITSDataset, \
    SpectralCubeFITSDataset, \
    SkyDataFITSDataset, \
    EventsFITSDataset
from yt.utilities.answer_testing import utils


# Test data
grs = "radio_fits/grs-50-cube.fits"
vf = "UnigridData/velocity_field_20.fits"
acis = "xray_fits/acisf05356N003_evt2.fits"
A2052 = "xray_fits/A2052_merged_0.3-2_match-core_tmap_bgecorr.fits"


axes = [0, 1, 2]
center = "c"
ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
weights = [None, "ones"]


def pytest_generate_tests(metafunc):
    param_tests = ['test_grs', 'test_velocity_field', 'test_acis', 'test_2052']
    if metafunc.function.__name__ == 'test_grs':
        fields = ("temperature",)
        ids = ['temperature']
    if metafunc.function.__name__ == 'test_velocity_field':
        fields = ("velocity_x","velocity_y","velocity_z")
        ids = ['vx', 'vy', 'vz']
    if metafunc.function.__name__ == 'test_acis':
        fields = ("counts_0.1-2.0", "counts_2.0-5.0")
        ids='counts.1-2', 'counts2-5']
    if metafunc.function.__name__ == 'test_2052':
        fields = ("flux",)
        ids = ['flux']
    if metafunc.function.__name__ in param_tests:
        metafunc.parametrize('f', fields, ids=ids)
        metafunc.parametrize('a', axes, ids=['0', '1', '2'])
        metafunc.parametrize('d', ds_objs, ids=['None', 'sphere'])
        metafunc.parametrize('w', weights, ids=['None', 'ones'])


@pytest.fixture(scope='class')
def ds_grs():
    ds = utils.data_dir_load(grs, cls=SpectralCubeFITSDataset, kwargs={"nan_mask":0.0})
    assert str(ds) == "grs-50-cube.fits"
    return ds

@pytest.fixture(scope='class')
def ds_vf():
    ds = utils.data_dir_load(vf, cls=FITSDataset)
    assert str(ds) == "velocity_field_20.fits"
    return ds

@pytest.fixture(scope='class')
def ds_acis():
    ds = utils.data_dir_load(acis, cls=EventsFITSDataset)
    from yt.frontends.fits.misc import setup_counts_fields
    ebounds = [(0.1, 2.0), (2.0, 5.0)]
    setup_counts_fields(ds, ebounds)
    assert_equal(str(ds), "acisf05356N003_evt2.fits")
    return ds

@pytest.fixture(scope='class')
def ds_A2052():
    ds = utils.data_dir_load(A2052, cls=SkyDataFITSDataset)
    assert str(ds) == "A2052_merged_0.3-2_match-core_tmap_bgecorr.fits"
    return ds
