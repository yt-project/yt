"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest

import yt
from yt.frontends.gizmo.api import GizmoDataset


# Test data
gmhd = "gizmo_mhd_mwdisk/gizmo_mhd_mwdisk.hdf5"
g64 = "gizmo_64/output/snap_N64L16_135.hdf5"


fields = 
    {
        ("gas", "density"): None,
        ("gas", "temperature"): ('gas', 'density'),
        ("gas", "metallicity"): ('gas', 'density'),
        ("gas", "O_metallicity"): ('gas', 'density'),
        ('gas', 'velocity_magnitude'): None,
        ("deposit", "all_count"): None,
        ("deposit", "all_cic"): None,
        ("deposit", "PartType0_density"): None
    }


def pytest_generate_tests(metafunc):
    if metafunc.function.__name__ == 'test_gizmo_64':
        metafunc.parametrize('f, w', [(k, v) for k, v in fields.items()],
            ids=['dens', 'temp', 'metallicity', 'O_metallicity', 'velocity_magnitude',
            'all_count', 'all_cic', 'PartType0Dens'])

@pytest.fixture(scope='class')
def ds_gmhd():
    gmhd_bbox = [[-400, 400]] * 3
    ds = yt.load(gmhd, bounding_box=gmhd_bbox, unit_system='code')
    return ds

@pytest.fixture(scope='class')
def ds_g64():
    ds = yt.load(g64)
    assert isinstance(ds, GizmoDataset)
    return ds
