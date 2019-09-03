"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import os

import pytest

import yt
from yt.utilities.answer_testing import utils


# Test data
gmhd = "gizmo_mhd_mwdisk/gizmo_mhd_mwdisk.hdf5"


@pytest.fixture(scope='class')
def ds_gmhd():
    gmhd_bbox = [[-400, 400]] * 3
    ds = yt.load(gmhd, bounding_box=gmhd_bbox, unit_system='code')
    return ds
