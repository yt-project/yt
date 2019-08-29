"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import os

import pytest

from yt.utilities.answer_testing import utils


# Test data
ahf_halos = 'ahf_halos/snap_N64L16_135.parameter'


@pytest.fixture(scope='class')
def ds_ahf_halos():
    ds = utils.data_dir_load(ahf_halos, kwargs={'hubble_constant' : 0.7})
    assert str(ds) == os.path.basename(ahf_halos)
    return ds
