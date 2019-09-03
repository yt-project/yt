"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import os

import pytest

from yt.utilities.answer_testing import utils


# Test data
rt = "RadTube/plt00500"
star = "StarParticles/plrd01000"
LyA = "Nyx_LyA/plt00000"
RT_particles = "RT_particles/plt00050"
plasma = "PlasmaAcceleration/plt00030_v2"


@pytest.fixture(scope='class')
def ds_rt():
    ds = utils.data_dir_load(rt)
    assert str(ds) == "plt00500"
    return ds

@pytest.fixture(scope='class')
def ds_star():
    ds = utils.data_dir_load(star)
    assert str(ds) == "plrd01000"
    return ds

@pytest.fixture(scope='class')
def ds_LyA():
    ds = utils.data_dir_load(LyA)
    assert str(ds) == "plt00000"
    return ds

@pytest.fixture(scope='class')
def ds_RT_particles():
    ds = utils.data_dir_load(RT_particles)
    assert str(ds) == "plt00050"
    return ds

@pytest.fixture(scope='class')
def ds_plasma():
    ds = utils.data_dir_load(plasma)
    assert str(ds) == "plt00030_v2"
    return ds
