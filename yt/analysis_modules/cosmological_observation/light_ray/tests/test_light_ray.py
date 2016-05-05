"""
Unit test for the light_ray analysis module
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from yt.testing import \
    assert_allclose_units, requires_file, requires_module, \
    assert_almost_equal
from yt.analysis_modules.absorption_spectrum.absorption_line import \
    voigt_old, voigt_scipy
from yt.analysis_modules.absorption_spectrum.api import AbsorptionSpectrum
from yt.analysis_modules.cosmological_observation.api import LightRay
from yt.utilities.answer_testing.framework import \
    GenericArrayTest, \
    requires_answer_testing
import tempfile
import os
import shutil
from yt.utilities.on_demand_imports import \
    _h5py as h5
from yt import load


COSMO_PLUS = "enzo_cosmology_plus/AMRCosmology.enzo"
COSMO_PLUS_SINGLE = "enzo_cosmology_plus/RD0009/RD0009"


@requires_file(COSMO_PLUS)
def test_light_ray_cosmo():
    """
    This test generates a cosmological light ray
    """
    # Set up in a temp dir
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    lr = LightRay(COSMO_PLUS, 'Enzo', 0.0, 0.03)

    lr.make_light_ray(seed=1234567,
                      fields=['temperature', 'density', 'H_number_density'],
                      data_filename='lightray.h5')

    # clean up
    os.chdir(curdir)
    shutil.rmtree(tmpdir)


@requires_file(COSMO_PLUS_SINGLE)
def test_light_ray_non_cosmo():
    """
    This test generates a non-cosmological light ray
    """

    # Set up in a temp dir
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    lr = LightRay(COSMO_PLUS_SINGLE)

    ray_start = [0,0,0]
    ray_end = [1,1,1]
    lr.make_light_ray(start_position=ray_start, end_position=ray_end,
                      fields=['temperature', 'density', 'H_number_density'],
                      data_filename='lightray.h5')

    # clean up
    os.chdir(curdir)
    shutil.rmtree(tmpdir)

@requires_file(COSMO_PLUS_SINGLE)
def test_light_ray_non_cosmo_from_dataset():
    """
    This test generates a non-cosmological light ray created from an already
    loaded dataset
    """

    # Set up in a temp dir
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    ds = load(COSMO_PLUS_SINGLE)
    lr = LightRay(ds)

    ray_start = [0,0,0]
    ray_end = [1,1,1]
    lr.make_light_ray(start_position=ray_start, end_position=ray_end,
                      fields=['temperature', 'density', 'H_number_density'],
                      data_filename='lightray.h5')

    # clean up
    os.chdir(curdir)
    shutil.rmtree(tmpdir)

