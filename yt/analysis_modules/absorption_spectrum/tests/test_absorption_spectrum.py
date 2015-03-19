"""
Unit test for the AbsorptionSpectrum analysis module
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import yt
from yt.testing import *
from yt.analysis_modules.absorption_spectrum.api import AbsorptionSpectrum
from yt.analysis_modules.cosmological_observation.api import LightRay
import tempfile
import os
import shutil

COSMO_PLUS = "enzo_cosmology_plus/AMRCosmology.enzo"

@requires_file(COSMO_PLUS)
def test_absorption_spectrum():
    """
    This test is simply following the description in the docs for how to
    generate an absorption spectrum from a cosmological light ray for one
    of the sample datasets
    """

    # Set up in a temp dir
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    lr = LightRay(COSMO_PLUS, 'Enzo', 0.0, 0.1)

    lr.make_light_ray(seed=1234567,
                      fields=['temperature', 'density', 'H_number_density'],
                      get_los_velocity=True,
                      data_filename='lightray.h5')

    sp = AbsorptionSpectrum(900.0, 1800.0, 10000)

    my_label = 'HI Lya'
    field = 'H_number_density'
    wavelength = 1215.6700 # Angstromss
    f_value = 4.164E-01
    gamma = 6.265e+08
    mass = 1.00794

    sp.add_line(my_label, field, wavelength, f_value, gamma, mass, label_threshold=1.e10)

    my_label = 'HI Lya'
    field = 'H_number_density'
    wavelength = 912.323660 # Angstroms
    normalization = 1.6e17
    index = 3.0

    sp.add_continuum(my_label, field, wavelength, normalization, index)

    wavelength, flux = sp.make_spectrum('lightray.h5', output_file='spectrum.txt',
                                        line_list_file='lines.txt',
                                        use_peculiar_velocity=True)

    # clean up
    os.chdir(curdir)
    shutil.rmtree(tmpdir)


@requires_file(COSMO_PLUS)
@requires_module("astropy")
def test_absorption_spectrum_fits():
    """
    This test is simply following the description in the docs for how to
    generate an absorption spectrum from a cosmological light ray for one
    of the sample datasets.  Outputs to fits file if astropy is installed.
    """

    # Set up in a temp dir
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    lr = LightRay(COSMO_PLUS, 'Enzo', 0.0, 0.1)

    lr.make_light_ray(seed=1234567,
                      fields=['temperature', 'density', 'H_number_density'],
                      get_los_velocity=True,
                      data_filename='lightray.h5')

    sp = AbsorptionSpectrum(900.0, 1800.0, 10000)

    my_label = 'HI Lya'
    field = 'H_number_density'
    wavelength = 1215.6700 # Angstromss
    f_value = 4.164E-01
    gamma = 6.265e+08
    mass = 1.00794

    sp.add_line(my_label, field, wavelength, f_value, gamma, mass, label_threshold=1.e10)

    my_label = 'HI Lya'
    field = 'H_number_density'
    wavelength = 912.323660 # Angstroms
    normalization = 1.6e17
    index = 3.0

    sp.add_continuum(my_label, field, wavelength, normalization, index)

    wavelength, flux = sp.make_spectrum('lightray.h5', output_file='spectrum.fits',
                                        line_list_file='lines.txt',
                                        use_peculiar_velocity=True)

    # clean up
    os.chdir(curdir)
    shutil.rmtree(tmpdir)
