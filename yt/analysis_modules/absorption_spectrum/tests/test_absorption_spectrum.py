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

import numpy as np
from yt.testing import \
    assert_allclose_units, requires_file, requires_module, \
    assert_almost_equal, assert_array_almost_equal
from yt.analysis_modules.absorption_spectrum.absorption_line import \
    voigt_old, voigt_scipy
from yt.analysis_modules.absorption_spectrum.api import AbsorptionSpectrum
from yt.analysis_modules.cosmological_observation.api import LightRay
from yt.config import ytcfg
import tempfile
import os
import shutil
from yt.utilities.on_demand_imports import \
    _h5py as h5

test_dir = ytcfg.get("yt", "test_data_dir")

COSMO_PLUS = "enzo_cosmology_plus/AMRCosmology.enzo"
COSMO_PLUS_SINGLE = "enzo_cosmology_plus/RD0009/RD0009"
HI_SPECTRUM_COSMO = "absorption_spectrum_data/enzo_lyman_alpha_cosmo_spec.h5"
HI_SPECTRUM_COSMO_FILE = os.path.join(test_dir, HI_SPECTRUM_COSMO)
HI_SPECTRUM = "absorption_spectrum_data/enzo_lyman_alpha_spec.h5"
HI_SPECTRUM_FILE = os.path.join(test_dir, HI_SPECTRUM)

@requires_file(COSMO_PLUS)
@requires_file(HI_SPECTRUM_COSMO)
def test_absorption_spectrum_cosmo():
    """
    This test generates an absorption spectrum from a cosmological light ray
    """

    # Set up in a temp dir
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    lr = LightRay(COSMO_PLUS, 'Enzo', 0.0, 0.03)

    lr.make_light_ray(seed=1234567,
                      fields=['temperature', 'density', 'H_number_density'],
                      data_filename='lightray.h5')

    sp = AbsorptionSpectrum(900.0, 1800.0, 10000)

    my_label = 'HI Lya'
    field = 'H_number_density'
    wavelength = 1215.6700  # Angstromss
    f_value = 4.164E-01
    gamma = 6.265e+08
    mass = 1.00794

    sp.add_line(my_label, field, wavelength, f_value,
                gamma, mass, label_threshold=1.e10)

    my_label = 'HI Lya'
    field = 'H_number_density'
    wavelength = 912.323660  # Angstroms
    normalization = 1.6e17
    index = 3.0

    sp.add_continuum(my_label, field, wavelength, normalization, index)

    wavelength, flux = sp.make_spectrum('lightray.h5',
                                        output_file='spectrum.h5',
                                        line_list_file='lines.txt',
                                        use_peculiar_velocity=True)

    # load just-generated hdf5 file of spectral data (for consistency)
    f_new = h5.File('spectrum.h5', 'r')

    # load standard data for comparison
    f_old = h5.File(HI_SPECTRUM_COSMO_FILE, 'r')

    # compare between standard data and current data for each array saved 
    # (wavelength, flux, tau)
    for key in f_old.keys():
        assert_array_almost_equal(f_new[key].value, f_old[key].value, 10)

    # clean up
    os.chdir(curdir)
    shutil.rmtree(tmpdir)

@requires_file(COSMO_PLUS_SINGLE)
@requires_file(HI_SPECTRUM)
def test_absorption_spectrum_non_cosmo():
    """
    This test generates an absorption spectrum from a non-cosmological light ray
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

    sp = AbsorptionSpectrum(1200.0, 1300.0, 10001)

    my_label = 'HI Lya'
    field = 'H_number_density'
    wavelength = 1215.6700  # Angstromss
    f_value = 4.164E-01
    gamma = 6.265e+08
    mass = 1.00794

    sp.add_line(my_label, field, wavelength, f_value,
                gamma, mass, label_threshold=1.e10)

    wavelength, flux = sp.make_spectrum('lightray.h5',
                                        output_file='spectrum.h5',
                                        line_list_file='lines.txt',
                                        use_peculiar_velocity=True)

    # load just-generated hdf5 file of spectral data (for consistency)
    f_new = h5.File('spectrum.h5', 'r')

    # load standard data for comparison
    f_old = h5.File(HI_SPECTRUM_FILE, 'r')

    # compare between standard data and current data for each array saved 
    # (wavelength, flux, tau)
    for key in f_old.keys():
        assert_array_almost_equal(f_new[key].value, f_old[key].value, 10)

    # clean up
    os.chdir(curdir)
    shutil.rmtree(tmpdir)

@requires_file(COSMO_PLUS_SINGLE)
def test_equivalent_width_conserved():
    """
    This tests that the equivalent width of the optical depth is conserved 
    regardless of the bin width employed in wavelength space.
    Unresolved lines should still deposit optical depth into the spectrum.
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

    my_label = 'HI Lya'
    field = 'H_number_density'
    wave = 1215.6700  # Angstromss
    f_value = 4.164E-01
    gamma = 6.265e+08
    mass = 1.00794

    lambda_min= 1200
    lambda_max= 1300
    lambda_bin_widths = [1e-3, 1e-2, 1e-1, 1e0, 1e1]
    total_tau = []

    for lambda_bin_width in lambda_bin_widths:
        n_lambda = ((lambda_max - lambda_min)/ lambda_bin_width) + 1
        sp = AbsorptionSpectrum(lambda_min=lambda_min, lambda_max=lambda_max, 
                                n_lambda=n_lambda)
        sp.add_line(my_label, field, wave, f_value, gamma, mass)
        wavelength, flux = sp.make_spectrum('lightray.h5')
        total_tau.append((lambda_bin_width * sp.tau_field).sum())
        
    # assure that the total tau values are all within 1e-5 of each other
    for tau in total_tau:
        assert_almost_equal(tau, total_tau[0], 5)

    # clean up
    os.chdir(curdir)
    shutil.rmtree(tmpdir)


@requires_file(COSMO_PLUS_SINGLE)
@requires_module("astropy")
def test_absorption_spectrum_fits():
    """
    This test generates an absorption spectrum and saves it as a fits file.
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

    sp = AbsorptionSpectrum(900.0, 1800.0, 10000)

    my_label = 'HI Lya'
    field = 'H_number_density'
    wavelength = 1215.6700  # Angstromss
    f_value = 4.164E-01
    gamma = 6.265e+08
    mass = 1.00794

    sp.add_line(my_label, field, wavelength, f_value,
                gamma, mass, label_threshold=1.e10)

    my_label = 'HI Lya'
    field = 'H_number_density'
    wavelength = 912.323660  # Angstroms
    normalization = 1.6e17
    index = 3.0

    sp.add_continuum(my_label, field, wavelength, normalization, index)

    wavelength, flux = sp.make_spectrum('lightray.h5',
                                        output_file='spectrum.fits',
                                        line_list_file='lines.txt',
                                        use_peculiar_velocity=True)

    # clean up
    os.chdir(curdir)
    shutil.rmtree(tmpdir)


@requires_module("scipy")
def test_voigt_profiles():
    a = 1.7e-4
    x = np.linspace(5.0, -3.6, 60)
    yield assert_allclose_units, voigt_old(a, x), voigt_scipy(a, x), 1e-8
