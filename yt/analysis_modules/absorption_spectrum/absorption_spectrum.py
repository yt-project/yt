"""
AbsorptionSpectrum class and member functions.

Author: Britton Smith <brittonsmith@gmail.com>
Affiliation: Michigan State University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008-2011 Britton Smith.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import h5py
import numpy as np

from absorption_line import tau_profile

from yt.funcs import get_pbar
from yt.utilities.physical_constants import \
    amu_cgs, boltzmann_constant_cgs, \
    speed_of_light_cgs, km_per_cm

speed_of_light_kms = speed_of_light_cgs * km_per_cm

class AbsorptionSpectrum(object):
    r"""Create an absorption spectrum object.

    Parameters
    ----------

    lambda_min : float
       lower wavelength bound in angstroms.
    lambda_max : float
       upper wavelength bound in angstroms.
    n_lambda : float
       number of wavelength bins.
    """

    def __init__(self, lambda_min, lambda_max, n_lambda):
        self.n_lambda = n_lambda
        self.tau_field = None
        self.flux_field = None
        self.spectrum_line_list = None
        self.lambda_bins = np.linspace(lambda_min, lambda_max, n_lambda)
        self.bin_width = (lambda_max - lambda_min) / float(n_lambda - 1)
        self.line_list = []
        self.continuum_list = []

    def add_line(self, label, field_name, wavelength,
                 f_value, gamma, atomic_mass,
                 label_threshold=None):
        r"""Add an absorption line to the list of lines included in the spectrum.

        Parameters
        ----------
        
        label : string
           label for the line.
        field_name : string
           field name from ray data for column densities.
        wavelength : float
           line rest wavelength in angstroms.
        f_value  : float
           line f-value.
        gamma : float
           line gamme value.
        atomic_mass : float
           mass of atom in amu.
        """
        self.line_list.append({'label': label, 'field_name': field_name,
                               'wavelength': wavelength, 'f_value': f_value,
                               'gamma': gamma, 'atomic_mass': atomic_mass,
                               'label_threshold': label_threshold})

    def add_continuum(self, label, field_name, wavelength,
                      normalization, index):
        """
        Add a continuum feature that follows a power-law.

        Parameters
        ----------

        label : string
           label for the feature.
        field_name : string
           field name from ray data for column densities.
        wavelength : float
           line rest wavelength in angstroms.
        normalization : float
           the column density normalization.
        index : float
           the power-law index for the wavelength dependence.
        """

        self.continuum_list.append({'label': label, 'field_name': field_name,
                                    'wavelength': wavelength,
                                    'normalization': normalization,
                                    'index': index})

    def make_spectrum(self, input_file, output_file='spectrum.h5',
                      line_list_file='lines.txt',
                      use_peculiar_velocity=True):
        """
        Make spectrum from ray data using the line list.

        Parameters
        ----------

        input_file : string
           path to input ray data.
        output_file : string
           path for output file.  File formats are chosen based on the filename extension.
           ``.h5`` for hdf5, ``.fits`` for fits, and everything else is ASCII.
        use_peculiar_velocity : bool
           if True, include line of sight velocity for shifting lines.
        """

        input_fields = ['dl', 'redshift', 'Temperature']
        field_data = {}
        if use_peculiar_velocity: input_fields.append('los_velocity')
        for feature in self.line_list + self.continuum_list:
            if not feature['field_name'] in input_fields:
                input_fields.append(feature['field_name'])

        input = h5py.File(input_file, 'r')
        for field in input_fields:
            field_data[field] = input[field].value
        input.close()

        self.tau_field = np.zeros(self.lambda_bins.size)
        self.spectrum_line_list = []

        self._add_lines_to_spectrum(field_data, use_peculiar_velocity)
        self._add_continua_to_spectrum(field_data, use_peculiar_velocity)

        self.flux_field = np.exp(-self.tau_field)

        if output_file.endswith('.h5'):
            self._write_spectrum_hdf5(output_file)
        elif output_file.endswith('.fits'):
            self._write_spectrum_fits(output_file)
        else:
            self._write_spectrum_ascii(output_file)
        self._write_spectrum_line_list(line_list_file)

        del field_data
        return (self.lambda_bins, self.flux_field)

    def _add_continua_to_spectrum(self, field_data, use_peculiar_velocity):
        """
        Add continuum features to the spectrum.
        """
        # Only add continuum features down to tau of 1.e-4.
        tau_min = 1.e-4

        for continuum in self.continuum_list:
            column_density = field_data[continuum['field_name']] * field_data['dl']
            delta_lambda = continuum['wavelength'] * field_data['redshift']
            if use_peculiar_velocity:
                # include factor of (1 + z) because our velocity is in proper frame.
                delta_lambda += continuum['wavelength'] * (1 + field_data['redshift']) * \
                    field_data['los_velocity'] / speed_of_light_cgs
            this_wavelength = delta_lambda + continuum['wavelength']
            right_index = np.digitize(this_wavelength, self.lambda_bins).clip(0, self.n_lambda)
            left_index = np.digitize((this_wavelength *
                                     np.power((tau_min * continuum['normalization'] /
                                               column_density), (1. / continuum['index']))),
                                    self.lambda_bins).clip(0, self.n_lambda)

            valid_continuua = np.where(((column_density /
                                         continuum['normalization']) > tau_min) &
                                       (right_index - left_index > 1))[0]
            pbar = get_pbar("Adding continuum feature - %s [%f A]: " % \
                                (continuum['label'], continuum['wavelength']),
                            valid_continuua.size)
            for i, lixel in enumerate(valid_continuua):
                line_tau = np.power((self.lambda_bins[left_index[lixel]:right_index[lixel]] /
                                     this_wavelength[lixel]), continuum['index']) * \
                                     column_density[lixel] / continuum['normalization']
                self.tau_field[left_index[lixel]:right_index[lixel]] += line_tau
                pbar.update(i)
            pbar.finish()

    def _add_lines_to_spectrum(self, field_data, use_peculiar_velocity):
        """
        Add the absorption lines to the spectrum.
        """
        # Only make voigt profile for slice of spectrum that is 10 times the line width.
        spectrum_bin_ratio = 5
        # Widen wavelength window until optical depth reaches a max value at the ends.
        max_tau = 0.001

        for line in self.line_list:
            column_density = field_data[line['field_name']] * field_data['dl']
            delta_lambda = line['wavelength'] * field_data['redshift']
            if use_peculiar_velocity:
                # include factor of (1 + z) because our velocity is in proper frame.
                delta_lambda += line['wavelength'] * (1 + field_data['redshift']) * \
                    field_data['los_velocity'] / speed_of_light_cgs
            thermal_b = km_per_cm * np.sqrt((2 * boltzmann_constant_cgs *
                                             field_data['Temperature']) /
                                            (amu_cgs * line['atomic_mass']))
            center_bins = np.digitize((delta_lambda + line['wavelength']),
                                      self.lambda_bins)

            # ratio of line width to bin width
            width_ratio = (line['wavelength'] + delta_lambda) * \
                thermal_b / speed_of_light_kms / self.bin_width

            # do voigt profiles for a subset of the full spectrum
            left_index  = (center_bins -
                           spectrum_bin_ratio * width_ratio).astype(int).clip(0, self.n_lambda)
            right_index = (center_bins +
                           spectrum_bin_ratio * width_ratio).astype(int).clip(0, self.n_lambda)

            # loop over all lines wider than the bin width
            valid_lines = np.where((width_ratio >= 1.0) &
                                   (right_index - left_index > 1))[0]
            pbar = get_pbar("Adding line - %s [%f A]: " % (line['label'], line['wavelength']),
                            valid_lines.size)
            for i, lixel in enumerate(valid_lines):
                my_bin_ratio = spectrum_bin_ratio
                while True:
                    lambda_bins, line_tau = \
                        tau_profile(line['wavelength'], line['f_value'],
                                    line['gamma'], thermal_b[lixel],
                                    column_density[lixel],
                                    delta_lambda=delta_lambda[lixel],
                                    lambda_bins=self.lambda_bins[left_index[lixel]:right_index[lixel]])
                    # Widen wavelength window until optical depth reaches a max value at the ends.
                    if (line_tau[0] < max_tau and line_tau[-1] < max_tau) or \
                      (left_index[lixel] <= 0 and right_index[lixel] >= self.n_lambda):
                        break
                    my_bin_ratio *= 2
                    left_index[lixel]  = (center_bins[lixel] -
                                          my_bin_ratio *
                                          width_ratio[lixel]).astype(int).clip(0, self.n_lambda)
                    right_index[lixel] = (center_bins[lixel] +
                                          my_bin_ratio *
                                          width_ratio[lixel]).astype(int).clip(0, self.n_lambda)
                self.tau_field[left_index[lixel]:right_index[lixel]] += line_tau
                if line['label_threshold'] is not None and \
                        column_density[lixel] >= line['label_threshold']:
                    if use_peculiar_velocity:
                        peculiar_velocity = km_per_cm * field_data['los_velocity'][lixel]
                    else:
                        peculiar_velocity = 0.0
                    self.spectrum_line_list.append({'label': line['label'],
                                                    'wavelength': (line['wavelength'] +
                                                                   delta_lambda[lixel]),
                                                    'column_density': column_density[lixel],
                                                    'b_thermal': thermal_b[lixel],
                                                    'redshift': field_data['redshift'][lixel],
                                                    'v_pec': peculiar_velocity})
                    pbar.update(i)
            pbar.finish()

            del column_density, delta_lambda, thermal_b, \
                center_bins, width_ratio, left_index, right_index

    def _write_spectrum_line_list(self, filename):
        """
        Write out list of spectral lines.
        """
        print "Writing spectral line list: %s." % filename
        self.spectrum_line_list.sort(key=lambda obj: obj['wavelength'])
        f = open(filename, 'w')
        f.write('#%-14s %-14s %-12s %-12s %-12s %-12s\n' %
                ('Wavelength', 'Line', 'N [cm^-2]', 'b [km/s]', 'z', 'v_pec [km/s]'))
        for line in self.spectrum_line_list:
            f.write('%-14.6f %-14ls %e %e %e %e.\n' % (line['wavelength'], line['label'],
                                                line['column_density'], line['b_thermal'],
                                                line['redshift'], line['v_pec']))
        f.close()

    def _write_spectrum_ascii(self, filename):
        """
        Write spectrum to an ascii file.
        """
        print "Writing spectrum to ascii file: %s." % filename
        f = open(filename, 'w')
        f.write("# wavelength[A] tau flux\n")
        for i in xrange(self.lambda_bins.size):
            f.write("%e %e %e\n" % (self.lambda_bins[i],
                                    self.tau_field[i], self.flux_field[i]))
        f.close()

    def _write_spectrum_fits(self, filename):
        """
        Write spectrum to a fits file.
        """
        try:
            import pyfits
        except:
            print "Could not import the pyfits module.  Please install pyfits."
            return

        print "Writing spectrum to fits file: %s." % filename
        col1 = pyfits.Column(name='wavelength', format='E', array=self.lambda_bins)
        col2 = pyfits.Column(name='flux', format='E', array=self.flux_field)
        cols = pyfits.ColDefs([col1, col2])
        tbhdu = pyfits.new_table(cols)
        tbhdu.writeto(filename, clobber=True)

    def _write_spectrum_hdf5(self, filename):
        """
        Write spectrum to an hdf5 file.

        """
        print "Writing spectrum to hdf5 file: %s." % filename
        output = h5py.File(filename, 'w')
        output.create_dataset('wavelength', data=self.lambda_bins)
        output.create_dataset('tau', data=self.tau_field)
        output.create_dataset('flux', data=self.flux_field)
        output.close()
