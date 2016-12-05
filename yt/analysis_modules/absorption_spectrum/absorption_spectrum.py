"""
AbsorptionSpectrum class and member functions.



"""

from __future__ import absolute_import

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.on_demand_imports import _h5py as h5py
import numpy as np

from .absorption_line import tau_profile

from yt.extern.six import string_types
from yt.convenience import load
from yt.funcs import get_pbar, mylog
from yt.units.yt_array import YTArray, YTQuantity
from yt.utilities.physical_constants import \
    boltzmann_constant_cgs, \
    speed_of_light_cgs
from yt.utilities.on_demand_imports import _astropy
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    _get_comm, \
    parallel_objects, \
    parallel_root_only

pyfits = _astropy.pyfits

class AbsorptionSpectrum(object):
    r"""Create an absorption spectrum object.

    Parameters
    ----------

    lambda_min : float
       lower wavelength bound in angstroms.
    lambda_max : float
       upper wavelength bound in angstroms.
    n_lambda : int
       number of wavelength bins.
    """

    def __init__(self, lambda_min, lambda_max, n_lambda):
        self.n_lambda = int(n_lambda)
        # lambda, flux, and tau are wavelength, flux, and optical depth
        self.lambda_min = lambda_min
        self.lambda_max = lambda_max
        self.lambda_field = YTArray(np.linspace(lambda_min, lambda_max,
                                    n_lambda), "angstrom")
        self.tau_field = None
        self.flux_field = None
        self.absorbers_list = None
        self.bin_width = YTQuantity((lambda_max - lambda_min) /
                                    float(n_lambda - 1), "angstrom")
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
                               'wavelength': YTQuantity(wavelength, "angstrom"),
                               'f_value': f_value,
                               'gamma': gamma,
                               'atomic_mass': YTQuantity(atomic_mass, "amu"),
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

    def make_spectrum(self, input_file, output_file=None,
                      line_list_file=None, output_absorbers_file=None,
                      use_peculiar_velocity=True,
                      subgrid_resolution=10, observing_redshift=0.,
                      njobs="auto"):
        """
        Make spectrum from ray data using the line list.

        Parameters
        ----------

        input_file : string or dataset
           path to input ray data or a loaded ray dataset
        output_file : optional, string
           Option to save a file containing the wavelength, flux, and optical
           depth fields.  File formats are chosen based on the filename
           extension. ``.h5`` for hdf5, ``.fits`` for fits, and everything
           else is ASCII.
           Default: None
        output_absorbers_file : optional, string
           Option to save a text file containing all of the absorbers and
           corresponding wavelength and redshift information.
           For parallel jobs, combining the lines lists can be slow so it
           is recommended to set to None in such circumstances.
           Default: None
        use_peculiar_velocity : optional, bool
           if True, include peculiar velocity for calculating doppler redshift
           to shift lines.  Requires similar flag to be set in LightRay
           generation.
           Default: True
        subgrid_resolution : optional, int
           When a line is being added that is unresolved (ie its thermal
           width is less than the spectral bin width), the voigt profile of
           the line is deposited into an array of virtual wavelength bins at
           higher resolution.  The optical depth from these virtual bins is
           integrated and then added to the coarser spectral wavelength bin.
           The subgrid_resolution value determines the ratio between the
           thermal width and the bin width of the virtual bins.  Increasing
           this value yields smaller virtual bins, which increases accuracy,
           but is more expensive.  A value of 10 yields accuracy to the 4th
           significant digit in tau.
           Default: 10
        observing_redshift : optional, float
           This is the redshift at which the observer is observing
           the absorption spectrum.
           Default: 0
        njobs : optional, int or "auto"
           the number of process groups into which the loop over
           absorption lines will be divided.  If set to -1, each
           absorption line will be deposited by exactly one processor.
           If njobs is set to a value less than the total number of
           available processors (N), then the deposition of an
           individual line will be parallelized over (N / njobs)
           processors.  If set to "auto", it will first try to
           parallelize over the list of lines and only parallelize
           the line deposition if there are more processors than
           lines.  This is the optimal strategy for parallelizing
           spectrum generation.
           Default: "auto"
        """
        if line_list_file is not None:
            mylog.info("'line_list_file' keyword is deprecated. Please use " \
                       "'output_absorbers_file'.")
            output_absorbers_file = line_list_file

        input_fields = ['dl', 'redshift', 'temperature']
        field_units = {"dl": "cm", "redshift": "", "temperature": "K"}
        if use_peculiar_velocity:
            input_fields.append('velocity_los')
            input_fields.append('redshift_eff')
            field_units["velocity_los"] = "cm/s"
            field_units["redshift_eff"] = ""
        if observing_redshift != 0.:
            input_fields.append('redshift_dopp')
            field_units["redshift_dopp"] = ""
        for feature in self.line_list + self.continuum_list:
            if not feature['field_name'] in input_fields:
                input_fields.append(feature['field_name'])
                field_units[feature["field_name"]] = "cm**-3"

        if isinstance(input_file, string_types):
            input_ds = load(input_file)
        else:
            input_ds = input_file
        field_data = input_ds.all_data()

        # temperature field required to calculate voigt profile widths
        if ('temperature' not in input_ds.derived_field_list) and \
           (('gas', 'temperature') not in input_ds.derived_field_list):
            raise RuntimeError(
                "('gas', 'temperature') field required to be present in %s "
                "for AbsorptionSpectrum to function." % input_file)

        self.tau_field = np.zeros(self.lambda_field.size)
        self.absorbers_list = []

        if njobs == "auto":
            comm = _get_comm(())
            njobs = min(comm.size, len(self.line_list))

        mylog.info("Creating spectrum")
        self._add_lines_to_spectrum(field_data, use_peculiar_velocity,
                                    output_absorbers_file,
                                    subgrid_resolution=subgrid_resolution,
                                    observing_redshift=observing_redshift,
                                    njobs=njobs)
        self._add_continua_to_spectrum(field_data, use_peculiar_velocity,
                                       observing_redshift=observing_redshift)

        self.flux_field = np.exp(-self.tau_field)

        if output_file is None:
            pass
        elif output_file.endswith('.h5'):
            self._write_spectrum_hdf5(output_file)
        elif output_file.endswith('.fits'):
            self._write_spectrum_fits(output_file)
        else:
            self._write_spectrum_ascii(output_file)
        if output_absorbers_file is not None:
            self._write_absorbers_file(output_absorbers_file)

        del field_data
        return (self.lambda_field, self.flux_field)

    def _apply_observing_redshift(self, field_data, use_peculiar_velocity,
                                 observing_redshift):
        """
        Change the redshifts of individual absorbers to account for the
        redshift at which the observer sits.

        The intermediate redshift that is seen by an observer
        at a redshift other than z=0 is z12, where z1 is the
        observing redshift and z2 is the emitted photon's redshift
        Hogg (2000) eq. 13:

        1 + z12 = (1 + z2) / (1 + z1)
        """
        if observing_redshift == 0.:
            # This is already assumed in the generation of the LightRay
            redshift = field_data['redshift']
            if use_peculiar_velocity:
                redshift_eff = field_data['redshift_eff']
        else:
            # The intermediate redshift that is seen by an observer
            # at a redshift other than z=0 is z12, where z1 is the
            # observing redshift and z2 is the emitted photon's redshift
            # Hogg (2000) eq. 13:
            # 1 + z12 = (1 + z2) / (1 + z1)
            redshift = ((1 + field_data['redshift']) / \
                        (1 + observing_redshift)) - 1.
            # Combining cosmological redshift and doppler redshift
            # into an effective redshift is found in Peacock's
            # Cosmological Physics eqn 3.75:
            # 1 + z_eff = (1 + z_cosmo) * (1 + z_doppler)
            if use_peculiar_velocity:
                redshift_eff = ((1 + redshift) * \
                                (1 + field_data['redshift_dopp'])) - 1.

        if not use_peculiar_velocity:
            redshift_eff = redshift

        return redshift, redshift_eff

    def _add_continua_to_spectrum(self, field_data, use_peculiar_velocity,
                                  observing_redshift=0.):
        """
        Add continuum features to the spectrum.  Continuua are recorded as
        a name, associated field, wavelength, normalization value, and index.
        Continuua are applied at and below the denoted wavelength, where the
        optical depth decreases as a power law of desired index.  For positive 
        index values, this means optical depth is highest at the denoted 
        wavelength, and it drops with shorter and shorter wavelengths.  
        Consequently, transmitted flux undergoes a discontinuous cutoff at the 
        denoted wavelength, and then slowly increases with decreasing wavelength 
        according to the power law.
        """
        # Change the redshifts of continuum sources to account for the
        # redshift at which the observer sits
        redshift, redshift_eff = self._apply_observing_redshift(field_data,
                                 use_peculiar_velocity, observing_redshift)

        # min_tau is the minimum optical depth value that warrants 
        # accounting for an absorber.  for a single absorber, noticeable 
        # continuum effects begin for tau = 1e-3 (leading to transmitted 
        # flux of e^-tau ~ 0.999).  but we apply a cutoff to remove
        # absorbers with insufficient column_density to contribute 
        # significantly to a continuum (see below).  because lots of 
        # low column density absorbers can add up to a significant
        # continuum effect, we normalize min_tau by the n_absorbers.
        n_absorbers = field_data['dl'].size
        min_tau = 1.e-3/n_absorbers

        for continuum in self.continuum_list:

            # Normalization is in cm**-2, so column density must be as well
            column_density = (field_data[continuum['field_name']] * 
                              field_data['dl']).in_units('cm**-2')
            if (column_density == 0).all():
                mylog.info("Not adding continuum %s: insufficient column density" % continuum['label'])
                continue

            # redshift_eff field combines cosmological and velocity redshifts
            if use_peculiar_velocity:
                delta_lambda = continuum['wavelength'] * redshift_eff
            else:
                delta_lambda = continuum['wavelength'] * redshift

            # right index of continuum affected area is wavelength itself
            this_wavelength = delta_lambda + continuum['wavelength']
            right_index = np.digitize(this_wavelength, 
                                      self.lambda_field).clip(0, self.n_lambda)
            # left index of continuum affected area wavelength at which 
            # optical depth reaches tau_min
            left_index = np.digitize((this_wavelength *
                              np.power((min_tau * continuum['normalization'] /
                                        column_density),
                                       (1. / continuum['index']))),
                              self.lambda_field).clip(0, self.n_lambda)

            # Only calculate the effects of continuua where normalized 
            # column_density is greater than min_tau
            # because lower column will not have significant contribution
            valid_continuua = np.where(((column_density /
                                         continuum['normalization']) > min_tau) &
                                       (right_index - left_index > 1))[0]
            if valid_continuua.size == 0:
                mylog.info("Not adding continuum %s: insufficient column density or out of range" %
                    continuum['label'])
                continue

            pbar = get_pbar("Adding continuum - %s [%f A]: " % \
                                (continuum['label'], continuum['wavelength']),
                            valid_continuua.size)

            # Tau value is (wavelength / continuum_wavelength)**index / 
            #              (column_dens / norm)
            # i.e. a power law decreasing as wavelength decreases

            # Step through the absorber list and add continuum tau for each to
            # the total optical depth for all wavelengths
            for i, lixel in enumerate(valid_continuua):
                cont_tau = \
                    np.power((self.lambda_field[left_index[lixel] :
                                                right_index[lixel]] /
                                   this_wavelength[lixel]), \
                              continuum['index']) * \
                    (column_density[lixel] / continuum['normalization'])
                self.tau_field[left_index[lixel]:right_index[lixel]] += cont_tau
                pbar.update(i)
            pbar.finish()

    def _add_lines_to_spectrum(self, field_data, use_peculiar_velocity,
                               output_absorbers_file, subgrid_resolution=10,
                               observing_redshift=0., njobs=-1):
        """
        Add the absorption lines to the spectrum.
        """

        # Change the redshifts of individual absorbers to account for the
        # redshift at which the observer sits
        redshift, redshift_eff = self._apply_observing_redshift(field_data,
                                 use_peculiar_velocity, observing_redshift)

        # Widen wavelength window until optical depth falls below this tau
        # value at the ends to assure that the wings of a line have been
        # fully resolved.
        min_tau = 1e-3

        # step through each ionic transition (e.g. HI, HII, MgII) specified
        # and deposit the lines into the spectrum
        for line in parallel_objects(self.line_list, njobs=njobs):
            column_density = field_data[line['field_name']] * field_data['dl']
            if (column_density < 0).any():
                mylog.warn("Setting negative densities for field %s to 0! Bad!" % line['field_name'])
                np.clip(column_density, 0, np.inf, out=column_density)
            if (column_density == 0).all():
                mylog.info("Not adding line %s: insufficient column density" % line['label'])
                continue

            # redshift_eff field combines cosmological and velocity redshifts
            # so delta_lambda gives the offset in angstroms from the rest frame
            # wavelength to the observed wavelength of the transition
            if use_peculiar_velocity:
                delta_lambda = line['wavelength'] * redshift_eff
            else:
                delta_lambda = line['wavelength'] * redshift
            # lambda_obs is central wavelength of line after redshift
            lambda_obs = line['wavelength'] + delta_lambda
            # the total number of absorbers per transition
            n_absorbers = len(lambda_obs)

            # we want to know the bin index in the lambda_field array
            # where each line has its central wavelength after being
            # redshifted.  however, because we don't know a priori how wide
            # a line will be (ie DLAs), we have to include bin indices
            # *outside* the spectral range of the AbsorptionSpectrum
            # object.  Thus, we find the "equivalent" bin index, which
            # may be <0 or >the size of the array.  In the end, we deposit
            # the bins that actually overlap with the AbsorptionSpectrum's
            # range in lambda.

            # this equation gives us the "equivalent" bin index for each line
            # if it were placed into the self.lambda_field array
            center_index = (lambda_obs.in_units('Angstrom').d - self.lambda_min) \
                            / self.bin_width.d
            center_index = np.ceil(center_index).astype('int')

            # thermal broadening b parameter
            thermal_b =  np.sqrt((2 * boltzmann_constant_cgs *
                                  field_data['temperature']) /
                                  line['atomic_mass'])

            # the actual thermal width of the lines
            thermal_width = (lambda_obs * thermal_b /
                             speed_of_light_cgs).convert_to_units("angstrom")

            # Sanitize units for faster runtime of the tau_profile machinery.
            lambda_0 = line['wavelength'].d  # line's rest frame; angstroms
            cdens = column_density.in_units("cm**-2").d # cm**-2
            thermb = thermal_b.in_cgs().d  # thermal b coefficient; cm / s
            dlambda = delta_lambda.d  # lambda offset; angstroms
            if use_peculiar_velocity:
                vlos = field_data['velocity_los'].in_units("km/s").d # km/s
            else:
                vlos = np.zeros(field_data['temperature'].size)

            # When we actually deposit the voigt profile, sometimes we will
            # have underresolved lines (ie lines with smaller widths than
            # the spectral bin size).  Here, we create virtual wavelength bins
            # small enough in width to well resolve each line, deposit the
            # voigt profile into them, then numerically integrate their tau
            # values and sum them to redeposit them into the actual spectral
            # bins.

            # virtual bins (vbins) will be:
            # 1) <= the bin_width; assures at least as good as spectral bins
            # 2) <= 1/10th the thermal width; assures resolving voigt profiles
            #   (actually 1/subgrid_resolution value, default is 1/10)
            # 3) a bin width will be divisible by vbin_width times a power of
            #    10; this will assure we don't get spikes in the deposited
            #    spectra from uneven numbers of vbins per bin
            resolution = thermal_width / self.bin_width
            n_vbins_per_bin = (10 ** (np.ceil( np.log10( subgrid_resolution / 
                               resolution) ).clip(0, np.inf) ) ).astype('int')
            vbin_width = self.bin_width.d / n_vbins_per_bin

            # a note to the user about which lines components are unresolved
            if (thermal_width < self.bin_width).any():
                mylog.info("%d out of %d line components will be " +
                            "deposited as unresolved lines.",
                            (thermal_width < self.bin_width).sum(),
                            n_absorbers)

            # provide a progress bar with information about lines processsed
            pbar = get_pbar("Adding line - %s [%f A]: " % \
                            (line['label'], line['wavelength']), n_absorbers)

            # for a given transition, step through each location in the
            # observed spectrum where it occurs and deposit a voigt profile
            for i in parallel_objects(np.arange(n_absorbers), njobs=-1):

                # if there is a ray element with temperature = 0 or column
                # density = 0, skip it
                if (thermal_b[i] == 0.) or (cdens[i] == 0.):
                    pbar.update(i)
                    continue

                # the virtual window into which the line is deposited initially
                # spans a region of 2 coarse spectral bins
                # (one on each side of the center_index) but the window
                # can expand as necessary.
                # it will continue to expand until the tau value in the far
                # edge of the wings is less than the min_tau value or it
                # reaches the edge of the spectrum
                window_width_in_bins = 2

                while True:
                    left_index = (center_index[i] - window_width_in_bins/2)
                    right_index = (center_index[i] + window_width_in_bins/2)
                    n_vbins = (right_index - left_index) * n_vbins_per_bin[i]

                    # the array of virtual bins in lambda space
                    vbins = \
                        np.linspace(self.lambda_min + self.bin_width.d * left_index,
                                    self.lambda_min + self.bin_width.d * right_index,
                                    n_vbins, endpoint=False)

                    # the virtual bins and their corresponding opacities
                    vbins, vtau = \
                        tau_profile(
                            lambda_0, line['f_value'], line['gamma'],
                            thermb[i], cdens[i],
                            delta_lambda=dlambda[i], lambda_bins=vbins)

                    # If tau has not dropped below min tau threshold by the
                    # edges (ie the wings), then widen the wavelength
                    # window and repeat process.
                    if (vtau[0] < min_tau and vtau[-1] < min_tau):
                        break
                    window_width_in_bins *= 2

                # numerically integrate the virtual bins to calculate a
                # virtual equivalent width; then sum the virtual equivalent
                # widths and deposit into each spectral bin
                vEW = vtau * vbin_width[i]
                EW = np.zeros(right_index - left_index)
                EW_indices = np.arange(left_index, right_index)
                for k, val in enumerate(EW_indices):
                    EW[k] = vEW[n_vbins_per_bin[i] * k: \
                                n_vbins_per_bin[i] * (k + 1)].sum()
                EW = EW/self.bin_width.d

                # only deposit EW bins that actually intersect the original
                # spectral wavelength range (i.e. lambda_field)

                # if EW bins don't intersect the original spectral range at all
                # then skip the deposition
                if ((left_index >= self.n_lambda) or \
                    (right_index < 0)):
                    pbar.update(i)
                    continue

                # otherwise, determine how much of the original spectrum
                # is intersected by the expanded line window to be deposited,
                # and deposit the Equivalent Width data into that intersecting
                # window in the original spectrum's tau
                else:
                    intersect_left_index = max(left_index, 0)
                    intersect_right_index = min(right_index, self.n_lambda-1)
                    self.tau_field[intersect_left_index:intersect_right_index] \
                        += EW[(intersect_left_index - left_index): \
                              (intersect_right_index - left_index)]


                # write out absorbers to file if the column density of
                # an absorber is greater than the specified "label_threshold"
                # of that absorption line
                if output_absorbers_file and \
                   line['label_threshold'] is not None and \
                   cdens[i] >= line['label_threshold']:

                    if use_peculiar_velocity:
                        peculiar_velocity = vlos[i]
                    else:
                        peculiar_velocity = 0.0
                    self.absorbers_list.append({'label': line['label'],
                                                'wavelength': (lambda_0 + dlambda[i]),
                                                'column_density': column_density[i],
                                                'b_thermal': thermal_b[i],
                                                'redshift': redshift[i],
                                                'redshift_eff': redshift_eff[i],
                                                'v_pec': peculiar_velocity})
                pbar.update(i)
            pbar.finish()

            del column_density, delta_lambda, lambda_obs, center_index, \
                thermal_b, thermal_width, cdens, thermb, dlambda, \
                vlos, resolution, vbin_width, n_vbins, n_vbins_per_bin

        comm = _get_comm(())
        self.tau_field = comm.mpi_allreduce(self.tau_field, op="sum")
        if output_absorbers_file:
            self.absorbers_list = comm.par_combine_object(
                self.absorbers_list, "cat", datatype="list")

    @parallel_root_only
    def _write_absorbers_file(self, filename):
        """
        Write out ASCII list of all substantial absorbers found in spectrum
        """
        if filename is None:
            return
        mylog.info("Writing absorber list: %s.", filename)
        self.absorbers_list.sort(key=lambda obj: obj['wavelength'])
        f = open(filename, 'w')
        f.write('#%-14s %-14s %-12s %-14s %-15s %-9s %-10s\n' %
                ('Wavelength', 'Line', 'N [cm^-2]', 'b [km/s]', 'z_cosmo', \
                 'z_eff', 'v_pec [km/s]'))
        for line in self.absorbers_list:
            f.write('%-14.6f %-14ls %e %e % e % e % e\n' % (line['wavelength'], \
                line['label'], line['column_density'], line['b_thermal'], \
                line['redshift'], line['redshift_eff'], line['v_pec']))
        f.close()

    @parallel_root_only
    def _write_spectrum_ascii(self, filename):
        """
        Write spectrum to an ascii file.
        """
        mylog.info("Writing spectrum to ascii file: %s.", filename)
        f = open(filename, 'w')
        f.write("# wavelength[A] tau flux\n")
        for i in range(self.lambda_field.size):
            f.write("%e %e %e\n" % (self.lambda_field[i],
                                    self.tau_field[i], self.flux_field[i]))
        f.close()

    @parallel_root_only
    def _write_spectrum_fits(self, filename):
        """
        Write spectrum to a fits file.
        """
        mylog.info("Writing spectrum to fits file: %s.", filename)
        col1 = pyfits.Column(name='wavelength', format='E', array=self.lambda_field)
        col2 = pyfits.Column(name='tau', format='E', array=self.tau_field)
        col3 = pyfits.Column(name='flux', format='E', array=self.flux_field)
        cols = pyfits.ColDefs([col1, col2, col3])
        tbhdu = pyfits.BinTableHDU.from_columns(cols)
        tbhdu.writeto(filename, clobber=True)

    @parallel_root_only
    def _write_spectrum_hdf5(self, filename):
        """
        Write spectrum to an hdf5 file.

        """
        mylog.info("Writing spectrum to hdf5 file: %s.", filename)
        output = h5py.File(filename, 'w')
        output.create_dataset('wavelength', data=self.lambda_field)
        output.create_dataset('tau', data=self.tau_field)
        output.create_dataset('flux', data=self.flux_field)
        output.close()
