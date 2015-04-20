"""
CosmologyTimeSeries class and member functions.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.convenience import \
    simulation
from yt.funcs import *
from yt.utilities.cosmology import \
    Cosmology

class CosmologySplice(object):
    """
    Class for splicing together datasets to extend over a
    cosmological distance.
    """

    def __init__(self, parameter_filename, simulation_type, find_outputs=False):
        self.parameter_filename = parameter_filename
        self.simulation_type = simulation_type
        self.simulation = simulation(parameter_filename, simulation_type, 
                                     find_outputs=find_outputs)

        self.cosmology = Cosmology(
            hubble_constant=(self.simulation.hubble_constant),
            omega_matter=self.simulation.omega_matter,
            omega_lambda=self.simulation.omega_lambda)

    def create_cosmology_splice(self, near_redshift, far_redshift,
                                minimal=True, deltaz_min=0.0,
                                time_data=True, redshift_data=True):
        r"""Create list of datasets capable of spanning a redshift
        interval.

        For cosmological simulations, the physical width of the simulation
        box corresponds to some \Delta z, which varies with redshift.
        Using this logic, one can stitch together a series of datasets to
        create a continuous volume or length element from one redshift to
        another. This method will return such a list

        Parameters
        ----------
        near_redshift : float
            The nearest (lowest) redshift in the cosmology splice list.
        far_redshift : float
            The furthest (highest) redshift in the cosmology splice list.
        minimal : bool
            If True, the minimum number of datasets is used to connect the
            initial and final redshift.  If false, the list will contain as
            many entries as possible within the redshift
            interval.
            Default: True.
        deltaz_min : float
            Specifies the minimum delta z between consecutive datasets
            in the returned
            list.
            Default: 0.0.
        time_data : bool
            Whether or not to include time outputs when gathering
            datasets for time series.
            Default: True.
        redshift_data : bool
            Whether or not to include redshift outputs when gathering
            datasets for time series.
            Default: True.

        Examples
        --------

        >>> co = CosmologySplice("enzo_tiny_cosmology/32Mpc_32.enzo", "Enzo")
        >>> cosmo = co.create_cosmology_splice(1.0, 0.0)

        """

        if time_data and redshift_data:
            self.splice_outputs = self.simulation.all_outputs
        elif time_data:
            self.splice_outputs = self.simulation.all_time_outputs
        elif redshift_data:
            self.splice_outputs = self.simulation.all_redshift_outputs
        else:
            mylog.error('Both time_data and redshift_data are False.')
            return

        # Link datasets in list with pointers.
        # This is used for connecting datasets together.
        for i, output in enumerate(self.splice_outputs):
            if i == 0:
                output['previous'] = None
                output['next'] = self.splice_outputs[i + 1]
            elif i == len(self.splice_outputs) - 1:
                output['previous'] = self.splice_outputs[i - 1]
                output['next'] = None
            else:
                output['previous'] = self.splice_outputs[i - 1]
                output['next'] = self.splice_outputs[i + 1]

        # Calculate maximum delta z for each data dump.
        self._calculate_deltaz_max()

        # Calculate minimum delta z for each data dump.
        self._calculate_deltaz_min(deltaz_min=deltaz_min)

        cosmology_splice = []
 
        if near_redshift == far_redshift:
            self.simulation.get_time_series(redshifts=[near_redshift])
            cosmology_splice.append({'time': self.simulation[0].current_time,
                                     'redshift': self.simulation[0].current_redshift,
                                     'filename': os.path.join(self.simulation[0].fullpath,
                                                              self.simulation[0].basename),
                                     'next': None})
            mylog.info("create_cosmology_splice: Using %s for z = %f ." %
                       (cosmology_splice[0]['filename'], near_redshift))
            return cosmology_splice
        
        # Use minimum number of datasets to go from z_i to z_f.
        if minimal:

            z_Tolerance = 1e-3
            z = far_redshift

            # fill redshift space with datasets
            while ((z > near_redshift) and
                   (np.abs(z - near_redshift) > z_Tolerance)):

                # For first data dump, choose closest to desired redshift.
                if (len(cosmology_splice) == 0):
                    # Sort data outputs by proximity to current redshift.
                    self.splice_outputs.sort(key=lambda obj:np.fabs(z - \
                        obj['redshift']))
                    cosmology_splice.append(self.splice_outputs[0])

                # Move forward from last slice in stack until z > z_max.
                else:
                    current_slice = cosmology_splice[-1]
                    while current_slice['next'] is not None and \
                            (z < current_slice['next']['redshift'] or \
                                 np.abs(z - current_slice['next']['redshift']) <
                                 z_Tolerance):
                        current_slice = current_slice['next']

                    if current_slice is cosmology_splice[-1]:
                        near_redshift = cosmology_splice[-1]['redshift'] - \
                          cosmology_splice[-1]['dz_max']
                        mylog.error("Cosmology splice incomplete due to insufficient data outputs.")
                        break
                    else:
                        cosmology_splice.append(current_slice)

                z = cosmology_splice[-1]['redshift'] - \
                  cosmology_splice[-1]['dz_max']

        # Make light ray using maximum number of datasets (minimum spacing).
        else:
            # Sort data outputs by proximity to current redsfhit.
            self.splice_outputs.sort(key=lambda obj:np.abs(far_redshift -
                                                           obj['redshift']))
            # For first data dump, choose closest to desired redshift.
            cosmology_splice.append(self.splice_outputs[0])

            nextOutput = cosmology_splice[-1]['next']
            while (nextOutput is not None):
                if (nextOutput['redshift'] <= near_redshift):
                    break
                if ((cosmology_splice[-1]['redshift'] - nextOutput['redshift']) >
                    cosmology_splice[-1]['dz_min']):
                    cosmology_splice.append(nextOutput)
                nextOutput = nextOutput['next']
            if (cosmology_splice[-1]['redshift'] -
                cosmology_splice[-1]['dz_max']) > near_redshift:
                mylog.error("Cosmology splice incomplete due to insufficient data outputs.")
                near_redshift = cosmology_splice[-1]['redshift'] - \
                  cosmology_splice[-1]['dz_max']

        mylog.info("create_cosmology_splice: Used %d data dumps to get from z = %f to %f." %
                   (len(cosmology_splice), far_redshift, near_redshift))
        
        # change the 'next' and 'previous' pointers to point to the correct outputs for the created
        # splice
        for i, output in enumerate(cosmology_splice):
            if len(cosmology_splice) == 1:
                output['previous'] = None
                output['next'] = None
            elif i == 0:
                output['previous'] = None
                output['next'] = cosmology_splice[i + 1]
            elif i == len(cosmology_splice) - 1:
                output['previous'] = cosmology_splice[i - 1]
                output['next'] = None
            else:
                output['previous'] = cosmology_splice[i - 1]
                output['next'] = cosmology_splice[i + 1]
        
        self.splice_outputs.sort(key=lambda obj: obj['time'])
        return cosmology_splice

    def plan_cosmology_splice(self, near_redshift, far_redshift,
                              decimals=3, filename=None,
                              start_index=0):
        r"""Create imaginary list of redshift outputs to maximally
        span a redshift interval.

        If you want to run a cosmological simulation that will have just
        enough data outputs to create a cosmology splice,
        this method will calculate a list of redshifts outputs that will
        minimally connect a redshift interval.

        Parameters
        ----------
        near_redshift : float
            The nearest (lowest) redshift in the cosmology splice list.
        far_redshift : float
            The furthest (highest) redshift in the cosmology splice list.
        decimals : int
            The decimal place to which the output redshift will be rounded.
            If the decimal place in question is nonzero, the redshift will
            be rounded up to
            ensure continuity of the splice.  Default: 3.
        filename : string
            If provided, a file will be written with the redshift outputs in
            the form in which they should be given in the enzo dataset.
            Default: None.
        start_index : int
            The index of the first redshift output.  Default: 0.

        Examples
        --------
        >>> from yt.analysis_modules.api import CosmologySplice
        >>> my_splice = CosmologySplice('enzo_tiny_cosmology/32Mpc_32.enzo', 'Enzo')
        >>> my_splice.plan_cosmology_splice(0.0, 0.1, filename='redshifts.out')

        """

        z = far_redshift
        outputs = []

        while z > near_redshift:
            rounded = np.round(z, decimals=decimals)
            if rounded - z < 0:
                rounded += np.power(10.0, (-1.0*decimals))
            z = rounded

            deltaz_max = self._deltaz_forward(z, self.simulation.box_size)
            outputs.append({'redshift': z, 'dz_max': deltaz_max})
            z -= deltaz_max

        mylog.info("%d data dumps will be needed to get from z = %f to %f." %
                   (len(outputs), near_redshift, far_redshift))

        if filename is not None:
            self.simulation._write_cosmology_outputs(filename, outputs,
                                                     start_index,
                                                     decimals=decimals)
        return outputs

    def _calculate_deltaz_max(self):
        r"""Calculate delta z that corresponds to full box length going
        from z to (z - delta z).
        """

        d_Tolerance = 1e-4
        max_Iterations = 100

        target_distance = self.simulation.box_size

        for output in self.splice_outputs:
            z = output['redshift']

            # Calculate delta z that corresponds to the length of the box
            # at a given redshift using Newton's method.
            z1 = z
            z2 = z1 - 0.1 # just an initial guess
            distance1 = self.simulation.quan(0.0, "Mpccm / h")
            distance2 = self.cosmology.comoving_radial_distance(z2, z)
            iteration = 1

            while ((np.abs(distance2-target_distance)/distance2) > d_Tolerance):
                m = (distance2 - distance1) / (z2 - z1)
                z1 = z2
                distance1 = distance2
                z2 = ((target_distance - distance2) / m.in_units("Mpccm / h")) + z2
                distance2 = self.cosmology.comoving_radial_distance(z2, z)
                iteration += 1
                if (iteration > max_Iterations):
                    mylog.error("calculate_deltaz_max: Warning - max iterations " +
                                "exceeded for z = %f (delta z = %f)." %
                                (z, np.abs(z2 - z)))
                    break
            output['dz_max'] = np.abs(z2 - z)
            
    def _calculate_deltaz_min(self, deltaz_min=0.0):
        r"""Calculate delta z that corresponds to a single top grid pixel
        going from z to (z - delta z).
        """

        d_Tolerance = 1e-4
        max_Iterations = 100

        target_distance = self.simulation.box_size / \
          self.simulation.domain_dimensions[0]

        for output in self.splice_outputs:
            z = output['redshift']

            # Calculate delta z that corresponds to the length of a
            # top grid pixel at a given redshift using Newton's method.
            z1 = z
            z2 = z1 - 0.01 # just an initial guess
            distance1 = self.simulation.quan(0.0, "Mpccm / h")
            distance2 = self.cosmology.comoving_radial_distance(z2, z)
            iteration = 1

            while ((np.abs(distance2 - target_distance) / distance2) > d_Tolerance):
                m = (distance2 - distance1) / (z2 - z1)
                z1 = z2
                distance1 = distance2
                z2 = ((target_distance - distance2) / m.in_units("Mpccm / h")) + z2
                distance2 = self.cosmology.comoving_radial_distance(z2, z)
                iteration += 1
                if (iteration > max_Iterations):
                    mylog.error("calculate_deltaz_max: Warning - max iterations " +
                                "exceeded for z = %f (delta z = %f)." %
                                (z, np.abs(z2 - z)))
                    break
            # Use this calculation or the absolute minimum specified by the user.
            output['dz_min'] = max(np.abs(z2 - z), deltaz_min)

    def _deltaz_forward(self, z, target_distance):
        r"""Calculate deltaz corresponding to moving a comoving distance
        starting from some redshift.
        """

        d_Tolerance = 1e-4
        max_Iterations = 100

        # Calculate delta z that corresponds to the length of the
        # box at a given redshift.
        z1 = z
        z2 = z1 - 0.1 # just an initial guess
        distance1 = self.cosmology.quan(0.0, "Mpccm / h")
        distance2 = self.cosmology.comoving_radial_distance(z2, z)
        iteration = 1

        while ((np.abs(distance2 - target_distance)/distance2) > d_Tolerance):
            m = (distance2 - distance1) / (z2 - z1)
            z1 = z2
            distance1 = distance2
            z2 = ((target_distance - distance2) / m.in_units("Mpccm / h")) + z2
            distance2 = self.cosmology.comoving_radial_distance(z2, z)
            iteration += 1
            if (iteration > max_Iterations):
                mylog.error("deltaz_forward: Warning - max iterations " +
                            "exceeded for z = %f (delta z = %f)." %
                            (z, np.abs(z2 - z)))
                break
        return np.abs(z2 - z)
