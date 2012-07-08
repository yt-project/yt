"""
CosmologyTimeSeries class and member functions.

Author: Britton Smith <brittonsmith@gmail.com>
Affiliation: Michigan State University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008-2012 Britton Smith.  All Rights Reserved.

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

import numpy as na

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

    def __init__(self, parameter_filename, simulation_type):
        self.parameter_filename = parameter_filename
        self.simulation_type = simulation_type
        self.simulation = simulation(parameter_filename, simulation_type)

        self.cosmology = Cosmology(
            HubbleConstantNow=(100.0 * self.simulation.hubble_constant),
            OmegaMatterNow=self.simulation.omega_matter,
            OmegaLambdaNow=self.simulation.omega_lambda)

    def create_cosmology_splice(self, near_redshift, far_redshift,
                                minimal=True, deltaz_min=0.0):
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

        Examples
        --------
        >>> cosmo = es.create_cosmology_splice(1.0, 0.0, minimal=True,
                                               deltaz_min=0.0)

        """

        # Link datasets in list with pointers.
        # This is used for connecting datasets together.
        for i, output in enumerate(self.simulation.all_outputs):
            if i == 0:
                output['previous'] = None
                output['next'] = self.simulation.all_outputs[i + 1]
            elif i == len(self.simulation.all_outputs) - 1:
                output['previous'] = self.simulation.all_outputs[i - 1]
                output['next'] = None
            else:
                output['previous'] = self.simulation.all_outputs[i - 1]
                output['next'] = self.simulation.all_outputs[i + 1]

        # Calculate maximum delta z for each data dump.
        self._calculate_deltaz_max()

        # Calculate minimum delta z for each data dump.
        self._calculate_deltaz_min(deltaz_min=deltaz_min)

        cosmology_splice = []

        # Use minimum number of datasets to go from z_i to z_f.
        if minimal:

            z_Tolerance = 1e-3
            z = far_redshift

            # fill redshift space with datasets
            while ((z > near_redshift) and
                   (na.fabs(z - near_redshift) > z_Tolerance)):

                # For first data dump, choose closest to desired redshift.
                if (len(cosmology_splice) == 0):
                    # Sort data outputs by proximity to current redsfhit.
                    self.simulation.all_outputs.sort(key=lambda obj:na.fabs(z - \
                        obj['redshift']))
                    cosmology_splice.append(self.simulation.all_outputs[0])

                # Move forward from last slice in stack until z > z_max.
                else:
                    current_slice = cosmology_splice[-1]
                    while current_slice['next'] is not None and \
                            (z < current_slice['next']['redshift'] or \
                                 na.abs(z - current_slice['next']['redshift']) <
                                 z_Tolerance):
                        current_slice = current_slice['next']

                    if current_slice is cosmology_splice[-1]:
                        near_redshift = cosmology_splice[-1]['redshift'] - \
                          cosmology_splice[-1]['deltazMax']
                        mylog.error("Cosmology splice incomplete due to insufficient data outputs.")
                        break
                    else:
                        cosmology_splice.append(current_slice)

                z = cosmology_splice[-1]['redshift'] - \
                  cosmology_splice[-1]['deltazMax']

        # Make light ray using maximum number of datasets (minimum spacing).
        else:
            # Sort data outputs by proximity to current redsfhit.
            self.simulation.all_outputs.sort(key=lambda obj:na.fabs(far_redshift -
                                                                    obj['redshift']))
            # For first data dump, choose closest to desired redshift.
            cosmology_splice.append(self.simulation.all_outputs[0])

            nextOutput = cosmology_splice[-1]['next']
            while (nextOutput is not None):
                if (nextOutput['redshift'] <= near_redshift):
                    break
                if ((cosmology_splice[-1]['redshift'] - nextOutput['redshift']) >
                    cosmology_splice[-1]['deltazMin']):
                    cosmology_splice.append(nextOutput)
                nextOutput = nextOutput['next']
            if (cosmology_splice[-1]['redshift'] -
                cosmology_splice[-1]['deltazMax']) > near_redshift:
                mylog.error("Cosmology splice incomplete due to insufficient data outputs.")
                near_redshift = cosmology_splice[-1]['redshift'] - \
                  cosmology_splice[-1]['deltazMax']

        mylog.info("create_cosmology_splice: Used %d data dumps to get from z = %f to %f." %
                   (len(cosmology_splice), far_redshift, near_redshift))

        self.simulation.all_outputs.sort(key=lambda obj: obj['time'])
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
            the form in which they should be given in the enzo parameter file.
            Default: None.
        start_index : int
            The index of the first redshift output.  Default: 0.

        Examples
        --------
        >>> outputs = es.imagine_minimal_splice(0., 0.4,
                filename='outputs.out')

        """

        z = far_redshift
        outputs = []

        while z > near_redshift:
            rounded = na.round(z, decimals=decimals)
            if rounded - z < 0:
                rounded += na.power(10.0, (-1.0*decimals))
            z = rounded

            deltaz_max = self._deltaz_forward(z, self.simulation.box_size)
            outputs.append({'redshift': z, 'deltazMax': deltaz_max})
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

        for output in self.simulation.all_outputs:
            z = output['redshift']

            # Calculate delta z that corresponds to the length of the box
            # at a given redshift using Newton's method.
            z1 = z
            z2 = z1 - 0.1 # just an initial guess
            distance1 = 0.0
            iteration = 1

            # Convert comoving radial distance into Mpc / h,
            # since that's how box size is stored.
            distance2 = self.cosmology.ComovingRadialDistance(z2, z) * \
              self.simulation.hubble_constant

            while ((na.fabs(distance2-target_distance)/distance2) > d_Tolerance):
                m = (distance2 - distance1) / (z2 - z1)
                z1 = z2
                distance1 = distance2
                z2 = ((target_distance - distance2) / m) + z2
                distance2 = self.cosmology.ComovingRadialDistance(z2, z) * \
                  self.simulation.hubble_constant
                iteration += 1
                if (iteration > max_Iterations):
                    mylog.error("calculate_deltaz_max: Warning - max iterations exceeded for z = %f (delta z = %f)." %
                                (z, na.fabs(z2 - z)))
                    break
            output['deltazMax'] = na.fabs(z2 - z)

    def _calculate_deltaz_min(self, deltaz_min=0.0):
        r"""Calculate delta z that corresponds to a single top grid pixel
        going from z to (z - delta z).
        """

        d_Tolerance = 1e-4
        max_Iterations = 100

        target_distance = self.simulation.box_size / \
          self.simulation.domain_dimensions[0]

        for output in self.simulation.all_outputs:
            z = output['redshift']

            # Calculate delta z that corresponds to the length of a
            # top grid pixel at a given redshift using Newton's method.
            z1 = z
            z2 = z1 - 0.01 # just an initial guess
            distance1 = 0.0
            iteration = 1

            # Convert comoving radial distance into Mpc / h,
            # since that's how box size is stored.
            distance2 = self.cosmology.ComovingRadialDistance(z2, z) * \
              self.simulation.hubble_constant

            while ((na.fabs(distance2 - target_distance) / distance2) > d_Tolerance):
                m = (distance2 - distance1) / (z2 - z1)
                z1 = z2
                distance1 = distance2
                z2 = ((target_distance - distance2) / m) + z2
                distance2 = self.cosmology.ComovingRadialDistance(z2, z) * \
                  self.simulation.hubble_constant
                iteration += 1
                if (iteration > max_Iterations):
                    mylog.error("calculate_deltaz_max: Warning - max iterations exceeded for z = %f (delta z = %f)." %
                                (z, na.fabs(z2 - z)))
                    break
            # Use this calculation or the absolute minimum specified by the user.
            output['deltazMin'] = max(na.fabs(z2 - z), deltaz_min)

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
        distance1 = 0.0
        iteration = 1

        # Convert comoving radial distance into Mpc / h,
        # since that's how box size is stored.
        distance2 = self.cosmology.ComovingRadialDistance(z2, z) * \
          self.cosmology.HubbleConstantNow / 100.0

        while ((na.fabs(distance2 - target_distance)/distance2) > d_Tolerance):
            m = (distance2 - distance1) / (z2 - z1)
            z1 = z2
            distance1 = distance2
            z2 = ((target_distance - distance2) / m) + z2
            distance2 = self.cosmology.ComovingRadialDistance(z2, z) * \
              self.cosmology.HubbleConstantNow / 100.0
            iteration += 1
            if (iteration > max_Iterations):
                mylog.error("deltaz_forward: Warning - max iterations exceeded for z = %f (delta z = %f)." %
                            (z, na.fabs(z2 - z)))
                break
        return na.fabs(z2 - z)
