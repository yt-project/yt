"""
EnzoSimulation class and member functions.

Author: Britton Smith <brittons@origins.colorado.edu>
Affiliation: CASA/University of CO, Boulder
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008-2009 Britton Smith.  All Rights Reserved.

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

from yt.funcs import *

import numpy as na
import glob
import os

dt_Tolerance = 1e-3

from yt.utilities.cosmology import \
    Cosmology, \
    EnzoCosmology

from yt.convenience import \
    load

class EnzoSimulation(object):
    """
    Super class for performing the same operation over all data dumps in 
    a simulation from one redshift to another.
    """
    def __init__(self, EnzoParameterFile, initial_time=None, final_time=None, initial_redshift=None, final_redshift=None,
                 links=False, enzo_parameters=None, get_time_outputs=True, get_redshift_outputs=True, get_available_data=False,
                 get_data_by_force=False):
        """
        Initialize an EnzoSimulation object.
        :param initial_time (float): the initial time in code units for the dataset list.  Default: None.
        :param final_time (float): the final time in code units for the dataset list.  Default: None.
        :param initial_redshift (float): the initial (highest) redshift for the dataset list.  Only for 
               cosmological simulations.  Default: None.
        :param final_redshift (float): the final (lowest) redshift for the dataset list.  Only for cosmological 
               simulations.  Default: None.
        :param links (bool): if True, each entry in the dataset list will contain entries, previous and next, that 
               point to the previous and next entries on the dataset list.  Default: False.
        :param enzo_parameters (dict): a dictionary specify additional parameters to be retrieved from the 
               parameter file.  The format should be the name of the parameter as the key and the variable type as 
               the value.  For example, {'CosmologyComovingBoxSize':float}.  All parameter values will be stored in 
               the dictionary attribute, enzoParameters.  Default: None.
        :param get_time_outputs (bool): if False, the time datasets, specified in Enzo with the dtDataDump, will not 
               be added to the dataset list.  Default: True.
        :param get_redshift_outputs (bool): if False, the redshift datasets will not be added to the dataset list.  Default: True.
        :param get_available_data (bool): if True, only datasets that are found to exist at the file path are added 
               to the list.  Default: False.
        :param get_data_by_force (bool): if True, time data dumps are not calculated using dtDataDump.  Instead, the 
               the working directory is searched for directories that match the datadumpname keyword.  Each dataset 
               is loaded up to get the time and redshift manually.  This is useful with collapse simulations that use 
               OutputFirstTimeAtLevel or with simulations that make outputs based on cycle numbers.  Default: False.
        """
        self.EnzoParameterFile = EnzoParameterFile
        self.enzoParameters = {}
        self.redshift_outputs = []
        self.allOutputs = []
        self.InitialTime = initial_time
        self.FinalTime = final_time
        self.InitialRedshift = initial_redshift
        self.FinalRedshift = final_redshift
        self.links = links
        self.get_time_outputs = get_time_outputs
        self.get_redshift_outputs = get_redshift_outputs
        self.get_available_data = get_available_data

        # Add any extra parameters to parameter dict.
        if enzo_parameters is None: enzo_parameters = {}
        EnzoParameterDict.update(enzo_parameters)

        # Set some parameter defaults.
        self._set_parameter_defaults()

        # Read parameters.
        self._read_enzo_parameter_file()

        # Get all the appropriate datasets.
        self._get_all_outputs(brute_force=get_data_by_force)

    def _calculate_redshift_dump_times(self):
        "Calculates time from redshift of redshift dumps."

        for output in self.redshift_outputs:
            output['time'] = self.enzo_cosmology.ComputeTimeFromRedshift(output['redshift']) / \
                self.enzo_cosmology.TimeUnits

    def _calculate_time_dumps(self):
        "Calculates time dumps and their redshifts if cosmological."

        if self.enzoParameters['dtDataDump'] <= 0.0: return []
        time_outputs = []

        index = 0
        current_time = self.SimulationInitialTime
        while (current_time <= self.FinalTime) or \
                (abs(self.FinalTime - current_time) / self.FinalTime) < dt_Tolerance:
            filename = "%s/%s%04d/%s%04d" % (self.enzoParameters['GlobalDir'],
                                             self.enzoParameters['DataDumpDir'],index,
                                             self.enzoParameters['DataDumpName'],index)

            output = {'index':index,'filename':filename,'time':current_time}
            if self.enzoParameters['ComovingCoordinates']:
                t = self.enzo_cosmology.InitialTime + (self.enzoParameters['dtDataDump'] * self.enzo_cosmology.TimeUnits * index)
                output['redshift'] = self.enzo_cosmology.ComputeRedshiftFromTime(t)

            if not self.get_available_data or os.path.exists(filename):
                time_outputs.append(output)

            current_time += self.enzoParameters['dtDataDump']
            index += 1

        return time_outputs

    def _combine_data_outputs(self, brute_force=False):
        "Combines redshift and time data into one sorted list."

        # If True, get time dumps just by looking through the working directory.
        if brute_force:
            time_outputs = self._get_outputs_by_force()

        # Calculate time dumps based on dtDataDump
        elif self.enzoParameters.has_key('dtDataDump') and self.get_time_outputs:
            time_outputs = self._calculate_time_dumps()

        else:
            time_outputs =[]

        # Calculate times for redshift dumps.
        if self.enzoParameters['ComovingCoordinates'] and self.get_redshift_outputs:
            self._calculate_redshift_dump_times()
        else:
            self.redshift_outputs = []

        self.allOutputs = self.redshift_outputs + time_outputs
        self.allOutputs.sort(key=lambda obj:obj['time'])

        start_index = None
        end_index = None

        # Add links to next and previous dataset to each entry.
        for q in range(len(self.allOutputs)):
            if self.allOutputs[q].has_key('index'): del self.allOutputs[q]['index']

            if start_index is None:
                if self.allOutputs[q]['time'] >= self.InitialTime or \
                        abs(self.allOutputs[q]['time'] - self.InitialTime) < \
                        self.allOutputs[q]['time'] * dt_Tolerance:
                    start_index = q

            if end_index is None:
                if self.allOutputs[q]['time'] > self.FinalTime:
                    end_index = q - 1
                if abs(self.allOutputs[q]['time'] - self.FinalTime) < \
                        self.allOutputs[q]['time'] * dt_Tolerance:
                    end_index = q
                if q == len(self.allOutputs) - 1:
                    end_index = q

        for q in range(len(self.allOutputs)):
            if self.links and start_index is not None:
                if q <= start_index:
                    self.allOutputs[q]['previous'] = None
                else:
                    self.allOutputs[q]['previous'] = self.allOutputs[q-1]

                if q >= end_index:
                    self.allOutputs[q]['next'] = None
                else:
                    self.allOutputs[q]['next'] = self.allOutputs[q+1]

        del self.redshift_outputs

        if end_index is None:
            end_index = len(self.allOutputs)-1

        self.allOutputs = self.allOutputs[start_index:end_index+1]
        mylog.info("Loaded %s total data outputs." % len(self.allOutputs))

    def _get_all_outputs(self, brute_force=False):
        """
        Get all the datasets in the time/redshift interval requested, or search 
        the data directory for potential datasets.
        """

        # Check for sufficient starting/ending parameters.
        if self.InitialTime is None and self.InitialRedshift is None:
            if self.enzoParameters['ComovingCoordinates'] and \
               'CosmologyInitialRedshift' in self.enzoParameters:
                self.InitialRedshift = self.enzoParameters['CosmologyInitialRedshift']
            elif 'InitialTime' in self.enzoParameters:
                self.InitialTime = self.enzoParameters['InitialTime']
            else:
                mylog.error("Couldn't find parameter for initial time or redshift from parameter file.")
                return None

        if self.FinalTime is None and self.FinalRedshift is None:
            if self.enzoParameters['ComovingCoordinates'] and \
               'CosmologyFinalRedshift' in self.enzoParameters:
                self.FinalRedshift = self.enzoParameters['CosmologyFinalRedshift']
            elif 'StopTime' in self.enzoParameters:
                self.FinalTime = self.enzoParameters['StopTime']
            else:
                mylog.error("Couldn't find parameter for final time or redshift from parameter file.")
                return None

        # Convert initial/final redshifts to times.
        if self.enzoParameters['ComovingCoordinates']:
            # Instantiate a cosmology calculator.
            self.cosmology = Cosmology(HubbleConstantNow = 
                                       (100.0 * self.enzoParameters['CosmologyHubbleConstantNow']),
                                       OmegaMatterNow = self.enzoParameters['CosmologyOmegaMatterNow'],
                                       OmegaLambdaNow = self.enzoParameters['CosmologyOmegaLambdaNow'])

            # Instantiate EnzoCosmology object for units and time conversions.
            self.enzo_cosmology = EnzoCosmology(HubbleConstantNow = 
                                                (100.0 * self.enzoParameters['CosmologyHubbleConstantNow']),
                                                OmegaMatterNow = self.enzoParameters['CosmologyOmegaMatterNow'],
                                                OmegaLambdaNow = self.enzoParameters['CosmologyOmegaLambdaNow'],
                                                InitialRedshift = self.enzoParameters['CosmologyInitialRedshift'])
            if self.InitialRedshift is not None:
                self.InitialTime = self.enzo_cosmology.ComputeTimeFromRedshift(self.InitialRedshift) / \
                    self.enzo_cosmology.TimeUnits
            if self.FinalRedshift is not None:
                self.FinalTime = self.enzo_cosmology.ComputeTimeFromRedshift(self.FinalRedshift) / \
                    self.enzo_cosmology.TimeUnits

        # Get initial time of simulation.
        if self.enzoParameters['ComovingCoordinates'] and \
                'CosmologyInitialRedshift' in self.enzoParameters:
            self.SimulationInitialTime = self.enzo_cosmology.InitialTime / self.enzo_cosmology.TimeUnits
        elif 'InitialTime' in self.enzoParameters:
            self.SimulationInitialTime = self.enzoParameters['InitialTime']
        else:
            self.SimulationInitialTime = 0.0

        # Combine all data dumps in to ordered list.
        self._combine_data_outputs(brute_force=brute_force)

    def _get_outputs_by_force(self):
        """
        Search for directories matching the data dump keywords.
        If found, get dataset times by brute force py opening the pf.
        """

        # look for time dumps.
        potential_time_outputs = glob.glob("%s/%s*" % (self.enzoParameters['GlobalDir'], 
                                                       self.enzoParameters['DataDumpDir']))
        time_outputs = []
        mylog.info("Checking validity of %d potential time dumps." % 
                   len(potential_time_outputs))

        for output in potential_time_outputs:
            index = output[output.find(self.enzoParameters['DataDumpDir']) + 
                           len(self.enzoParameters['DataDumpDir']):]
            filename = "%s/%s%s/%s%s" % (self.enzoParameters['GlobalDir'],
                                         self.enzoParameters['DataDumpDir'], index,
                                         self.enzoParameters['DataDumpName'], index)
            if os.path.exists(filename):
                pf = load(filename)
                if pf is not None:
                    time_outputs.append({'filename': filename, 'time': pf.current_time})
                    if self.enzoParameters['ComovingCoordinates']:
                        time_outputs[-1]['redshift'] = pf.current_redshift
                del pf
        mylog.info("Located %d time dumps." % len(time_outputs))
        return time_outputs

    def _read_enzo_parameter_file(self):
        "Reads an Enzo parameter file looking for cosmology and output parameters."
        lines = open(self.EnzoParameterFile).readlines()
        for line in lines:
            if line.find("#") >= 0: # Keep the commented lines
                line=line[:line.find("#")]
            line=line.strip()
            if len(line) < 2:
                continue
            try:
                param, vals = map(str,line.split("="))
                param = param.strip()
                vals = vals.strip().split("#", 1)[0].split("//", 1)[0]
            except ValueError:
                continue
            if EnzoParameterDict.has_key(param):
                t = map(EnzoParameterDict[param], vals.split())
                if len(t) == 1:
                    self.enzoParameters[param] = t[0]
                else:
                    self.enzoParameters[param] = t
            elif param.startswith("CosmologyOutputRedshift["):
                index = param[param.find("[")+1:param.find("]")]
                self.redshift_outputs.append({'index':int(index), 'redshift':float(vals)})

        # Add filenames to redshift outputs.
        tempRedshiftList = []
        for output in self.redshift_outputs:
            output["filename"] = "%s/%s%04d/%s%04d" % (self.enzoParameters['GlobalDir'],
                                                       self.enzoParameters['RedshiftDumpDir'],output['index'],
                                                       self.enzoParameters['RedshiftDumpName'],output['index'])
            if not self.get_available_data or os.path.exists(output["filename"]):
                tempRedshiftList.append(output)
        self.redshift_outputs = tempRedshiftList
        del tempRedshiftList

    def _set_parameter_defaults(self):
        "Set some default parameters to avoid problems if they are not in the parameter file."
        self.enzoParameters['GlobalDir'] = "."
        self.enzoParameters['RedshiftDumpName'] = "RedshiftOutput"
        self.enzoParameters['RedshiftDumpDir'] = "RD"
        self.enzoParameters['DataDumpName'] = "data"
        self.enzoParameters['DataDumpDir'] = "DD"
        self.enzoParameters['ComovingCoordinates'] = 0

    def imagine_minimal_splice(self, initial_redshift, final_redshift, decimals=3, filename=None, 
                               redshift_output_string='CosmologyOutputRedshift', start_index=0):
        """
        Create imaginary list of redshift outputs to maximally span a redshift interval.
        :param decimals (int): The decimal place to which the output redshift will be rounded.  
               If the decimal place in question is nonzero, the redshift will be rounded up to 
               ensure continuity of the splice.  Default: 3.
        :param filename (str): If provided, a file will be written with the redshift outputs in 
               the form in which they should be given in the enzo parameter file.  Default: None.
        :param redshift_output_string (str): The parameter accompanying the redshift outputs in the 
               enzo parameter file.  Default: "CosmologyOutputRedshift".
        :param start_index (int): The index of the first redshift output.  Default: 0.
        """

        z = initial_redshift
        outputs = []

        while z > final_redshift:
            rounded = na.round(z, decimals=decimals)
            if rounded - z < 0:
                rounded += na.power(10.0,(-1.0*decimals))
            z = rounded

            deltaz_max = _deltaz_forward(self.cosmology, z, self.enzoParameters['CosmologyComovingBoxSize'])
            outputs.append({'redshift': z, 'deltazMax': deltaz_max})
            z -= deltaz_max

        mylog.info("imagine_maximal_splice: Needed %d data dumps to get from z = %f to %f." %
                   (len(outputs), initial_redshift, final_redshift))

        if filename is not None:
            mylog.info("Writing redshift dump list to %s." % filename)
            f = open(filename,'w')
            for q, output in enumerate(outputs):
                z_string = "%%s[%%d] = %%.%df" % decimals
                f.write(("%s[%d] = %." + str(decimals) + "f\n") % (redshift_output_string, (q+start_index), output['redshift']))
            f.close()

        return outputs

    def create_cosmology_splice(self, minimal=True, deltaz_min=0.0, initial_redshift=None, final_redshift=None):
        """
        Create list of datasets to be used for LightCones or LightRays.
        :param minimal (bool): if True, the minimum number of datasets is used to connect the initial and final 
               redshift.  If false, the list will contain as many entries as possible within the redshift 
               interval.  Default: True.
        :param deltaz_min (float): specifies the minimum delta z between consecutive datasets in the returned 
               list.  Default: 0.0.
        :param initial_redshift (float): the initial (highest) redshift in the cosmology splice list.  If none 
               given, the highest redshift dataset present will be used.  Default: None.
        :param final_redshift (float): the final (lowest) redshift in the cosmology splice list.  If none given, 
               the lowest redshift dataset present will be used.  Default: None.
        """

        if initial_redshift is None: initial_redshift = self.InitialRedshift
        if final_redshift is None: final_redshift = self.FinalRedshift

        # Calculate maximum delta z for each data dump.
        self._calculate_deltaz_max()

        # Calculate minimum delta z for each data dump.
        self._calculate_deltaz_min(deltaz_min=deltaz_min)

        cosmology_splice = []

        # Use minimum number of datasets to go from z_i to z_f.
        if minimal:

            z_Tolerance = 1e-3
            z = initial_redshift

            # fill redshift space with datasets
            while ((z > final_redshift) and 
                   (na.fabs(z - final_redshift) > z_Tolerance)):

                # For first data dump, choose closest to desired redshift.
                if (len(cosmology_splice) == 0):
                    # Sort data outputs by proximity to current redsfhit.
                    self.allOutputs.sort(key=lambda obj:na.fabs(z - obj['redshift']))
                    cosmology_splice.append(self.allOutputs[0])

                # Move forward from last slice in stack until z > z_max.
                else:
                    current_slice = cosmology_splice[-1]
                    while current_slice['next'] is not None and \
                            (z < current_slice['next']['redshift'] or \
                                 na.abs(z - current_slice['next']['redshift']) < z_Tolerance):
                        current_slice = current_slice['next']

                    if current_slice is cosmology_splice[-1]:
                        final_redshift = cosmology_splice[-1]['redshift'] - cosmology_splice[-1]['deltazMax']
                        mylog.error("Cosmology splice incomplete due to insufficient data outputs.")
                        break
                    else:
                        cosmology_splice.append(current_slice)

                z = cosmology_splice[-1]['redshift'] - cosmology_splice[-1]['deltazMax']

        # Make light ray using maximum number of datasets (minimum spacing).
        else:
            # Sort data outputs by proximity to current redsfhit.
            self.allOutputs.sort(key=lambda obj:na.fabs(initial_redshift - obj['redshift']))
            # For first data dump, choose closest to desired redshift.
            cosmology_splice.append(self.allOutputs[0])

            nextOutput = cosmology_splice[-1]['next']
            while (nextOutput is not None):
                if (nextOutput['redshift'] <= final_redshift):
                    break
                if ((cosmology_splice[-1]['redshift'] - nextOutput['redshift']) > cosmology_splice[-1]['deltazMin']):
                    cosmology_splice.append(nextOutput)
                nextOutput = nextOutput['next']
            if (cosmology_splice[-1]['redshift'] - cosmology_splice[-1]['deltazMax']) > final_redshift:
                mylog.error("Cosmology splice incomplete due to insufficient data outputs.")
                final_redshift = cosmology_splice[-1]['redshift'] - cosmology_splice[-1]['deltazMax']

        mylog.info("create_cosmology_splice: Used %d data dumps to get from z = %f to %f." % 
                   (len(cosmology_splice),initial_redshift,final_redshift))

        self.allOutputs.sort(key=lambda obj: obj['time'])
        return cosmology_splice

    def get_data_by_redshift(self, redshifts, tolerance=None):
        """
        : param redshifts: a list of redshifts.
        : tolerance: if not None, do not return a dataset unless the redshift is within the tolerance value.
        Get datasets for a list of redshifts.
        """

        redshifts = ensure_list(redshifts)
        my_datasets = []
        for redshift in redshifts:
            self.allOutputs.sort(key=lambda obj:na.fabs(redshift - obj['redshift']))
            if (tolerance is None or na.abs(redshift - self.allOutputs[0]['redshift']) <= tolerance) \
                    and self.allOutputs[0] not in my_datasets:
                my_datasets.append(self.allOutputs[0])
            else:
                mylog.error("No dataset added for z = %f." % redshift)

        self.allOutputs.sort(key=lambda obj: obj['time'])
        return my_datasets

    def get_data_by_time(self, times, tolerance=None):
        """
        : param redshifts: a list of times.
        : tolerance: if not None, do not return a dataset unless the redshift is within the tolerance value.
        Get datasets for a list of times.
        """

        times = ensure_list(times)
        my_datasets = []
        for my_time in times:
            self.allOutputs.sort(key=lambda obj:na.fabs(my_time - obj['time']))
            if (tolerance is None or na.abs(my_time - self.allOutputs[0]['time']) <= tolerance) \
                    and self.allOutputs[0] not in my_datasets:
                my_datasets.append(self.allOutputs[0])
            else:
                mylog.error("No dataset added for z = %f." % my_time)

        self.allOutputs.sort(key=lambda obj: obj['time'])
        return my_datasets

    def _calculate_deltaz_max(self):
        "Calculate delta z that corresponds to full box length going from z to (z - delta z)."

        d_Tolerance = 1e-4
        max_Iterations = 100

        targetDistance = self.enzoParameters['CosmologyComovingBoxSize']

        for output in self.allOutputs:
            z = output['redshift']

            # Calculate delta z that corresponds to the length of the box at a given redshift.
            # Use Newton's method to calculate solution.
            z1 = z
            z2 = z1 - 0.1 # just an initial guess
            distance1 = 0.0
            iteration = 1

            # Convert comoving radial distance into Mpc / h, since that's how box size is stored.
            distance2 = self.cosmology.ComovingRadialDistance(z2,z) * self.enzoParameters['CosmologyHubbleConstantNow']

            while ((na.fabs(distance2-targetDistance)/distance2) > d_Tolerance):
                m = (distance2 - distance1) / (z2 - z1)
                z1 = z2
                distance1 = distance2
                z2 = ((targetDistance - distance2) / m) + z2
                distance2 = self.cosmology.ComovingRadialDistance(z2,z) * self.enzoParameters['CosmologyHubbleConstantNow']
                iteration += 1
                if (iteration > max_Iterations):
                    mylog.error("calculate_deltaz_max: Warning - max iterations exceeded for z = %f (delta z = %f)." % 
                                (z,na.fabs(z2-z)))
                    break
            output['deltazMax'] = na.fabs(z2-z)

    def _calculate_deltaz_min(self, deltaz_min=0.0):
        "Calculate delta z that corresponds to a single top grid pixel going from z to (z - delta z)."

        d_Tolerance = 1e-4
        max_Iterations = 100

        targetDistance = self.enzoParameters['CosmologyComovingBoxSize'] / self.enzoParameters['TopGridDimensions'][0]

        for output in self.allOutputs:
            z = output['redshift']

            # Calculate delta z that corresponds to the length of a top grid pixel at a given redshift.
            # Use Newton's method to calculate solution.
            z1 = z
            z2 = z1 - 0.01 # just an initial guess
            distance1 = 0.0
            iteration = 1

            # Convert comoving radial distance into Mpc / h, since that's how box size is stored.
            distance2 = self.cosmology.ComovingRadialDistance(z2,z) * self.enzoParameters['CosmologyHubbleConstantNow']

            while ((na.fabs(distance2 - targetDistance) / distance2) > d_Tolerance):
                m = (distance2 - distance1) / (z2 - z1)
                z1 = z2
                distance1 = distance2
                z2 = ((targetDistance - distance2) / m) + z2
                distance2 = self.cosmology.ComovingRadialDistance(z2,z) * self.enzoParameters['CosmologyHubbleConstantNow']
                iteration += 1
                if (iteration > max_Iterations):
                    mylog.error("calculate_deltaz_max: Warning - max iterations exceeded for z = %f (delta z = %f)." % 
                                (z,na.fabs(z2-z)))
                    break
            # Use this calculation or the absolute minimum specified by the user.
            output['deltazMin'] = max(na.fabs(z2-z),deltaz_min)

EnzoParameterDict = {"CosmologyOmegaMatterNow": float,
                     "CosmologyOmegaLambdaNow": float,
                     "CosmologyHubbleConstantNow": float,
                     "CosmologyInitialRedshift": float,
                     "CosmologyFinalRedshift": float,
                     "ComovingCoordinates": int,
                     "InitialTime": float,
                     "StopTime": float,
                     "TopGridDimensions": float,
                     "dtDataDump": float,
                     "RedshiftDumpName": str,
                     "RedshiftDumpDir":  str,
                     "DataDumpName": str,
                     "DataDumpDir": str,
                     "GlobalDir" : str}

def _deltaz_forward(cosmology, z, target_distance):
    "Calculate deltaz corresponding to moving a comoving distance starting from some redshift."

    d_Tolerance = 1e-4
    max_Iterations = 100

    # Calculate delta z that corresponds to the length of the box at a given redshift.
    # Use Newton's method to calculate solution.
    z1 = z
    z2 = z1 - 0.1 # just an initial guess
    distance1 = 0.0
    iteration = 1

    # Convert comoving radial distance into Mpc / h, since that's how box size is stored.
    distance2 = cosmology.ComovingRadialDistance(z2,z) * cosmology.HubbleConstantNow / 100.0

    while ((na.fabs(distance2-target_distance)/distance2) > d_Tolerance):
        m = (distance2 - distance1) / (z2 - z1)
        z1 = z2
        distance1 = distance2
        z2 = ((target_distance - distance2) / m) + z2
        distance2 = cosmology.ComovingRadialDistance(z2,z) * cosmology.HubbleConstantNow / 100.0
        iteration += 1
        if (iteration > max_Iterations):
            mylog.error("deltaz_forward: Warning - max iterations exceeded for z = %f (delta z = %f)." % 
                        (z,na.fabs(z2-z)))
            break
    return na.fabs(z2-z)
