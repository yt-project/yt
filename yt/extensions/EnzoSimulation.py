import yt.lagos as lagos
from yt.logger import lagosLogger as mylog
import numpy as na

dt_Tolerance = 1e-3

class EnzoSimulation(object):
    """
    Super class for performing the same operation over all data dumps in 
    a simulation from one redshift to another.
    """
    def __init__(self, EnzoParameterFile, initial_time=None, final_time=None, initial_redshift=None, final_redshift=None,
                 links=False, enzo_parameters={}, get_time_outputs=True, get_redshift_outputs=True):
        self.EnzoParameterFile = EnzoParameterFile
        self.enzoParameters = {}
        self.redshiftOutputs = []
        self.timeOutputs = []
        self.allOutputs = []
        self.InitialTime = initial_time
        self.FinalTime = final_time
        self.InitialRedshift = initial_redshift
        self.FinalRedshift = final_redshift
        self.links = links
        self.get_time_outputs = get_time_outputs
        self.get_redshift_outputs = get_redshift_outputs

        # Add any extra parameters to parameter dict.
        EnzoParameterDict.update(enzo_parameters)

        # Set some parameter defaults.
        self._SetParameterDefaults()

        # Read parameters.
        self._ReadEnzoParameterFile()

        # Check for sufficient starting/ending parameters.
        if self.InitialTime is None and self.InitialRedshift is None:
            if self.enzoParameters['ComovingCoordinates'] and \
                    self.enzoParameters.has_key('CosmologyInitialRedshift'):
                self.InitialRedshift = self.enzoParameters['CosmologyInitialRedshift']
            elif self.enzoParameters.has_key('InitialTime'):
                self.InitialTime = self.enzoParameters['InitialTime']
            else:
                mylog.error("Couldn't find parameter for initial time or redshift from parameter file.")
                return None
        if self.FinalTime is None and self.FinalRedshift is None:
            if self.enzoParameters['ComovingCoordinates'] and \
                    self.enzoParameters.has_key('CosmologyFinalRedshift'):
                self.FinalRedshift = self.enzoParameters['CosmologyFinalRedshift']
            elif self.enzoParameters.has_key('StopTime'):
                self.FinalTime = self.enzoParameters['StopTime']
            else:
                mylog.error("Couldn't find parameter for final time or redshift from parameter file.")
                return None

        # Convert initial/final redshifts to times.
        if self.enzoParameters['ComovingCoordinates']:
            # Instantiate a cosmology calculator.
            self.cosmology = lagos.Cosmology(HubbleConstantNow = 
                                             (100.0 * self.enzoParameters['CosmologyHubbleConstantNow']),
                                             OmegaMatterNow = self.enzoParameters['CosmologyOmegaMatterNow'],
                                             OmegaLambdaNow = self.enzoParameters['CosmologyOmegaLambdaNow'])

            # Instantiate EnzoCosmology object for units and time conversions.
            self.enzo_cosmology = lagos.EnzoCosmology(HubbleConstantNow = 
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
                self.enzoParameters.has_key('CosmologyInitialRedshift'):
            self.SimulationInitialTime = self.enzo_cosmology.InitialTime / self.enzo_cosmology.TimeUnits
        elif self.enzoParameters.has_key('InitialTime'):
            self.SimulationInitialTime = self.enzoParameters['InitialTime']
        else:
            self.SimulationInitialTime = 0.0

        # Calculate redshifts for dt data dumps.
        if self.enzoParameters.has_key('dtDataDump'):
            self._CalculateTimeDumps()

        # Calculate times for redshift dumps.
        if self.enzoParameters['ComovingCoordinates']:
            self._CalculateRedshiftDumpTimes()

        # Combine all data dumps.
        self._CombineDataOutputs()

    def _CalculateRedshiftDumpTimes(self):
        "Calculates time from redshift of redshift dumps."

        for output in self.redshiftOutputs:
            output['time'] = self.enzo_cosmology.ComputeTimeFromRedshift(output['redshift']) / \
                self.enzo_cosmology.TimeUnits

    def _CalculateTimeDumps(self):
        "Calculates time dumps and their redshifts if cosmological."

        index = 0
        current_time = self.SimulationInitialTime
        while (current_time <= self.FinalTime) or \
                (abs(self.FinalTime - current_time) / self.FinalTime) < dt_Tolerance:
            filename = "%s/%s%04d/%s%04d" % (self.enzoParameters['GlobalDir'],
                                             self.enzoParameters['DataDumpDir'],index,
                                             self.enzoParameters['DataDumpName'],index)
                                             
            self.timeOutputs.append({'index':index,'filename':filename,'time':current_time})
            if self.enzoParameters['ComovingCoordinates']:
                t = self.enzo_cosmology.InitialTime + (self.enzoParameters['dtDataDump'] * self.enzo_cosmology.TimeUnits * index)
                self.timeOutputs[-1]['redshift'] = self.enzo_cosmology.ComputeRedshiftFromTime(t)

            current_time += self.enzoParameters['dtDataDump']
            index += 1

    def _CombineDataOutputs(self):
        "Combines redshift and time data into one sorted list."

        if not self.get_time_outputs: self.timeOutputs = []
        if not self.get_redshift_outputs: self.redshiftOutputs = []
        self.allOutputs = self.redshiftOutputs + self.timeOutputs
        self.allOutputs.sort(key=lambda obj:obj['time'])

        start_index = None
        end_index = None
        for q in range(len(self.allOutputs)):
            del self.allOutputs[q]['index']

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

            if self.links and start_index is not None:
                if q == start_index:
                    self.allOutputs[q]['previous'] = None
                    self.allOutputs[q]['next'] = self.allOutputs[q+1]
                elif q == end_index:
                    self.allOutputs[q]['previous'] = self.allOutputs[q-1]                
                    self.allOutputs[q]['next'] = None
                elif end_index is None:
                    self.allOutputs[q]['previous'] = self.allOutputs[q-1]
                    self.allOutputs[q]['next'] = self.allOutputs[q+1]

        del self.redshiftOutputs
        del self.timeOutputs

        if end_index is None:
            end_index = len(self.allOutputs)-1

        self.allOutputs = self.allOutputs[start_index:end_index+1]

    def _ReadEnzoParameterFile(self):
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
                self.redshiftOutputs.append({'index':int(index),
                                             'redshift':float(vals)})

        # Add filenames to redshift outputs.
        for output in self.redshiftOutputs:
            output["filename"] = "%s/%s%04d/%s%04d" % (self.enzoParameters['GlobalDir'],
                                                       self.enzoParameters['RedshiftDumpDir'],output['index'],
                                                       self.enzoParameters['RedshiftDumpName'],output['index'])

    def _SetParameterDefaults(self):
        "Set some default parameters to avoid problems if they are not in the parameter file."
        self.enzoParameters['GlobalDir'] = "."
        self.enzoParameters['RedshiftDumpName'] = "RD"
        self.enzoParameters['RedshiftDumpDir'] = "RD"
        self.enzoParameters['DataDumpName'] = "DD"
        self.enzoParameters['DataDumpDir'] = "DD"
        self.enzoParameters['ComovingCoordinates'] = 0

    def _create_cosmology_splice(self, minimal=True, deltaz_min=0.0):
        "Create list of datasets to be used for LightCones or LightRays."

        # Calculate maximum delta z for each data dump.
        self._calculate_deltaz_max()

        # Calculate minimum delta z for each data dump.
        self._calculate_deltaz_min(deltaz_min=deltaz_min)

        cosmology_splice = []

        # Use minimum number of datasets to go from z_i to z_f.
        if minimal:

            z_Tolerance = 1e-4
            z = self.InitialRedshift

            # fill redshift space with datasets
            while ((z > self.FinalRedshift) and 
                   (na.fabs(z - self.FinalRedshift) > z_Tolerance)):
                # Sort data outputs by proximity to current redsfhit.
                self.allOutputs.sort(key=lambda obj:na.fabs(z - obj['redshift']))
                # For first data dump, choose closest to desired redshift.
                if (len(cosmology_splice) == 0):
                    cosmology_splice.append(self.allOutputs[0])
                # Start with data dump closest to desired redshift and move backward 
                # until one is within max delta z of last output in solution list.
                else:
                    output = self.allOutputs[0]
                    while (z > output['redshift']):
                        output = output['previous']
                        if (output is None):
                            mylog.error("CalculateLightRaySolution: search for data output went off the end of the stack.")
                            mylog.error("Could not calculate light ray solution.")
                            return
                        if (output['redshift'] == cosmology_splice[-1]['redshift']):
                            mylog.error("CalculateLightRaySolution: No data dump between z = %f and %f." % \
                                ((cosmology_splice[-1]['redshift'] - cosmology_splice[-1]['deltazMax']),
                                 cosmology_splice[-1]['redshift']))
                            mylog.error("Could not calculate light ray solution.")
                            return
                    cosmology_splice.append(output)
                z = cosmology_splice[-1]['redshift'] - cosmology_splice[-1]['deltazMax']

        # Make light ray using maximum number of datasets (minimum spacing).
        else:
            # Sort data outputs by proximity to current redsfhit.
            self.allOutputs.sort(key=lambda obj:na.fabs(self.InitialRedshift - obj['redshift']))
            # For first data dump, choose closest to desired redshift.
            cosmology_splice.append(self.allOutputs[0])

            nextOutput = cosmology_splice[-1]['next']
            while (nextOutput is not None):
                if (nextOutput['redshift'] <= self.FinalRedshift):
                    break
                if ((cosmology_splice[-1]['redshift'] - nextOutput['redshift']) > cosmology_splice[-1]['deltazMin']):
                    cosmology_splice.append(nextOutput)
                nextOutput = nextOutput['next']

        mylog.info("create_cosmology_splice: Used %d data dumps to get from z = %f to %f." % 
                   (len(cosmology_splice),self.InitialRedshift,self.FinalRedshift))

        return cosmology_splice

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
                    if self.verbose: mylog.error("calculate_deltaz_max: Warning - max iterations exceeded for z = %f (delta z = %f)." % 
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
                    if self.verbose: mylog.error("calculate_deltaz_max: Warning - max iterations exceeded for z = %f (delta z = %f)." % 
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

