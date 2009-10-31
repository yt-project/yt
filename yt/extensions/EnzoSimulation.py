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
                 links=False, enzo_parameters=None, get_time_outputs=True, get_redshift_outputs=True):
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
        if enzo_parameters is None: enzo_parameters = {}
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

        # Combine all data dumps.
        self._CombineDataOutputs()

    def _CalculateRedshiftDumpTimes(self):
        "Calculates time from redshift of redshift dumps."

        for output in self.redshiftOutputs:
            output['time'] = self.enzo_cosmology.ComputeTimeFromRedshift(output['redshift']) / \
                self.enzo_cosmology.TimeUnits

    def _CalculateTimeDumps(self):
        "Calculates time dumps and their redshifts if cosmological."

        if self.enzoParameters['dtDataDump'] <= 0.0: return

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

        # Calculate redshifts for dt data dumps.
        if self.enzoParameters.has_key('dtDataDump') and self.get_time_outputs:
            self._CalculateTimeDumps()

        # Calculate times for redshift dumps.
        if self.enzoParameters['ComovingCoordinates'] and self.get_redshift_outputs:
            self._CalculateRedshiftDumpTimes()

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
                if q == len(self.allOutputs) - 1:
                    end_index = q

            if self.links and start_index is not None:
                if q == start_index:
                    self.allOutputs[q]['previous'] = None
                else:
                    self.allOutputs[q]['previous'] = self.allOutputs[q-1]

                if q == end_index:
                    self.allOutputs[q]['next'] = None
                else:
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

    def imagine_minimal_splice(self, initial_redshift, final_redshift, decimals=3, filename=None, 
                               redshift_output_string='CosmologyOutputRedshift', start_index=0):
        "Create imaginary list of redshift outputs to maximally span a redshift interval."

        z = initial_redshift
        outputs = []

        while z > final_redshift:
            rounded = na.round(z, decimals=decimals)
            if rounded - z < 0:
                rounded += na.power(10.0,(-1.0*decimals))
            z = rounded
            deltaz_max = deltaz_forward(self.cosmology, z, self.enzoParameters['CosmologyComovingBoxSize'])
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

    def _create_cosmology_splice(self, minimal=True, deltaz_min=0.0, initial_redshift=None, final_redshift=None):
        "Create list of datasets to be used for LightCones or LightRays."

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

def deltaz_forward(cosmology, z, target_distance):
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
