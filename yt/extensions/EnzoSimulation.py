import yt.lagos as lagos
from yt.logger import lagosLogger as mylog
import numpy as na

dt_Tolerance = 1e-3

class EnzoSimulation(object):
    """
    Super class for performing the same operation over all data dumps in 
    a simulation from one redshift to another.
    """
    def __init__(self,EnzoParameterFile,initial_time=None,final_time=None,initial_redshift=None,final_redshift=None,
                 links=False,enzo_parameters={}):
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
            self.cosmology = lagos.EnzoCosmology(HubbleConstantNow = 
                                                 (100.0 * self.enzoParameters['CosmologyHubbleConstantNow']),
                                                 OmegaMatterNow = self.enzoParameters['CosmologyOmegaMatterNow'],
                                                 OmegaLambdaNow = self.enzoParameters['CosmologyOmegaLambdaNow'],
                                                 InitialRedshift = self.enzoParameters['CosmologyInitialRedshift'])
            if self.InitialRedshift is not None:
                self.InitialTime = self.cosmology.ComputeTimeFromRedshift(self.InitialRedshift) / self.cosmology.TimeUnits
            if self.FinalRedshift is not None:
                self.FinalTime = self.cosmology.ComputeTimeFromRedshift(self.FinalRedshift) / self.cosmology.TimeUnits

        # Get initial time of simulation.
        if self.enzoParameters['ComovingCoordinates'] and \
                self.enzoParameters.has_key('CosmologyInitialRedshift'):
            self.SimulationInitialTime = self.cosmology.InitialTime / self.cosmology.TimeUnits
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
            output['time'] = self.cosmology.ComputeTimeFromRedshift(output['redshift']) / self.cosmology.TimeUnits

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
                t = self.cosmology.InitialTime + (self.enzoParameters['dtDataDump'] * self.cosmology.TimeUnits * index)
                self.timeOutputs[-1]['redshift'] = self.cosmology.ComputeRedshiftFromTime(t)

            current_time += self.enzoParameters['dtDataDump']
            index += 1

    def _CombineDataOutputs(self):
        "Combines redshift and time data into one sorted list."
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
        self.enzoParameters['RedshiftDumpName'] = "RedshiftOutput"
        self.enzoParameters['RedshiftDumpDir'] = "RD"
        self.enzoParameters['DataDumpName'] = "DD"
        self.enzoParameters['DataDumpDir'] = "DD"
        self.enzoParameters['ComovingCoordinates'] = 0

    def __iter__(self):
        for output in self.allOutputs:
            yield lagos.EnzoStaticOutput(output['filename'])

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

