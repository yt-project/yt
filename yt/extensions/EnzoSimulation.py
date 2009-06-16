import yt.lagos as lagos
import numpy as na

class EnzoSimulation(object):
    """
    Super class for performing the same operation over all data dumps in 
    a simulation from one redshift to another.
    """
    def __init__(self,EnzoParameterFile,initial_time=None,final_time=None,initial_redshift=None,final_redshift=None,
                 links=False):
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

        ec = lagos.EnzoCosmology(HubbleConstantNow = \
                                     (100.0 * self.enzoParameters['CosmologyHubbleConstantNow']),
                                 OmegaMatterNow = self.enzoParameters['CosmologyOmegaMatterNow'],
                                 OmegaLambdaNow = self.enzoParameters['CosmologyOmegaLambdaNow'],
                                 InitialRedshift = self.enzoParameters['CosmologyInitialRedshift'])

        for output in self.redshiftOutputs:
            output['time'] = ec.ComputeTimeFromRedshift(output['redshift']) / ec.TimeUnits

    def _CalculateTimeDumps(self):
        "Calculates time dumps and their redshifts if cosmological."

        initial_time = self.InitialTime
        final_time = self.FinalTime

        if self.enzoParameters['ComovingCoordinates']:
            ec = lagos.EnzoCosmology(HubbleConstantNow = \
                                         (100.0 * self.enzoParameters['CosmologyHubbleConstantNow']),
                                     OmegaMatterNow = self.enzoParameters['CosmologyOmegaMatterNow'],
                                     OmegaLambdaNow = self.enzoParameters['CosmologyOmegaLambdaNow'],
                                     InitialRedshift = self.enzoParameters['CosmologyInitialRedshift'])
            if initial_time is None:
                initial_time = ec.ComputeTimeFromRedshift(self.InitialRedshift) / ec.TimeUnits
            if final_time is None:
                final_time = ec.ComputeTimeFromRedshift(self.FinalRedshift) / ec.TimeUnits

        index = 0
        t_Tolerance = 1e-3
        current_time = initial_time
        while ((current_time <= final_time) and 
               (na.fabs(current_time - final_time) > t_Tolerance)):
            filename = "%s/%s%04d/%s%04d" % (self.enzoParameters['GlobalDir'],
                                             self.enzoParameters['DataDumpDir'],index,
                                             self.enzoParameters['DataDumpName'],index)
                                             
            self.timeOutputs.append({'index':index,'filename':filename,'time':current_time})
            if self.enzoParameters['ComovingCoordinates']:
                t = ec.InitialTime + (self.enzoParameters['dtDataDump'] * ec.TimeUnits * index)
                self.timeOutputs[-1]['redshift'] = ec.ComputeRedshiftFromTime(t)

            current_time += self.enzoParameters['dtDataDump']
            index += 1

    def _CombineDataOutputs(self):
        "Combines redshift and time data into one sorted list."
        self.allOutputs = self.redshiftOutputs + self.timeOutputs
        self.allOutputs.sort(reverse=True,key=lambda obj:obj['redshift'])
        start_index = None
        end_index = None
        for q in range(len(self.allOutputs)):
            del self.allOutputs[q]['index']
            # set beginning and end indices from redshift limits
            if start_index is None:
                if self.InitialRedshift is not None and \
                        self.allOutputs[q]['redshift'] <= self.InitialRedshift:
                    start_index = q
                else:
                    self.allOutputs[q]['time'] >= self.InitialTime
            if end_index is None:
                if self.FinalRedshift is not None:
                    if self.allOutputs[q]['redshift'] < self.FinalRedshift:
                        end_index = q - 1
                    elif self.allOutputs[q]['redshift'] == self.FinalRedshift:
                        end_index = q
                else:
                    if self.allOutputs[q]['time'] > self.FinalTime:
                        end_index = q - 1
                    elif self.allOutputs[q]['time'] == self.FinalTime:
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
                vals = vals.strip()
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
        self.enzoParameters['GlobalDir'] = ""
        self.enzoParameters['RedshiftDumpName'] = "RD"
        self.enzoParameters['RedshiftDumpDir'] = "RD"
        self.enzoParameters['DataDumpName'] = "DD"
        self.enzoParameters['DataDumpDir'] = "DD"
        self.enzoParameters['ComovingCoordinates'] = 0

EnzoParameterDict = {"CosmologyCurrentRedshift": float,
                     "CosmologyComovingBoxSize": float,
                     "CosmologyOmegaMatterNow": float,
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

