import yt.lagos as lagos
import numpy as na
import yt.extensions.HaloProfiler as HP
from EnzoSimulation import *

class SimulationHaloProfiler(EnzoSimulation):
    def __init__(self,EnzoParameterFile,HaloProfilerParameterFile,**kwargs):
        EnzoSimulation.__init__(self,EnzoParameterFile,**kwargs)
        self.HaloProfilerParameterFile = HaloProfilerParameterFile

    def runHaloProfiler(self):
        for output in self.allOutputs:
            halo_profiler = HP.HaloProfiler(output['filename'],self.HaloProfilerParameterFile)
            halo_profiler.makeProfiles()
            #halo_profiler.makeProjections()
            del halo_profiler
