from yt.lagos import *

class TimeSeriesData(object):
    def __init__(self, name):
        self.outputs = []

    def __iter__(self):
        # We can make this fancier, but this works
        return self.outputs.__iter__()

_enzo_header = "DATASET WRITTEN "

class EnzoTimeSeries(TimeSeriesData):
    def __init__(self, name, output_log = "OutputLog"):
        TimeSeriesData.__init__(self, name)
        for line in open(output_log):
            if not line.startswith(_enzo_header): continue
            fn = line[len(_enzo_header):].strip()
            self._insert(EnzoStaticOutput(fn))

    def __getitem__(self, key):
        return self.outputs[key]
        if isinstance(key, types.SliceType):
            if isinstance(key.start, types.FloatType):
                self.get_range(key.start, key.stop)
        
    def _insert(self, pf):
        # We get handed an instantiated parameter file
        # Here we'll figure out a couple things about it, and then stick it
        # inside our list.
        self.outputs.append(pf)
        
    def eval(self, tasks):
        tasks = ensure_list(tasks)
        return_values = []
        for pf in self:
            return_values.append([])
            for task in tasks:
                return_values[-1].append(
                    task.eval(pf.h.all_data()))
        return return_values

class TimeSeriesDataObject(object):
    def __init__(self, data_object_cls, *args, **kwargs):
        pass
