from yt.lagos import *
import inspect

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
        
    def eval(self, tasks, obj=None):
        if obj == None: obj = TimeSeriesDataObject(self, "all_data")
        tasks = ensure_list(tasks)
        return_values = []
        for pf in self:
            return_values.append([])
            for task in tasks:
                style = inspect.getargspec(task.eval)[0][1]
                if style == 'pf': arg = pf
                elif style == 'data_object': arg = obj.get(pf)
                return_values[-1].append(task.eval(arg))
        return return_values

class TimeSeriesDataObject(object):
    def __init__(self, time_series, data_object_name, *args, **kwargs):
        self.time_series = weakref.proxy(time_series)
        self.data_object_name = data_object_name
        self._args = args
        self._kwargs = kwargs

    def eval(self, tasks):
        self.time_series.eval(tasks, self)

    def get(self, pf):
        # We get the type name, which corresponds to an attribute of the
        # hierarchy
        cls = getattr(pf.h, self.data_object_name)
        return cls(*self._args, **self._kwargs)
