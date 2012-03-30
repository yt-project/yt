"""
Time series analysis functions.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2010-2011 Matthew Turk.  All Rights Reserved.

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

import inspect, functools, weakref

from yt.funcs import *
from yt.convenience import load
from .data_containers import data_object_registry
from .analyzer_objects import create_quantity_proxy, \
    analysis_task_registry, AnalysisTask
from .derived_quantities import quantity_info
from yt.utilities.exceptions import YTException

class AnalysisTaskProxy(object):
    def __init__(self, time_series):
        self.time_series = time_series

    def __getitem__(self, key):
        task_cls = analysis_task_registry[key]
        @wraps(task_cls.__init__)
        def func(*args, **kwargs):
            task = task_cls(*args, **kwargs)
            return self.time_series.eval(task)
        return func

    def keys(self):
        return analysis_task_registry.keys()

    def __contains__(self, key):
        return key in analysis_task_registry

def get_pf_prop(propname):
    def _eval(params, pf):
        return getattr(pf, propname)
    cls = type(propname, (AnalysisTask,),
                dict(eval = _eval, _params = tuple()))
    return cls

attrs = ("refine_by", "dimensionality", "current_time",
         "domain_dimensions", "domain_left_edge",
         "domain_right_edge", "unique_identifier",
         "current_redshift", "cosmological_simulation",
         "omega_matter", "omega_lambda", "hubble_constant")

class TimeSeriesParametersContainer(object):
    def __init__(self, data_object):
        self.data_object = data_object

    def __getattr__(self, attr):
        if attr in attrs:
            return self.data_object.eval(get_pf_prop(attr)())
        raise AttributeError(attr)

class TimeSeriesData(object):
    def __init__(self, outputs = None):
        if outputs is None: outputs = []
        self.outputs = outputs
        self.tasks = AnalysisTaskProxy(self)
        self.params = TimeSeriesParametersContainer(self)
        for type_name in data_object_registry:
            setattr(self, type_name, functools.partial(
                TimeSeriesDataObject, self, type_name))

    def __iter__(self):
        # We can make this fancier, but this works
        return self.outputs.__iter__()

    def __getitem__(self, key):
        if isinstance(key, types.SliceType):
            if isinstance(key.start, types.FloatType):
                return self.get_range(key.start, key.stop)
        return self.outputs[key]
        
    def _insert(self, pf):
        # We get handed an instantiated parameter file
        # Here we'll figure out a couple things about it, and then stick it
        # inside our list.
        self.outputs.append(pf)
        
    def eval(self, tasks, obj=None):
        tasks = ensure_list(tasks)
        return_values = []
        for pf in self:
            return_values.append([])
            for task in tasks:
                try:
                    style = inspect.getargspec(task.eval)[0][1]
                    if style == 'pf':
                        arg = pf
                    elif style == 'data_object':
                        if obj == None:
                            obj = TimeSeriesDataObject(self, "all_data")
                        arg = obj.get(pf)
                    rv = task.eval(arg)
                # We catch and store YT-originating exceptions
                # This fixes the standard problem of having a sphere that's too
                # small.
                except YTException as rv:
                    pass
                return_values[-1].append(rv)
        return return_values

    @classmethod
    def from_filenames(cls, filename_list):
        outputs = []
        for fn in filename_list:
            outputs.append(load(fn))
        obj = cls(outputs)
        return obj

    @classmethod
    def from_output_log(cls, output_log,
                        line_prefix = "DATASET WRITTEN"):
        outputs = []
        for line in open(output_log):
            if not line.startswith(line_prefix): continue
            cut_line = line[len(line_prefix):].strip()
            fn = cut_line.split()[0]
            outputs.append(load(fn))
        obj = cls(outputs)
        return obj

class TimeSeriesQuantitiesContainer(object):
    def __init__(self, data_object, quantities):
        self.data_object = data_object
        self.quantities = quantities

    def __getitem__(self, key):
        if key not in self.quantities: raise KeyError(key)
        q = self.quantities[key]
        def run_quantity_wrapper(quantity, quantity_name):
            @wraps(quantity_info[quantity_name][1])
            def run_quantity(*args, **kwargs):
                to_run = quantity(*args, **kwargs)
                return self.data_object.eval(to_run)
            return run_quantity
        return run_quantity_wrapper(q, key)

class TimeSeriesDataObject(object):
    def __init__(self, time_series, data_object_name, *args, **kwargs):
        self.time_series = weakref.proxy(time_series)
        self.data_object_name = data_object_name
        self._args = args
        self._kwargs = kwargs
        qs = dict([(qn, create_quantity_proxy(qv)) for qn, qv in quantity_info.items()])
        self.quantities = TimeSeriesQuantitiesContainer(self, qs)

    def eval(self, tasks):
        return self.time_series.eval(tasks, self)

    def get(self, pf):
        # We get the type name, which corresponds to an attribute of the
        # hierarchy
        cls = getattr(pf.h, self.data_object_name)
        return cls(*self._args, **self._kwargs)
