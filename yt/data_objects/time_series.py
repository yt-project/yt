"""
Time series analysis functions.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2010 Matthew Turk.  All Rights Reserved.

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
from .analyzer_objects import create_quantity_proxy
from .derived_quantities import quantity_info

class TimeSeriesData(object):
    def __init__(self, name):
        self.outputs = []

    def __iter__(self):
        # We can make this fancier, but this works
        return self.outputs.__iter__()

class EnzoTimeSeries(TimeSeriesData):
    _enzo_header = "DATASET WRITTEN "
    def __init__(self, name, **kwargs):
        TimeSeriesData.__init__(self, name)
        output_list = kwargs.pop('output_list', None)
        output_log = kwargs.pop('output_log', None)
        if output_list: self._populate_output_list(output_list)
        if output_log: self._populate_output_log(output_log)
        for type_name in data_object_registry:
            setattr(self, type_name, functools.partial(
                TimeSeriesDataObject, self, type_name))

    def _populate_output_list(self, output_list):
        for output in output_list:
            self._insert(EnzoStaticOutput(output))

    def _populate_output_log(self, output_log):
        for line in open(output_log):
            if not line.startswith(self._enzo_header): continue
            fn = line[len(self._enzo_header):].strip()
            self._insert(load(fn))

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
