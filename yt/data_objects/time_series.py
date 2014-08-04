"""
Time series analysis functions.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import inspect, functools, weakref, glob, types, os

from yt.funcs import *
from yt.extern.six import add_metaclass
from yt.convenience import load
from yt.config import ytcfg
from .data_containers import data_object_registry
from .analyzer_objects import create_quantity_proxy, \
    analysis_task_registry, AnalysisTask
from yt.units.yt_array import YTArray, YTQuantity
from yt.utilities.exceptions import YTException
from yt.utilities.parallel_tools.parallel_analysis_interface \
    import parallel_objects, parallel_root_only
from yt.utilities.parameter_file_storage import \
    simulation_time_series_registry
     
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

def get_ds_prop(propname):
    def _eval(params, ds):
        return getattr(ds, propname)
    cls = type(propname, (AnalysisTask,),
                dict(eval = _eval, _params = tuple()))
    return cls

def get_filenames_from_glob_pattern(filenames):
    file_list = glob.glob(filenames)
    if len(file_list) == 0:
        data_dir = ytcfg.get("yt", "test_data_dir")
        pattern = os.path.join(data_dir, filenames)
        td_filenames = glob.glob(pattern)
        if len(td_filenames) > 0:
            file_list = td_filenames
        else:
            raise YTOutputNotIdentified(filenames, {})
    return sorted(file_list)

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
            return self.data_object.eval(get_ds_prop(attr)())
        raise AttributeError(attr)

class DatasetSeries(object):
    r"""The DatasetSeries object is a container of multiple datasets,
    allowing easy iteration and computation on them.

    DatasetSeries objects are designed to provide easy ways to access,
    analyze, parallelize and visualize multiple datasets sequentially.  This is
    primarily expressed through iteration, but can also be constructed via
    analysis tasks (see :ref:`time-series-analysis`).

    Parameters
    ----------
    filenames : list or pattern
        This can either be a list of filenames (such as ["DD0001/DD0001",
        "DD0002/DD0002"]) or a pattern to match, such as
        "DD*/DD*.index").  If it's the former, they will be loaded in
        order.  The latter will be identified with the glob module and then
        sorted.
    parallel : True, False or int
        This parameter governs the behavior when .piter() is called on the
        resultant DatasetSeries object.  If this is set to False, the time
        series will not iterate in parallel when .piter() is called.  If
        this is set to either True or an integer, it will be iterated with
        1 or that integer number of processors assigned to each parameter
        file provided to the loop.
    setup_function : callable, accepts a ds
        This function will be called whenever a dataset is loaded.

    Examples
    --------

    >>> ts = DatasetSeries(
            "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0[0-6][0-9]0")
    >>> for ds in ts:
    ...     SlicePlot(ds, "x", "Density").save()
    ...
    >>> def print_time(ds):
    ...     print ds.current_time
    ...
    >>> ts = DatasetSeries(
    ...     "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0[0-6][0-9]0",
    ...      setup_function = print_time)
    ...
    >>> for ds in ts:
    ...     SlicePlot(ds, "x", "Density").save()

    """
    def __new__(cls, outputs, *args, **kwargs):
        if isinstance(outputs, basestring):
            outputs = get_filenames_from_glob_pattern(outputs)
        ret = super(DatasetSeries, cls).__new__(cls, *args, **kwargs)
        try:
            ret._pre_outputs = outputs[:]
        except TypeError:
            raise YTOutputNotIdentified(outputs, {})
        return ret

    def __init__(self, outputs, parallel = True, setup_function = None,
                 **kwargs):
        # This is needed to properly set _pre_outputs for Simulation subclasses.
        if iterable(outputs) and not isinstance(outputs, basestring):
            self._pre_outputs = outputs[:]
        self.tasks = AnalysisTaskProxy(self)
        self.params = TimeSeriesParametersContainer(self)
        if setup_function is None:
            setup_function = lambda a: None
        self._setup_function = setup_function
        for type_name in data_object_registry:
            setattr(self, type_name, functools.partial(
                DatasetSeriesObject, self, type_name))
        self.parallel = parallel
        self.kwargs = kwargs

    def __iter__(self):
        # We can make this fancier, but this works
        for o in self._pre_outputs:
            if isinstance(o, types.StringTypes):
                ds = load(o, **self.kwargs)
                self._setup_function(ds)
                yield ds
            else:
                yield o

    def __getitem__(self, key):
        if isinstance(key, types.SliceType):
            if isinstance(key.start, types.FloatType):
                return self.get_range(key.start, key.stop)
            # This will return a sliced up object!
            return DatasetSeries(self._pre_outputs[key], self.parallel)
        o = self._pre_outputs[key]
        if isinstance(o, types.StringTypes):
            o = load(o, **self.kwargs)
            self._setup_function(o)
        return o

    def __len__(self):
        return len(self._pre_outputs)

    @property
    def outputs(self):
        return self._pre_outputs

    def piter(self, storage = None):
        r"""Iterate over time series components in parallel.

        This allows you to iterate over a time series while dispatching
        individual components of that time series to different processors or
        processor groups.  If the parallelism strategy was set to be
        multi-processor (by "parallel = N" where N is an integer when the
        DatasetSeries was created) this will issue each dataset to an
        N-processor group.  For instance, this would allow you to start a 1024
        processor job, loading up 100 datasets in a time series and creating 8
        processor groups of 128 processors each, each of which would be
        assigned a different dataset.  This could be accomplished as shown in
        the examples below.  The *storage* option is as seen in
        :func:`~yt.utilities.parallel_tools.parallel_analysis_interface.parallel_objects`
        which is a mechanism for storing results of analysis on an individual
        dataset and then combining the results at the end, so that the entire
        set of processors have access to those results.

        Note that supplying a *store* changes the iteration mechanism; see
        below.

        Parameters
        ----------
        storage : dict
            This is a dictionary, which will be filled with results during the
            course of the iteration.  The keys will be the dataset
            indices and the values will be whatever is assigned to the *result*
            attribute on the storage during iteration.

        Examples
        --------
        Here is an example of iteration when the results do not need to be
        stored.  One processor will be assigned to each dataset.

        >>> ts = DatasetSeries("DD*/DD*.index")
        >>> for ds in ts.piter():
        ...    SlicePlot(ds, "x", "Density").save()
        ...
        
        This demonstrates how one might store results:

        >>> def print_time(ds):
        ...     print ds.current_time
        ...
        >>> ts = DatasetSeries("DD*/DD*.index",
        ...             setup_function = print_time )
        ...
        >>> my_storage = {}
        >>> for sto, ds in ts.piter(storage=my_storage):
        ...     v, c = ds.find_max("density")
        ...     sto.result = (v, c)
        ...
        >>> for i, (v, c) in sorted(my_storage.items()):
        ...     print "% 4i  %0.3e" % (i, v)
        ...

        This shows how to dispatch 4 processors to each dataset:

        >>> ts = DatasetSeries("DD*/DD*.index",
        ...                     parallel = 4)
        >>> for ds in ts.piter():
        ...     ProjectionPlot(ds, "x", "Density").save()
        ...

        """
        dynamic = False
        if self.parallel == False:
            njobs = 1
        else:
            if self.parallel == True: njobs = -1
            else: njobs = self.parallel
        return parallel_objects(self, njobs=njobs, storage=storage,
                                dynamic=dynamic)

    def eval(self, tasks, obj=None):
        tasks = ensure_list(tasks)
        return_values = {}
        for store, ds in self.piter(return_values):
            store.result = []
            for task in tasks:
                try:
                    style = inspect.getargspec(task.eval)[0][1]
                    if style == 'ds':
                        arg = ds
                    elif style == 'data_object':
                        if obj == None:
                            obj = DatasetSeriesObject(self, "all_data")
                        arg = obj.get(ds)
                    rv = task.eval(arg)
                # We catch and store YT-originating exceptions
                # This fixes the standard problem of having a sphere that's too
                # small.
                except YTException as rv:
                    pass
                store.result.append(rv)
        return [v for k, v in sorted(return_values.items())]

    @classmethod
    def from_filenames(cls, filenames, parallel = True, setup_function = None,
                       **kwargs):
        r"""Create a time series from either a filename pattern or a list of
        filenames.

        This method provides an easy way to create a
        :class:`~yt.data_objects.time_series.DatasetSeries`, given a set of
        filenames or a pattern that matches them.  Additionally, it can set the
        parallelism strategy.

        Parameters
        ----------
        filenames : list or pattern
            This can either be a list of filenames (such as ["DD0001/DD0001",
            "DD0002/DD0002"]) or a pattern to match, such as
            "DD*/DD*.index").  If it's the former, they will be loaded in
            order.  The latter will be identified with the glob module and then
            sorted.
        parallel : True, False or int
            This parameter governs the behavior when .piter() is called on the
            resultant DatasetSeries object.  If this is set to False, the time
            series will not iterate in parallel when .piter() is called.  If
            this is set to either True or an integer, it will be iterated with
            1 or that integer number of processors assigned to each parameter
            file provided to the loop.
        setup_function : callable, accepts a ds
            This function will be called whenever a dataset is loaded.

        Examples
        --------

        >>> def print_time(ds):
        ...     print ds.current_time
        ...
        >>> ts = DatasetSeries.from_filenames(
        ...     "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0[0-6][0-9]0",
        ...      setup_function = print_time)
        ...
        >>> for ds in ts:
        ...     SlicePlot(ds, "x", "Density").save()

        """
        
        if isinstance(filenames, types.StringTypes):
            filenames = get_filenames_from_glob_pattern(filenames)
        obj = cls(filenames[:], parallel = parallel,
                  setup_function = setup_function, **kwargs)
        return obj

    @classmethod
    def from_output_log(cls, output_log,
                        line_prefix = "DATASET WRITTEN",
                        parallel = True):
        filenames = []
        for line in open(output_log):
            if not line.startswith(line_prefix): continue
            cut_line = line[len(line_prefix):].strip()
            fn = cut_line.split()[0]
            filenames.append(fn)
        obj = cls(filenames, parallel = parallel)
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

class DatasetSeriesObject(object):
    def __init__(self, time_series, data_object_name, *args, **kwargs):
        self.time_series = weakref.proxy(time_series)
        self.data_object_name = data_object_name
        self._args = args
        self._kwargs = kwargs
        qs = dict([(qn, create_quantity_proxy(qv)) for qn, qv in quantity_info.items()])
        self.quantities = TimeSeriesQuantitiesContainer(self, qs)

    def eval(self, tasks):
        return self.time_series.eval(tasks, self)

    def get(self, ds):
        # We get the type name, which corresponds to an attribute of the
        # index
        cls = getattr(ds, self.data_object_name)
        return cls(*self._args, **self._kwargs)

class RegisteredSimulationTimeSeries(type):
    def __init__(cls, name, b, d):
        type.__init__(cls, name, b, d)
        code_name = name[:name.find('Simulation')]
        if code_name:
            simulation_time_series_registry[code_name] = cls
            mylog.debug("Registering simulation: %s as %s", code_name, cls)

@add_metaclass(RegisteredSimulationTimeSeries)
class SimulationTimeSeries(DatasetSeries):
    def __init__(self, parameter_filename, find_outputs=False):
        """
        Base class for generating simulation time series types.
        Principally consists of a *parameter_filename*.
        """

        if not os.path.exists(parameter_filename):
            raise IOError(parameter_filename)
        self.parameter_filename = parameter_filename
        self.basename = os.path.basename(parameter_filename)
        self.directory = os.path.dirname(parameter_filename)
        self.parameters = {}

        # Set some parameter defaults.
        self._set_parameter_defaults()
        # Read the simulation dataset.
        self._parse_parameter_file()
        # Set units
        self._set_units()
        # Figure out the starting and stopping times and redshift.
        self._calculate_simulation_bounds()
        # Get all possible datasets.
        self._get_all_outputs(find_outputs=find_outputs)
        
        self.print_key_parameters()

    def __repr__(self):
        return self.parameter_filename

    _arr = None
    @property
    def arr(self):
        if self._arr is not None:
            return self._arr
        self._arr = functools.partial(YTArray, registry = self.unit_registry)
        return self._arr
    
    _quan = None
    @property
    def quan(self):
        if self._quan is not None:
            return self._quan
        self._quan = functools.partial(YTQuantity,
                registry = self.unit_registry)
        return self._quan
    
    @parallel_root_only
    def print_key_parameters(self):
        """
        Print out some key parameters for the simulation.
        """
        for a in ["domain_dimensions", "domain_left_edge",
                  "domain_right_edge", "initial_time", "final_time",
                  "stop_cycle", "cosmological_simulation"]:
            if not hasattr(self, a):
                mylog.error("Missing %s in dataset definition!", a)
                continue
            v = getattr(self, a)
            mylog.info("Parameters: %-25s = %s", a, v)
        if hasattr(self, "cosmological_simulation") and \
           getattr(self, "cosmological_simulation"):
            for a in ["box_size", "omega_lambda",
                      "omega_matter", "hubble_constant",
                      "initial_redshift", "final_redshift"]:
                if not hasattr(self, a):
                    mylog.error("Missing %s in dataset definition!", a)
                    continue
                v = getattr(self, a)
                mylog.info("Parameters: %-25s = %s", a, v)
        mylog.info("Total datasets: %d." % len(self.all_outputs))

