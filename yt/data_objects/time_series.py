import functools
import glob
import inspect
import os
import weakref
from functools import wraps
from typing import Optional, Type

import numpy as np
from more_itertools import always_iterable

from yt.config import ytcfg
from yt.data_objects.analyzer_objects import AnalysisTask, create_quantity_proxy
from yt.data_objects.particle_trajectories import ParticleTrajectories
from yt.funcs import is_sequence, mylog
from yt.units.yt_array import YTArray, YTQuantity
from yt.utilities.exceptions import YTException
from yt.utilities.object_registries import (
    analysis_task_registry,
    data_object_registry,
    derived_quantity_registry,
    simulation_time_series_registry,
)
from yt.utilities.parallel_tools.parallel_analysis_interface import (
    communication_system,
    parallel_objects,
    parallel_root_only,
)


class AnalysisTaskProxy:
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

    cls = type(propname, (AnalysisTask,), dict(eval=_eval, _params=tuple()))
    return cls


attrs = (
    "refine_by",
    "dimensionality",
    "current_time",
    "domain_dimensions",
    "domain_left_edge",
    "domain_right_edge",
    "unique_identifier",
    "current_redshift",
    "cosmological_simulation",
    "omega_matter",
    "omega_lambda",
    "omega_radiation",
    "hubble_constant",
)


class TimeSeriesParametersContainer:
    def __init__(self, data_object):
        self.data_object = data_object

    def __getattr__(self, attr):
        if attr in attrs:
            return self.data_object.eval(get_ds_prop(attr)())
        raise AttributeError(attr)


class DatasetSeries:
    r"""The DatasetSeries object is a container of multiple datasets,
    allowing easy iteration and computation on them.

    DatasetSeries objects are designed to provide easy ways to access,
    analyze, parallelize and visualize multiple datasets sequentially.  This is
    primarily expressed through iteration, but can also be constructed via
    analysis tasks (see :ref:`time-series-analysis`).

    Note that contained datasets are lazily loaded and weakly referenced. This means
    that in order to perform follow-up operations on data it's best to define handles on
    these datasets during iteration.

    Parameters
    ----------
    outputs : list of filenames, or pattern
        A list of filenames, for instance ["DD0001/DD0001", "DD0002/DD0002"],
        or a glob pattern (i.e. containing wildcards '[]?!*') such as "DD*/DD*.index".
        In the latter case, results are sorted automatically.
        Filenames and patterns can be of type str, os.Pathlike or bytes.
    parallel : True, False or int
        This parameter governs the behavior when .piter() is called on the
        resultant DatasetSeries object.  If this is set to False, the time
        series will not iterate in parallel when .piter() is called.  If
        this is set to either True, one processor will be allocated for
        each iteration of the loop. If this is set to an integer, the loop
        will be parallelized over this many workgroups. It the integer
        value is less than the total number of available processors,
        more than one processor will be allocated to a given loop iteration,
        causing the functionality within the loop to be run in parallel.
    setup_function : callable, accepts a ds
        This function will be called whenever a dataset is loaded.
    mixed_dataset_types : True or False, default False
        Set to True if the DatasetSeries will load different dataset types, set
        to False if loading dataset of a single type as this will result in a
        considerable speed up from not having to figure out the dataset type.

    Examples
    --------

    >>> ts = DatasetSeries(
    ...     "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0[0-6][0-9]0"
    ... )
    >>> for ds in ts:
    ...     SlicePlot(ds, "x", ("gas", "density")).save()
    ...
    >>> def print_time(ds):
    ...     print(ds.current_time)
    ...
    >>> ts = DatasetSeries(
    ...     "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0[0-6][0-9]0",
    ...     setup_function=print_time,
    ... )
    ...
    >>> for ds in ts:
    ...     SlicePlot(ds, "x", ("gas", "density")).save()

    """
    # this annotation should really be Optional[Type[Dataset]]
    # but we cannot import the yt.data_objects.static_output.Dataset
    # class here without creating a circular import for now
    _dataset_cls: Optional[Type] = None

    def __init_subclass__(cls, *args, **kwargs):
        super().__init_subclass__(*args, **kwargs)
        code_name = cls.__name__[: cls.__name__.find("Simulation")]
        if code_name:
            simulation_time_series_registry[code_name] = cls
            mylog.debug("Registering simulation: %s as %s", code_name, cls)

    def __new__(cls, outputs, *args, **kwargs):
        try:
            outputs = cls._get_filenames_from_glob_pattern(outputs)
        except TypeError:
            pass
        ret = super().__new__(cls)
        ret._pre_outputs = outputs[:]
        ret.kwargs = {}
        return ret

    def __init__(
        self,
        outputs,
        parallel=True,
        setup_function=None,
        mixed_dataset_types=False,
        **kwargs,
    ):
        # This is needed to properly set _pre_outputs for Simulation subclasses.
        self._mixed_dataset_types = mixed_dataset_types
        if is_sequence(outputs) and not isinstance(outputs, str):
            self._pre_outputs = outputs[:]
        self.tasks = AnalysisTaskProxy(self)
        self.params = TimeSeriesParametersContainer(self)
        if setup_function is None:

            def _null(x):
                return None

            setup_function = _null
        self._setup_function = setup_function
        for type_name in data_object_registry:
            setattr(
                self, type_name, functools.partial(DatasetSeriesObject, self, type_name)
            )
        self.parallel = parallel
        self.kwargs = kwargs

    @staticmethod
    def _get_filenames_from_glob_pattern(outputs):
        """
        Helper function to DatasetSeries.__new__
        handle a special case where "outputs" is assumed to be really a pattern string
        """
        pattern = outputs
        epattern = os.path.expanduser(pattern)
        data_dir = ytcfg.get("yt", "test_data_dir")
        # if no match if found from the current work dir,
        # we try to match the pattern from the test data dir
        file_list = glob.glob(epattern) or glob.glob(os.path.join(data_dir, epattern))
        if not file_list:
            raise FileNotFoundError(f"No match found for pattern : {pattern}")
        return sorted(file_list)

    def __getitem__(self, key):
        if isinstance(key, slice):
            if isinstance(key.start, float):
                return self.get_range(key.start, key.stop)
            # This will return a sliced up object!
            return DatasetSeries(
                self._pre_outputs[key], parallel=self.parallel, **self.kwargs
            )
        o = self._pre_outputs[key]
        if isinstance(o, (str, os.PathLike)):
            o = self._load(o, **self.kwargs)
            self._setup_function(o)
        return o

    def __len__(self):
        return len(self._pre_outputs)

    @property
    def outputs(self):
        return self._pre_outputs

    def piter(self, storage=None, dynamic=False):
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
        dynamic : boolean
            This governs whether or not dynamic load balancing will be
            enabled.  This requires one dedicated processor; if this
            is enabled with a set of 128 processors available, only
            127 will be available to iterate over objects as one will
            be load balancing the rest.


        Examples
        --------
        Here is an example of iteration when the results do not need to be
        stored.  One processor will be assigned to each dataset.

        >>> ts = DatasetSeries("DD*/DD*.index")
        >>> for ds in ts.piter():
        ...     SlicePlot(ds, "x", ("gas", "density")).save()
        ...

        This demonstrates how one might store results:

        >>> def print_time(ds):
        ...     print(ds.current_time)
        ...
        >>> ts = DatasetSeries("DD*/DD*.index", setup_function=print_time)
        ...
        >>> my_storage = {}
        >>> for sto, ds in ts.piter(storage=my_storage):
        ...     v, c = ds.find_max(("gas", "density"))
        ...     sto.result = (v, c)
        ...
        >>> for i, (v, c) in sorted(my_storage.items()):
        ...     print("% 4i  %0.3e" % (i, v))
        ...

        This shows how to dispatch 4 processors to each dataset:

        >>> ts = DatasetSeries("DD*/DD*.index", parallel=4)
        >>> for ds in ts.piter():
        ...     ProjectionPlot(ds, "x", ("gas", "density")).save()
        ...

        """
        if not self.parallel:
            njobs = 1
        elif not dynamic:
            if self.parallel:
                njobs = -1
            else:
                njobs = self.parallel
        else:
            my_communicator = communication_system.communicators[-1]
            nsize = my_communicator.size
            if nsize == 1:
                self.parallel = False
                dynamic = False
                njobs = 1
            else:
                njobs = nsize - 1

        for output in parallel_objects(
            self._pre_outputs, njobs=njobs, storage=storage, dynamic=dynamic
        ):
            if storage is not None:
                sto, output = output

            if isinstance(output, str):
                ds = self._load(output, **self.kwargs)
                self._setup_function(ds)
            else:
                ds = output

            if storage is not None:
                next_ret = (sto, ds)
            else:
                next_ret = ds

            yield next_ret

    def eval(self, tasks, obj=None):
        return_values = {}
        for store, ds in self.piter(return_values):
            store.result = []
            for task in always_iterable(tasks):
                try:
                    style = inspect.getargspec(task.eval)[0][1]
                    if style == "ds":
                        arg = ds
                    elif style == "data_object":
                        if obj is None:
                            obj = DatasetSeriesObject(self, "all_data")
                        arg = obj.get(ds)
                    rv = task.eval(arg)
                # We catch and store YT-originating exceptions
                # This fixes the standard problem of having a sphere that's too
                # small.
                except YTException:
                    pass
                store.result.append(rv)
        return [v for k, v in sorted(return_values.items())]

    @classmethod
    def from_output_log(cls, output_log, line_prefix="DATASET WRITTEN", parallel=True):
        filenames = []
        for line in open(output_log):
            if not line.startswith(line_prefix):
                continue
            cut_line = line[len(line_prefix) :].strip()
            fn = cut_line.split()[0]
            filenames.append(fn)
        obj = cls(filenames, parallel=parallel)
        return obj

    def _load(self, output_fn, *, hint: Optional[str] = None, **kwargs):
        from yt.loaders import load

        if self._dataset_cls is not None:
            return self._dataset_cls(output_fn, **kwargs)
        elif self._mixed_dataset_types:
            return load(output_fn, hint=hint, **kwargs)
        ds = load(output_fn, hint=hint, **kwargs)
        self._dataset_cls = ds.__class__
        return ds

    def particle_trajectories(
        self, indices, fields=None, suppress_logging=False, ptype=None
    ):
        r"""Create a collection of particle trajectories in time over a series of
        datasets.

        Parameters
        ----------
        indices : array_like
            An integer array of particle indices whose trajectories we
            want to track. If they are not sorted they will be sorted.
        fields : list of strings, optional
            A set of fields that is retrieved when the trajectory
            collection is instantiated. Default: None (will default
            to the fields 'particle_position_x', 'particle_position_y',
            'particle_position_z')
        suppress_logging : boolean
            Suppress yt's logging when iterating over the simulation time
            series. Default: False
        ptype : str, optional
            Only use this particle type. Default: None, which uses all particle type.

        Examples
        --------
        >>> my_fns = glob.glob("orbit_hdf5_chk_00[0-9][0-9]")
        >>> my_fns.sort()
        >>> fields = [
        ...     ("all", "particle_position_x"),
        ...     ("all", "particle_position_y"),
        ...     ("all", "particle_position_z"),
        ...     ("all", "particle_velocity_x"),
        ...     ("all", "particle_velocity_y"),
        ...     ("all", "particle_velocity_z"),
        ... ]
        >>> ds = load(my_fns[0])
        >>> init_sphere = ds.sphere(ds.domain_center, (0.5, "unitary"))
        >>> indices = init_sphere[("all", "particle_index")].astype("int")
        >>> ts = DatasetSeries(my_fns)
        >>> trajs = ts.particle_trajectories(indices, fields=fields)
        >>> for t in trajs:
        ...     print(
        ...         t[("all", "particle_velocity_x")].max(),
        ...         t[("all", "particle_velocity_x")].min(),
        ...     )

        Notes
        -----
        This function will fail if there are duplicate particle ids or if some of the
        particle disappear.
        """
        return ParticleTrajectories(
            self, indices, fields=fields, suppress_logging=suppress_logging, ptype=ptype
        )


class TimeSeriesQuantitiesContainer:
    def __init__(self, data_object, quantities):
        self.data_object = data_object
        self.quantities = quantities

    def __getitem__(self, key):
        if key not in self.quantities:
            raise KeyError(key)
        q = self.quantities[key]

        def run_quantity_wrapper(quantity, quantity_name):
            @wraps(derived_quantity_registry[quantity_name][1])
            def run_quantity(*args, **kwargs):
                to_run = quantity(*args, **kwargs)
                return self.data_object.eval(to_run)

            return run_quantity

        return run_quantity_wrapper(q, key)


class DatasetSeriesObject:
    def __init__(self, time_series, data_object_name, *args, **kwargs):
        self.time_series = weakref.proxy(time_series)
        self.data_object_name = data_object_name
        self._args = args
        self._kwargs = kwargs
        qs = {
            qn: create_quantity_proxy(qv)
            for qn, qv in derived_quantity_registry.items()
        }
        self.quantities = TimeSeriesQuantitiesContainer(self, qs)

    def eval(self, tasks):
        return self.time_series.eval(tasks, self)

    def get(self, ds):
        # We get the type name, which corresponds to an attribute of the
        # index
        cls = getattr(ds, self.data_object_name)
        return cls(*self._args, **self._kwargs)


class SimulationTimeSeries(DatasetSeries):
    def __init__(self, parameter_filename, find_outputs=False):
        """
        Base class for generating simulation time series types.
        Principally consists of a *parameter_filename*.
        """

        if not os.path.exists(parameter_filename):
            raise FileNotFoundError(parameter_filename)
        self.parameter_filename = parameter_filename
        self.basename = os.path.basename(parameter_filename)
        self.directory = os.path.dirname(parameter_filename)
        self.parameters = {}
        self.key_parameters = []

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

    def _set_parameter_defaults(self):
        pass

    def _parse_parameter_file(self):
        pass

    def _set_units(self):
        pass

    def _calculate_simulation_bounds(self):
        pass

    def _get_all_outputs(**kwargs):
        pass

    def __repr__(self):
        return self.parameter_filename

    _arr = None

    @property
    def arr(self):
        if self._arr is not None:
            return self._arr
        self._arr = functools.partial(YTArray, registry=self.unit_registry)
        return self._arr

    _quan = None

    @property
    def quan(self):
        if self._quan is not None:
            return self._quan
        self._quan = functools.partial(YTQuantity, registry=self.unit_registry)
        return self._quan

    @parallel_root_only
    def print_key_parameters(self):
        """
        Print out some key parameters for the simulation.
        """
        if self.simulation_type == "grid":
            for a in ["domain_dimensions", "domain_left_edge", "domain_right_edge"]:
                self._print_attr(a)
        for a in ["initial_time", "final_time", "cosmological_simulation"]:
            self._print_attr(a)
        if getattr(self, "cosmological_simulation", False):
            for a in [
                "box_size",
                "omega_matter",
                "omega_lambda",
                "omega_radiation",
                "hubble_constant",
                "initial_redshift",
                "final_redshift",
            ]:
                self._print_attr(a)
        for a in self.key_parameters:
            self._print_attr(a)
        mylog.info("Total datasets: %d.", len(self.all_outputs))

    def _print_attr(self, a):
        """
        Print the attribute or warn about it missing.
        """
        if not hasattr(self, a):
            mylog.error("Missing %s in dataset definition!", a)
            return
        v = getattr(self, a)
        mylog.info("Parameters: %-25s = %s", a, v)

    def _get_outputs_by_key(self, key, values, tolerance=None, outputs=None):
        r"""
        Get datasets at or near to given values.

        Parameters
        ----------
        key : str
            The key by which to retrieve outputs, usually 'time' or
            'redshift'.
        values : array_like
            A list of values, given as floats.
        tolerance : float
            If not None, do not return a dataset unless the value is
            within the tolerance value.  If None, simply return the
            nearest dataset.
            Default: None.
        outputs : list
            The list of outputs from which to choose.  If None,
            self.all_outputs is used.
            Default: None.

        Examples
        --------
        >>> datasets = es.get_outputs_by_key("redshift", [0, 1, 2], tolerance=0.1)

        """

        if not isinstance(values, YTArray):
            if isinstance(values, tuple) and len(values) == 2:
                values = self.arr(*values)
            else:
                values = self.arr(values)
        values = values.in_base()

        if outputs is None:
            outputs = self.all_outputs
        my_outputs = []
        if not outputs:
            return my_outputs
        for value in values:
            outputs.sort(key=lambda obj, value=value: np.abs(value - obj[key]))
            if (
                tolerance is None or np.abs(value - outputs[0][key]) <= tolerance
            ) and outputs[0] not in my_outputs:
                my_outputs.append(outputs[0])
            else:
                mylog.error("No dataset added for %s = %f.", key, value)

        outputs.sort(key=lambda obj: obj["time"])
        return my_outputs
