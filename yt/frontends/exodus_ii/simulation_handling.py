import glob

from yt.data_objects.time_series import DatasetSeries
from yt.funcs import only_on_root
from yt.loaders import load
from yt.utilities.exceptions import YTUnidentifiedDataType
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.parallel_tools.parallel_analysis_interface import parallel_objects


class ExodusIISimulation(DatasetSeries):
    r"""
    Initialize an ExodusII Simulation object.

    Upon creation, the input directory is searched for valid ExodusIIDatasets.
    The get_time_series can be used to generate a DatasetSeries object.

    simulation_directory : str
        The directory that contain the simulation data.

    Examples
    --------
    >>> import yt
    >>> sim = yt.load_simulation("demo_second", "ExodusII")
    >>> sim.get_time_series()
    >>> for ds in sim:
    ...     print(ds.current_time)

    """

    def __init__(self, simulation_directory, find_outputs=False):
        self.simulation_directory = simulation_directory
        fn_pattern = f"{self.simulation_directory}/*"
        potential_outputs = glob.glob(fn_pattern)
        self.all_outputs = self._check_for_outputs(potential_outputs)
        self.all_outputs.sort(key=lambda obj: obj["filename"])

    def __iter__(self):
        for o in self._pre_outputs:
            fn, step = o
            ds = load(fn, step=step)
            self._setup_function(ds)
            yield ds

    def __getitem__(self, key):
        if isinstance(key, slice):
            if isinstance(key.start, float):
                return self.get_range(key.start, key.stop)
            # This will return a sliced up object!
            return DatasetSeries(self._pre_outputs[key], self.parallel)
        o = self._pre_outputs[key]
        fn, step = o
        o = load(fn, step=step)
        self._setup_function(o)
        return o

    def get_time_series(self, parallel=False, setup_function=None):
        r"""
        Instantiate a DatasetSeries object for a set of outputs.

        If no additional keywords given, a DatasetSeries object will be
        created with all potential datasets created by the simulation.

        Fine-level filtering is currently not implemented.

        """

        all_outputs = self.all_outputs
        ds_list = []
        for output in all_outputs:
            num_steps = output["num_steps"]
            fn = output["filename"]
            for step in range(num_steps):
                ds_list.append((fn, step))
        super().__init__(ds_list, parallel=parallel, setup_function=setup_function)

    def _check_for_outputs(self, potential_outputs):
        r"""
        Check a list of files to see if they are valid datasets.
        """

        only_on_root(
            mylog.info, "Checking %d potential outputs.", len(potential_outputs)
        )

        my_outputs = {}
        for my_storage, output in parallel_objects(
            potential_outputs, storage=my_outputs
        ):
            try:
                ds = load(output)
            except (FileNotFoundError, YTUnidentifiedDataType):
                mylog.error("Failed to load %s", output)
                continue
            my_storage.result = {"filename": output, "num_steps": ds.num_steps}

        my_outputs = [
            my_output for my_output in my_outputs.values() if my_output is not None
        ]
        return my_outputs
