import glob
import os
from yt.extern.six import add_metaclass
from yt.convenience import \
    load
from yt.funcs import \
    only_on_root
from yt.utilities.exceptions import \
    YTOutputNotIdentified
from yt.utilities.logger import ytLogger as \
    mylog
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_objects
from yt.data_objects.time_series import \
    DatasetSeries, \
    RegisteredSimulationTimeSeries
from yt.frontends.exodus_ii.api import ExodusIIDataset


@add_metaclass(RegisteredSimulationTimeSeries)
class ExodusIISimulation(DatasetSeries):
    r"""
    Initialize an ExodusII Simulation object.

    Upon creation, the input directoy is searched for valid ExodusIIDatasets.
    The get_time_series can be used to generate a DatasetSeries object.

    simulation_directory : str
        The directory that contain the simulation data.
    
    Examples
    --------
    >>> import yt
    >>> sim = yt.simulation("demo_second", "ExodusII")
    >>> sim.get_time_series()
    >>> for ds in sim:
    ...     print ds.current_time

    """
    
    def __init__(self, simulation_directory, find_outputs=False):
        self.simulation_directory = simulation_directory
        fn_pattern = "%s/*" % self.simulation_directory
        potential_outputs = glob.glob(fn_pattern)
        self.all_outputs = self._check_for_outputs(potential_outputs)
        self.all_outputs.sort(key=lambda obj: obj["filename"])

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
            num_steps = output['num_steps']
            fn = output['filename']
            for step in range(num_steps):
                ds = ExodusIIDataset(fn, step=step)
                ds_list.append(ds)
        super(ExodusIISimulation, self).__init__(ds_list, 
                                                 parallel=parallel, 
                                                 setup_function=setup_function)
        
    def _check_for_outputs(self, potential_outputs):
        r"""
        Check a list of files to see if they are valid datasets.
        """

        only_on_root(mylog.info, "Checking %d potential outputs.",
                     len(potential_outputs))

        my_outputs = {}
        for my_storage, output in parallel_objects(potential_outputs,
                                                   storage=my_outputs):
            if os.path.exists(output):
                try:
                    ds = load(output)
                    if ds is not None:
                        num_steps = ds.num_steps
                        my_storage.result = {"filename": output,
                                             "num_steps": num_steps}
                except YTOutputNotIdentified:
                    mylog.error("Failed to load %s", output)
        my_outputs = [my_output for my_output in my_outputs.values() \
                      if my_output is not None]
        return my_outputs
