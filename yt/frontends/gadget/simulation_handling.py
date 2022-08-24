import glob
import os

import numpy as np
from unyt import dimensions, unyt_array
from unyt.unit_registry import UnitRegistry

from yt.data_objects.time_series import DatasetSeries, SimulationTimeSeries
from yt.funcs import only_on_root
from yt.loaders import load
from yt.utilities.cosmology import Cosmology
from yt.utilities.exceptions import (
    InvalidSimulationTimeSeries,
    MissingParameter,
    NoStoppingCondition,
    YTUnidentifiedDataType,
)
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.parallel_tools.parallel_analysis_interface import parallel_objects


class GadgetSimulation(SimulationTimeSeries):
    r"""
    Initialize an Gadget Simulation object.

    Upon creation, the parameter file is parsed and the time and redshift
    are calculated and stored in all_outputs.  A time units dictionary is
    instantiated to allow for time outputs to be requested with physical
    time units.  The get_time_series can be used to generate a
    DatasetSeries object.

    parameter_filename : str
        The simulation parameter file.
    find_outputs : bool
        If True, the OutputDir directory is searched for datasets.
        Time and redshift information are gathered by temporarily
        instantiating each dataset.  This can be used when simulation
        data was created in a non-standard way, making it difficult
        to guess the corresponding time and redshift information.
        Default: False.

    Examples
    --------
    >>> import yt
    >>> gs = yt.load_simulation("my_simulation.par", "Gadget")
    >>> gs.get_time_series()
    >>> for ds in gs:
    ...     print(ds.current_time)

    """

    def __init__(self, parameter_filename, find_outputs=False):
        self.simulation_type = "particle"
        self.dimensionality = 3
        SimulationTimeSeries.__init__(
            self, parameter_filename, find_outputs=find_outputs
        )

    def _set_units(self):
        self.unit_registry = UnitRegistry()
        self.time_unit = self.quan(1.0, "s")
        if self.cosmological_simulation:
            # Instantiate Cosmology object for units and time conversions.
            self.cosmology = Cosmology(
                hubble_constant=self.hubble_constant,
                omega_matter=self.omega_matter,
                omega_lambda=self.omega_lambda,
                unit_registry=self.unit_registry,
            )
            if "h" in self.unit_registry:
                self.unit_registry.modify("h", self.hubble_constant)
            else:
                self.unit_registry.add(
                    "h", self.hubble_constant, dimensions.dimensionless
                )
            # Comoving lengths
            for my_unit in ["m", "pc", "AU"]:
                new_unit = f"{my_unit}cm"
                # technically not true, but should be ok
                self.unit_registry.add(
                    new_unit,
                    self.unit_registry.lut[my_unit][0],
                    dimensions.length,
                    "\\rm{%s}/(1+z)" % my_unit,
                    prefixable=True,
                )
            self.length_unit = self.quan(
                self.unit_base["UnitLength_in_cm"],
                "cmcm / h",
                registry=self.unit_registry,
            )
            self.mass_unit = self.quan(
                self.unit_base["UnitMass_in_g"], "g / h", registry=self.unit_registry
            )
            self.box_size = self.box_size * self.length_unit
            self.domain_left_edge = self.domain_left_edge * self.length_unit
            self.domain_right_edge = self.domain_right_edge * self.length_unit
            self.unit_registry.add(
                "unitary",
                float(self.box_size.in_base()),
                self.length_unit.units.dimensions,
            )
        else:
            # Read time from file for non-cosmological sim
            self.time_unit = self.quan(
                self.unit_base["UnitLength_in_cm"]
                / self.unit_base["UnitVelocity_in_cm_per_s"],
                "s",
            )
            self.unit_registry.add("code_time", 1.0, dimensions.time)
            self.unit_registry.modify("code_time", self.time_unit)
            # Length
            self.length_unit = self.quan(self.unit_base["UnitLength_in_cm"], "cm")
            self.unit_registry.add("code_length", 1.0, dimensions.length)
            self.unit_registry.modify("code_length", self.length_unit)

    def get_time_series(
        self,
        initial_time=None,
        final_time=None,
        initial_redshift=None,
        final_redshift=None,
        times=None,
        redshifts=None,
        tolerance=None,
        parallel=True,
        setup_function=None,
    ):

        """
        Instantiate a DatasetSeries object for a set of outputs.

        If no additional keywords given, a DatasetSeries object will be
        created with all potential datasets created by the simulation.

        Outputs can be gather by specifying a time or redshift range
        (or combination of time and redshift), with a specific list of
        times or redshifts), or by simply searching all subdirectories
        within the simulation directory.

        initial_time : tuple of type (float, str)
            The earliest time for outputs to be included.  This should be
            given as the value and the string representation of the units.
            For example, (5.0, "Gyr").  If None, the initial time of the
            simulation is used.  This can be used in combination with
            either final_time or final_redshift.
            Default: None.
        final_time : tuple of type (float, str)
            The latest time for outputs to be included.  This should be
            given as the value and the string representation of the units.
            For example, (13.7, "Gyr"). If None, the final time of the
            simulation is used.  This can be used in combination with either
            initial_time or initial_redshift.
            Default: None.
        times : tuple of type (float array, str)
            A list of times for which outputs will be found and the units
            of those values.  For example, ([0, 1, 2, 3], "s").
            Default: None.
        initial_redshift : float
            The earliest redshift for outputs to be included.  If None,
            the initial redshift of the simulation is used.  This can be
            used in combination with either final_time or
            final_redshift.
            Default: None.
        final_redshift : float
            The latest redshift for outputs to be included.  If None,
            the final redshift of the simulation is used.  This can be
            used in combination with either initial_time or
            initial_redshift.
            Default: None.
        redshifts : array_like
            A list of redshifts for which outputs will be found.
            Default: None.
        tolerance : float
            Used in combination with "times" or "redshifts" keywords,
            this is the tolerance within which outputs are accepted
            given the requested times or redshifts.  If None, the
            nearest output is always taken.
            Default: None.
        parallel : bool/int
            If True, the generated DatasetSeries will divide the work
            such that a single processor works on each dataset.  If an
            integer is supplied, the work will be divided into that
            number of jobs.
            Default: True.
        setup_function : callable, accepts a ds
            This function will be called whenever a dataset is loaded.

        Examples
        --------

        >>> import yt
        >>> gs = yt.load_simulation("my_simulation.par", "Gadget")

        >>> gs.get_time_series(initial_redshift=10, final_time=(13.7, "Gyr"))

        >>> gs.get_time_series(redshifts=[3, 2, 1, 0])

        >>> # after calling get_time_series
        >>> for ds in gs.piter():
        ...     p = ProjectionPlot(ds, "x", ("gas", "density"))
        ...     p.save()

        >>> # An example using the setup_function keyword
        >>> def print_time(ds):
        ...     print(ds.current_time)
        >>> gs.get_time_series(setup_function=print_time)
        >>> for ds in gs:
        ...     SlicePlot(ds, "x", "Density").save()

        """

        if (
            initial_redshift is not None or final_redshift is not None
        ) and not self.cosmological_simulation:
            raise InvalidSimulationTimeSeries(
                "An initial or final redshift has been given for a "
                + "noncosmological simulation."
            )

        my_all_outputs = self.all_outputs
        if not my_all_outputs:
            DatasetSeries.__init__(
                self, outputs=[], parallel=parallel, unit_base=self.unit_base
            )
            mylog.info("0 outputs loaded into time series.")
            return

        # Apply selection criteria to the set.
        if times is not None:
            my_outputs = self._get_outputs_by_key(
                "time", times, tolerance=tolerance, outputs=my_all_outputs
            )

        elif redshifts is not None:
            my_outputs = self._get_outputs_by_key(
                "redshift", redshifts, tolerance=tolerance, outputs=my_all_outputs
            )

        else:
            if initial_time is not None:
                if isinstance(initial_time, float):
                    initial_time = self.quan(initial_time, "code_time")
                elif isinstance(initial_time, tuple) and len(initial_time) == 2:
                    initial_time = self.quan(*initial_time)
                elif not isinstance(initial_time, unyt_array):
                    raise RuntimeError(
                        "Error: initial_time must be given as a float or "
                        + "tuple of (value, units)."
                    )
            elif initial_redshift is not None:
                my_initial_time = self.cosmology.t_from_z(initial_redshift)
            else:
                my_initial_time = self.initial_time

            if final_time is not None:
                if isinstance(final_time, float):
                    final_time = self.quan(final_time, "code_time")
                elif isinstance(final_time, tuple) and len(final_time) == 2:
                    final_time = self.quan(*final_time)
                elif not isinstance(final_time, unyt_array):
                    raise RuntimeError(
                        "Error: final_time must be given as a float or "
                        + "tuple of (value, units)."
                    )
                my_final_time = final_time.in_units("s")
            elif final_redshift is not None:
                my_final_time = self.cosmology.t_from_z(final_redshift)
            else:
                my_final_time = self.final_time

            my_initial_time.convert_to_units("s")
            my_final_time.convert_to_units("s")
            my_times = np.array([a["time"] for a in my_all_outputs])
            my_indices = np.digitize([my_initial_time, my_final_time], my_times)
            if my_initial_time == my_times[my_indices[0] - 1]:
                my_indices[0] -= 1
            my_outputs = my_all_outputs[my_indices[0] : my_indices[1]]

        init_outputs = []
        for output in my_outputs:
            if os.path.exists(output["filename"]):
                init_outputs.append(output["filename"])
        if len(init_outputs) == 0 and len(my_outputs) > 0:
            mylog.warning(
                "Could not find any datasets.  "
                "Check the value of OutputDir in your parameter file."
            )

        DatasetSeries.__init__(
            self,
            outputs=init_outputs,
            parallel=parallel,
            setup_function=setup_function,
            unit_base=self.unit_base,
        )
        mylog.info("%d outputs loaded into time series.", len(init_outputs))

    def _parse_parameter_file(self):
        """
        Parses the parameter file and establishes the various
        dictionaries.
        """

        self.unit_base = {}

        # Let's read the file
        lines = open(self.parameter_filename).readlines()
        comments = ["%", ";"]
        for line in (l.strip() for l in lines):
            for comment in comments:
                if comment in line:
                    line = line[0 : line.find(comment)]
            if len(line) < 2:
                continue
            param, vals = (i.strip() for i in line.split(None, 1))
            # First we try to decipher what type of value it is.
            vals = vals.split()
            # Special case approaching.
            if "(do" in vals:
                vals = vals[:1]
            if len(vals) == 0:
                pcast = str  # Assume NULL output
            else:
                v = vals[0]
                # Figure out if it's castable to floating point:
                try:
                    float(v)
                except ValueError:
                    pcast = str
                else:
                    if any("." in v or "e" in v for v in vals):
                        pcast = float
                    elif v == "inf":
                        pcast = str
                    else:
                        pcast = int
            # Now we figure out what to do with it.
            if param.startswith("Unit"):
                self.unit_base[param] = float(vals[0])
            if len(vals) == 0:
                vals = ""
            elif len(vals) == 1:
                vals = pcast(vals[0])
            else:
                vals = np.array([pcast(i) for i in vals])

            self.parameters[param] = vals

        # Domain dimensions for Gadget datasets are always 2x2x2 for octree
        self.domain_dimensions = np.array([2, 2, 2])

        if self.parameters["ComovingIntegrationOn"]:
            cosmo_attr = {
                "box_size": "BoxSize",
                "omega_lambda": "OmegaLambda",
                "omega_matter": "Omega0",
                "hubble_constant": "HubbleParam",
            }
            self.initial_redshift = 1.0 / self.parameters["TimeBegin"] - 1.0
            self.final_redshift = 1.0 / self.parameters["TimeMax"] - 1.0
            self.cosmological_simulation = 1
            for a, v in cosmo_attr.items():
                if v not in self.parameters:
                    raise MissingParameter(self.parameter_filename, v)
                setattr(self, a, self.parameters[v])
            self.domain_left_edge = np.array([0.0, 0.0, 0.0])
            self.domain_right_edge = (
                np.array([1.0, 1.0, 1.0]) * self.parameters["BoxSize"]
            )
        else:
            self.cosmological_simulation = 0
            self.omega_lambda = self.omega_matter = self.hubble_constant = 0.0

    def _find_data_dir(self):
        """
        Find proper location for datasets.  First look where parameter file
        points, but if this doesn't exist then default to the current
        directory.
        """
        if self.parameters["OutputDir"].startswith("/"):
            data_dir = self.parameters["OutputDir"]
        else:
            data_dir = os.path.join(self.directory, self.parameters["OutputDir"])
        if not os.path.exists(data_dir):
            mylog.info(
                "OutputDir not found at %s, instead using %s.", data_dir, self.directory
            )
            data_dir = self.directory
        self.data_dir = data_dir

    def _snapshot_format(self, index=None):
        """
        The snapshot filename for a given index.  Modify this for different
        naming conventions.
        """

        if self.parameters["NumFilesPerSnapshot"] > 1:
            suffix = ".0"
        else:
            suffix = ""
        if self.parameters["SnapFormat"] == 3:
            suffix += ".hdf5"
        if index is None:
            count = "*"
        else:
            count = "%03d" % index
        filename = f"{self.parameters['SnapshotFileBase']}_{count}{suffix}"
        return os.path.join(self.data_dir, filename)

    def _get_all_outputs(self, find_outputs=False):
        """
        Get all potential datasets and combine into a time-sorted list.
        """

        # Find the data directory where the outputs are
        self._find_data_dir()

        # Create the set of outputs from which further selection will be done.
        if find_outputs:
            self._find_outputs()
        else:
            if self.parameters["OutputListOn"]:
                a_values = [
                    float(a)
                    for a in open(
                        os.path.join(
                            self.data_dir, self.parameters["OutputListFilename"]
                        ),
                    ).readlines()
                ]
            else:
                a_values = [float(self.parameters["TimeOfFirstSnapshot"])]
                time_max = float(self.parameters["TimeMax"])
                while a_values[-1] < time_max:
                    if self.cosmological_simulation:
                        a_values.append(
                            a_values[-1] * self.parameters["TimeBetSnapshot"]
                        )
                    else:
                        a_values.append(
                            a_values[-1] + self.parameters["TimeBetSnapshot"]
                        )
                if a_values[-1] > time_max:
                    a_values[-1] = time_max

            if self.cosmological_simulation:
                self.all_outputs = [
                    {"filename": self._snapshot_format(i), "redshift": (1.0 / a - 1)}
                    for i, a in enumerate(a_values)
                ]

                # Calculate times for redshift outputs.
                for output in self.all_outputs:
                    output["time"] = self.cosmology.t_from_z(output["redshift"])
            else:
                self.all_outputs = [
                    {
                        "filename": self._snapshot_format(i),
                        "time": self.quan(a, "code_time"),
                    }
                    for i, a in enumerate(a_values)
                ]

            self.all_outputs.sort(key=lambda obj: obj["time"].to_ndarray())

    def _calculate_simulation_bounds(self):
        """
        Figure out the starting and stopping time and redshift for the simulation.
        """

        # Convert initial/final redshifts to times.
        if self.cosmological_simulation:
            self.initial_time = self.cosmology.t_from_z(self.initial_redshift)
            self.initial_time.units.registry = self.unit_registry
            self.final_time = self.cosmology.t_from_z(self.final_redshift)
            self.final_time.units.registry = self.unit_registry

        # If not a cosmology simulation, figure out the stopping criteria.
        else:
            if "TimeBegin" in self.parameters:
                self.initial_time = self.quan(self.parameters["TimeBegin"], "code_time")
            else:
                self.initial_time = self.quan(0.0, "code_time")

            if "TimeMax" in self.parameters:
                self.final_time = self.quan(self.parameters["TimeMax"], "code_time")
            else:
                self.final_time = None
            if "TimeMax" not in self.parameters:
                raise NoStoppingCondition(self.parameter_filename)

    def _find_outputs(self):
        """
        Search for directories matching the data dump keywords.
        If found, get dataset times py opening the ds.
        """
        potential_outputs = glob.glob(self._snapshot_format())
        self.all_outputs = self._check_for_outputs(potential_outputs)
        self.all_outputs.sort(key=lambda obj: obj["time"])
        only_on_root(mylog.info, "Located %d total outputs.", len(self.all_outputs))

        # manually set final time and redshift with last output
        if self.all_outputs:
            self.final_time = self.all_outputs[-1]["time"]
            if self.cosmological_simulation:
                self.final_redshift = self.all_outputs[-1]["redshift"]

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
            my_storage.result = {
                "filename": output,
                "time": ds.current_time.in_units("s"),
            }
            if ds.cosmological_simulation:
                my_storage.result["redshift"] = ds.current_redshift

        my_outputs = [
            my_output for my_output in my_outputs.values() if my_output is not None
        ]
        return my_outputs

    def _write_cosmology_outputs(self, filename, outputs, start_index, decimals=3):
        r"""
        Write cosmology output parameters for a cosmology splice.
        """

        mylog.info("Writing redshift output list to %s.", filename)
        f = open(filename, "w")
        for output in outputs:
            f.write(f"{1.0 / (1.0 + output['redshift']):f}\n")
        f.close()
