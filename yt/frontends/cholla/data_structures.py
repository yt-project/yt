import os
import weakref

import numpy as np

from yt.data_objects.index_subobjects.grid_patch import AMRGridPatch
from yt.data_objects.static_output import Dataset
from yt.funcs import get_pbar, setdefaultattr
from yt.geometry.api import Geometry
from yt.geometry.grid_geometry_handler import GridIndex
from yt.utilities.on_demand_imports import _h5py as h5py

from .fields import ChollaFieldInfo


def _split_fname_proc_suffix(filename: str):
    """Splits ``filename`` at the '.' separating the beginning part of the
    string from the process-id suffix, and returns both parts in a 2-tuple.

    When cholla is compiled with MPI and it directly writes data-files, each
    process appends a suffix to each filename that denotes the process-id. For
    example, the MPI-compiled version might write '0.h5.0'. If that function is
    passed such a string, then it returns ``('0.h5', '0')``.

    In cases where there is no suffix, the output is ``(filename, '')``. This
    might come up if the user concatenated the output files, which is common
    practice.
    """

    # at this time, we expect the suffix to be the minimum number of characters
    # that are necessary to represent the process id. For flexibility, we will
    # allow extra zero-padding

    sep_i = filename.rfind(".")
    suf_len = len(filename) - (sep_i + 1)
    if (sep_i == -1) or (suf_len == 0) or not filename[sep_i + 1 :].isdecimal():
        return (filename, "")
    elif (sep_i == 0) or ((sep_i - 1) == filename.rfind("/")):
        raise ValueError(
            f"can't split a process-suffix off of {filename!r} "
            "since the remaining filename would be empty"
        )
    else:
        return (filename[:sep_i], filename[sep_i + 1 :])


class ChollaGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, id, index, level, dims, filename):
        super().__init__(id, filename=filename, index=index)
        self.Parent = None
        self.Children = []
        self.Level = level
        self.ActiveDimensions = dims


class ChollaHierarchy(GridIndex):
    grid = ChollaGrid

    def __init__(self, ds, dataset_type="cholla"):
        self.dataset_type = dataset_type
        self.dataset = weakref.proxy(ds)
        # for now, the index file is the dataset!
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        # float type for the simulation edges and must be float64 now
        self.float_type = np.float64
        super().__init__(ds, dataset_type)

    def _detect_output_fields(self):
        with h5py.File(self.index_filename, mode="r") as h5f:
            self.field_list = [("cholla", k) for k in h5f.keys()]

    def _count_grids(self):
        # the number of grids is equal to the number of processes, unless the
        # dataset has been concatenated. But, when the dataset is concatenated
        # (a common post-processing step), the "nprocs" hdf5 attribute is
        # usually dropped.

        with h5py.File(self.index_filename, mode="r") as h5f:
            nprocs = h5f.attrs.get("nprocs", np.array([1, 1, 1]))[:].astype("=i8")
        self.num_grids = np.prod(nprocs)

        if self.num_grids > 1:
            # When there's more than 1 grid, we expect the user to
            # - have not changed the names of the output files
            # - have passed the file written by process 0 to ``yt.load``
            # Let's perform a sanity-check that self.index_filename has the
            # expected suffix for a file written by mpi-process 0
            if int(_split_fname_proc_suffix(self.index_filename)[1]) != 0:
                raise ValueError(
                    "the primary file associated with a "
                    "distributed cholla dataset must end in '.0'"
                )

    def _parse_index(self):
        self.grids = np.empty(self.num_grids, dtype="object")

        # construct an iterable over the pairs of grid-index and corresponding
        # filename
        if self.num_grids == 1:
            ind_fname_pairs = [(0, self.index_filename)]
        else:
            # index_fname should has the form f'{self.directory}/<prefix>.0'
            # strip off the '.0' and determine the contents of <prefix>
            pref, suf = _split_fname_proc_suffix(self.index_filename)
            assert int(suf) == 0  # sanity check!

            ind_fname_pairs = ((i, f"{pref}.{i}") for i in range(self.num_grids))

        dims_global = self.ds.domain_dimensions[:]
        pbar = get_pbar("Parsing Hierarchy", self.num_grids)

        # It would be nice if we could avoid reading in every hdf5 file during
        # this step... (to do this, Cholla could probably encode how the blocks
        # are sorted in an hdf5 attribute)

        for i, fname in ind_fname_pairs:
            if self.num_grids == 1:
                # if the file was concatenated, we might be missing attributes
                # that are accessed in the other branch. To avoid issues, we use
                # hardcoded values
                left_frac, right_frac, dims_local = 0.0, 1.0, dims_global
            else:
                with h5py.File(fname, "r") as f:
                    offset = f.attrs["offset"][:].astype("=i8")
                    dims_local = f.attrs["dims_local"][:].astype("=i8")
                left_frac = offset / dims_global
                right_frac = (offset + dims_local) / dims_global

            level = 0

            self.grids[i] = self.grid(
                i,
                index=self,
                level=level,
                dims=dims_local,
                filename=fname,
            )

            self.grid_left_edge[i] = left_frac
            self.grid_right_edge[i] = right_frac
            self.grid_dimensions[i] = dims_local
            self.grid_levels[i, 0] = level
            self.grid_particle_count[i, 0] = 0

            pbar.update(i + 1)
        pbar.finish()

        slope = self.ds.domain_width / self.ds.arr(np.ones(3), "code_length")
        self.grid_left_edge = self.grid_left_edge * slope + self.ds.domain_left_edge
        self.grid_right_edge = self.grid_right_edge * slope + self.ds.domain_left_edge

        self.max_level = 0

    def _populate_grid_objects(self):
        for i in range(self.num_grids):
            g = self.grids[i]
            g._prepare_grid()
            g._setup_dx()


class ChollaDataset(Dataset):
    _load_requirements = ["h5py"]
    _index_class = ChollaHierarchy
    _field_info_class = ChollaFieldInfo

    def __init__(
        self,
        filename,
        dataset_type="cholla",
        storage_filename=None,
        units_override=None,
        unit_system="cgs",
    ):
        self.fluid_types += ("cholla",)
        super().__init__(filename, dataset_type, units_override=units_override)
        self.storage_filename = storage_filename

    def _set_code_unit_attributes(self):
        # This is where quantities are created that represent the various
        # on-disk units.  These are the defaults, but if they are listed
        # in the HDF5 attributes for a file, which is loaded first, then those are
        # used instead.
        #
        if not self.length_unit:
            self.length_unit = self.quan(1.0, "pc")
        if not self.mass_unit:
            self.mass_unit = self.quan(1.0, "Msun")
        if not self.time_unit:
            self.time_unit = self.quan(1000, "yr")
        if not self.velocity_unit:
            self.velocity_unit = self.quan(1.0, "cm/s")
        if not self.magnetic_unit:
            self.magnetic_unit = self.quan(1.0, "gauss")

        for key, unit in self.__class__.default_units.items():
            setdefaultattr(self, key, self.quan(1, unit))

    def _parse_parameter_file(self):
        with h5py.File(self.parameter_filename, mode="r") as h5f:
            attrs = h5f.attrs
            self.parameters = dict(attrs.items())
            self.domain_left_edge = attrs["bounds"][:].astype("=f8")
            self.domain_right_edge = self.domain_left_edge + attrs["domain"][:].astype(
                "=f8"
            )
            self.dimensionality = len(attrs["dims"][:])
            self.domain_dimensions = attrs["dims"][:].astype("=i8")
            self.current_time = attrs["t"][:]
            self._periodicity = tuple(attrs.get("periodicity", (False, False, False)))
            self.gamma = attrs.get("gamma", 5.0 / 3.0)
            self.mu = attrs.get("mu", 1.0)
            self.refine_by = 1

            # If header specifies code units, default to those (in CGS)
            length_unit = attrs.get("length_unit", None)
            mass_unit = attrs.get("mass_unit", None)
            time_unit = attrs.get("time_unit", None)
            velocity_unit = attrs.get("velocity_unit", None)
            magnetic_unit = attrs.get("magnetic_unit", None)
            if length_unit:
                self.length_unit = self.quan(length_unit[0], "cm")
            if mass_unit:
                self.mass_unit = self.quan(mass_unit[0], "g")
            if time_unit:
                self.time_unit = self.quan(time_unit[0], "s")
            if velocity_unit:
                self.velocity_unit = self.quan(velocity_unit[0], "cm/s")
            if magnetic_unit:
                self.magnetic_unit = self.quan(magnetic_unit[0], "gauss")

            # this minimalistic implementation fills the requirements for
            # this frontend to run, change it to make it run _correctly_ !
            for key, unit in self.__class__.default_units.items():
                setdefaultattr(self, key, self.quan(1, unit))

        # CHOLLA cannot yet be run as a cosmological simulation
        self.cosmological_simulation = 0
        self.current_redshift = 0.0
        self.omega_lambda = 0.0
        self.omega_matter = 0.0
        self.hubble_constant = 0.0

        # CHOLLA datasets are always unigrid cartesian
        self.geometry = Geometry.CARTESIAN

    @classmethod
    def _is_valid(cls, filename: str, *args, **kwargs) -> bool:
        # This accepts a filename or a set of arguments and returns True or
        # False depending on if the file is of the type requested.
        if cls._missing_load_requirements():
            return False

        try:
            fileh = h5py.File(filename, mode="r")
        except OSError:
            return False

        try:
            attrs = fileh.attrs
        except AttributeError:
            return False
        else:
            return (
                "bounds" in attrs
                and "domain" in attrs
                and attrs.get("data_type") != "yt_light_ray"
            )
        finally:
            fileh.close()
