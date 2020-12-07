import os
import weakref
from pathlib import Path

import numpy as np

from yt.data_objects.index_subobjects.grid_patch import AMRGridPatch
from yt.data_objects.static_output import Dataset
from yt.frontends.idefix.dmpfile_io import read_idefix_dmpfile
from yt.frontends.idefix.inifile_io import read_idefix_inifile
from yt.funcs import setdefaultattr
from yt.geometry.grid_geometry_handler import GridIndex
from yt.utilities.logger import ytLogger

from .fields import IdefixFieldInfo


class IdefixGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, id, index, level):
        super(IdefixGrid, self).__init__(id, filename=index.index_filename, index=index)
        self.Parent = None
        self.Children = []
        self.Level = level

    def __repr__(self):
        return "IdefixGrid_%04i (%s)" % (self.id, self.ActiveDimensions)


class IdefixHierarchy(GridIndex):
    grid = IdefixGrid

    def __init__(self, ds, dataset_type="idefix"):
        self.dataset_type = dataset_type
        self.dataset = weakref.proxy(ds)
        # for now, the index file is the dataset!
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        # float type for the simulation edges and must be float64 now
        self.float_type = np.float64
        super(IdefixHierarchy, self).__init__(ds, dataset_type)

    def _detect_output_fields(self):
        # This needs to set a self.field_list that contains all the available,
        # on-disk fields. No derived fields should be defined here.
        # NOTE: Each should be a tuple, where the first element is the on-disk
        # fluid type or particle type.  Convention suggests that the on-disk
        # fluid type is usually the dataset_type and the on-disk particle type
        # (for a single population of particles) is "io".
        pass

    def _count_grids(self):
        # This needs to set self.num_grids (int)
        pass

    def _parse_index(self):
        # This needs to fill the following arrays, where N is self.num_grids:
        #   self.grid_left_edge         (N, 3) <= float64
        #   self.grid_right_edge        (N, 3) <= float64
        #   self.grid_dimensions        (N, 3) <= int
        #   self.grid_particle_count    (N, 1) <= int
        #   self.grid_levels            (N, 1) <= int
        #   self.grids                  (N, 1) <= grid objects
        #   self.max_level = self.grid_levels.max()
        pass

    def _populate_grid_objects(self):
        # the minimal form of this method is
        #
        # for g in self.grids:
        #     g._prepare_grid()
        #     g._setup_dx()
        #
        # This must also set:
        #   g.Children <= list of child grids
        #   g.Parent   <= parent grid
        # This is handled by the frontend because often the children must be identified.
        pass


class IdefixDataset(Dataset):
    _index_class = IdefixHierarchy
    _field_info_class = IdefixFieldInfo

    def __init__(
        self,
        filename,
        dataset_type="idefix",
        storage_filename=None,
        units_override=None,
        inifile=None,
    ):
        self.fluid_types += ("idefix",)
        self._inifile_name = inifile
        super(IdefixDataset, self).__init__(
            filename, dataset_type, units_override=units_override
        )
        self.storage_filename = storage_filename

        # idefix does not support grid refinement
        self.refine_by = 1

    def _set_code_unit_attributes(self):
        # This is where quantities are created that represent the various
        # on-disk units.  These are the currently available quantities which
        # should be set, along with examples of how to set them to standard
        # values.
        #
        # self.length_unit = self.quan(1.0, "cm")
        # self.mass_unit = self.quan(1.0, "g")
        # self.time_unit = self.quan(1.0, "s")
        # self.time_unit = self.quan(1.0, "s")
        #
        # These can also be set:
        # self.velocity_unit = self.quan(1.0, "cm/s")
        # self.magnetic_unit = self.quan(1.0, "gauss")
        for key, unit in self.__class__.default_units.items():
            setdefaultattr(self, key, self.quan(1, unit))

    def _parse_parameter_file(self):
        # first pass in the dmpfile: read everything except large arrays
        fprops, fdata = read_idefix_dmpfile(self.parameter_filename, skip_arrays=True)

        # parse the grid
        axes = ("x1", "x2", "x3")
        self.domain_dimensions = np.concatenate([fprops[k][-1] for k in axes])
        self.dimensionality = np.count_nonzero(self.domain_dimensions - 1)
        domain_half_width = np.array([fdata[k] for k in axes], dtype="float64") / 2
        self.domain_left_edge = -domain_half_width
        self.domain_right_edge = +domain_half_width

        self.current_time = fdata["time"]

        # ... this section needs to be improved
        if self._inifile_name is not None:
            self.parameters.update(read_idefix_inifile(self._inifile_name))
            if any(self.parameters["Grid"][f"X{i}-grid"] != "u" for i in "123"):
                ytLogger.warning(
                    "Currently, yt does not support non-uniform grid stepping."
                )
        self.periodicity = (True, True, True)
        self.geometry = "cartesian"
        # ...

        # idefix is never cosmological
        self.cosmological_simulation = 0
        self.current_redshift = 0.0
        self.omega_lambda = 0.0
        self.omega_matter = 0.0
        self.hubble_constant = 0.0

    @classmethod
    def _is_valid(self, fn, *args, **kwargs):
        # a stupid heuristic test
        name = Path(fn).name
        return name.startswith("dump") and name.endswith(".dmp")
