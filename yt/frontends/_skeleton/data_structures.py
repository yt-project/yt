import os
import weakref

import numpy as np

from yt.data_objects.index_subobjects.grid_patch import AMRGridPatch
from yt.data_objects.static_output import Dataset
from yt.funcs import setdefaultattr
from yt.geometry.grid_geometry_handler import GridIndex

from .fields import SkeletonFieldInfo


class SkeletonGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, id, index, level):
        super().__init__(id, filename=index.index_filename, index=index)
        self.Parent = None
        self.Children = []
        self.Level = level


class SkeletonHierarchy(GridIndex):
    grid = SkeletonGrid

    def __init__(self, ds, dataset_type="skeleton"):
        self.dataset_type = dataset_type
        self.dataset = weakref.proxy(ds)
        # for now, the index file is the dataset!
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        # float type for the simulation edges and must be float64 now
        self.float_type = np.float64
        super().__init__(ds, dataset_type)

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


class SkeletonDataset(Dataset):
    _index_class = SkeletonHierarchy
    _field_info_class = SkeletonFieldInfo

    def __init__(
        self,
        filename,
        dataset_type="skeleton",
        storage_filename=None,
        units_override=None,
        unit_system="cgs",
        default_species_fields=None,
    ):
        self.fluid_types += ("skeleton",)
        super().__init__(
            filename,
            dataset_type,
            units_override=units_override,
            unit_system=unit_system,
            default_species_fields=default_species_fields,
        )
        self.storage_filename = storage_filename
        # refinement factor between a grid and its subgrid
        # self.refine_by = 2

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
        #
        # If your frontend uses SI EM units, set magnetic units like this
        # instead:
        # self.magnetic_unit = self.quan(1.0, "T")

        # this minimalistic implementation fills the requirements for
        # this frontend to run, change it to make it run _correctly_ !
        for key, unit in self.__class__.default_units.items():
            setdefaultattr(self, key, self.quan(1, unit))

    def _parse_parameter_file(self):
        # This needs to set up the following items.  Note that these are all
        # assumed to be in code units; domain_left_edge and domain_right_edge
        # will be converted to YTArray automatically at a later time.
        # This includes the cosmological parameters.
        #
        #   self.unique_identifier      <= unique identifier for the dataset
        #                                  being read (e.g., UUID or ST_CTIME)
        #   self.parameters             <= dict full of code-specific items of use
        #   self.domain_left_edge       <= three-element array of float64
        #   self.domain_right_edge      <= three-element array of float64
        #   self.dimensionality         <= int
        #   self.domain_dimensions      <= three-element array of int64
        #   self.periodicity            <= three-element tuple of booleans
        #   self.current_time           <= simulation time in code units (float)
        #
        # We also set up cosmological information.  Set these to zero if
        # non-cosmological.
        #
        #   self.cosmological_simulation    <= int, 0 or 1
        #   self.current_redshift           <= float
        #   self.omega_lambda               <= float
        #   self.omega_matter               <= float
        #   self.hubble_constant            <= float

        # optional (the following have default implementations)
        #
        #   self.geometry  <= a lower case string
        #                     ("cartesian", "polar", "cylindrical"...)
        #                     (defaults to 'cartesian')

        # this attribute is required.
        # Change this value to a constant 0 if time is not relevant to your dataset.
        # Otherwise, parse its value in any appropriate fashion.
        self.current_time = -1

        # required. Change this if need be.
        self.cosmological_simulation = 0

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        # This accepts a filename or a set of arguments and returns True or
        # False depending on if the file is of the type requested.
        #
        # The functionality in this method should be unique enough that it can
        # differentiate the frontend from others. Sometimes this means looking
        # for specific fields or attributes in the dataset in addition to
        # looking at the file name or extension.
        return False
