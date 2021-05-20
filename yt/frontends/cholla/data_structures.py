import os
import weakref

import numpy as np

from yt.data_objects.index_subobjects.grid_patch import AMRGridPatch
from yt.data_objects.static_output import Dataset
from yt.funcs import setdefaultattr
from yt.geometry.grid_geometry_handler import GridIndex
from yt.utilities.on_demand_imports import _h5py as h5py

from .fields import ChollaFieldInfo


class ChollaGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, id, index, level):
        super().__init__(id, filename=index.index_filename, index=index)
        self.Parent = None
        self.Children = []
        self.Level = level

    def __repr__(self):
        return "ChollaGrid_%04i (%s)" % (self.id, self.ActiveDimensions)


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
        f = open(self.index_filename, "rb")
        self.field_list = [("cholla", k) for k in f.keys()]
        f.close()

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


class ChollaDataset(Dataset):
    _index_class = ChollaHierarchy
    _field_info_class = ChollaFieldInfo

    def __init__(
        self,
        filename,
        dataset_type="cholla",
        storage_filename=None,
        units_override=None,
    ):
        self.fluid_types += ("cholla",)
        super().__init__(filename, dataset_type, units_override=units_override)
        self.storage_filename = storage_filename

    def _set_code_unit_attributes(self):
        # This is where quantities are created that represent the various
        # on-disk units.  These are the currently available quantities which
        # should be set, along with examples of how to set them to standard
        # values.
        #
        self.length_unit = self.quan(1.0, "pc")
        self.mass_unit = self.quan(1.0, "Msun")
        self.time_unit = self.quan(1000, "yr")
        self.velocity_unit = self.quan(1.0, "cm/s")
        self.magnetic_unit = self.quan(1.0, "gauss")

        # this minimalistic implementation fills the requirements for
        # this frontend to run, change it to make it run _correctly_ !
        for key, unit in self.__class__.default_units.items():
            setdefaultattr(self, key, self.quan(1, unit))

    def _parse_parameter_file(self):

        h5f = h5py.File(self.parameter_filename, mode="r")
        attrs = h5f.attrs
        self.parameters = {k: v for (k, v) in attrs.items()}
        self.domain_left_edge = attrs["bounds"][:]
        self.domain_right_edge = attrs["domain"][:]
        self.dimensionality = len(attrs["dims"][:])
        self.domain_dimensions = attrs["dims"][:]
        self.current_time = attrs["t"][:]
        h5f.close()

        # CHOLLA cannot yet be run as a cosmological simulation
        self.cosmological_simulation = 0
        self.current_redshift = 0
        self.omega_lambda = 0
        self.omega_matter = 0
        self.hubble_constant = 0

        # CHOLLA datasets are always unigrid cartesian
        self.geometry = "cartesian"

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        # This accepts a filename or a set of arguments and returns True or
        # False depending on if the file is of the type requested.
        #
        # The functionality in this method should be unique enough that it can
        # differentiate the frontend from others. Sometimes this means looking
        # for specific fields or attributes in the dataset in addition to
        # looking at the file name or extension.
        try:
            fileh = h5py.File(filename, mode="r")
            attrs = fileh.attrs
            if "bounds" in attrs and "domain" in attrs:
                fileh.close()
                return True
            fileh.close()
        except Exception:
            pass
        return False
