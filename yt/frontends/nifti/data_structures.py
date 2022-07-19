import os
import weakref

import numpy as np

from yt.data_objects.index_subobjects.grid_patch import AMRGridPatch
from yt.data_objects.static_output import Dataset
from yt.funcs import setdefaultattr
from yt.geometry.grid_geometry_handler import GridIndex
from yt.utilities.on_demand_imports import _nibabel as nib, NotAModule
from .fields import NiftiFieldInfo


class NiftiGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, id, index, level):
        super().__init__(id, filename=index.index_filename, index=index)
        self.Parent = None
        self.Children = []
        self.Level = level


class NiftiHierarchy(GridIndex):
    grid = NiftiGrid

    def __init__(self, ds, dataset_type="nifti"):
        self.dataset_type = dataset_type
        self.dataset = weakref.proxy(ds)
        # for now, the index file is the dataset!
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        # float type for the simulation edges and must be float64 now
        self.float_type = np.float64
        super().__init__(ds, dataset_type)

    def _detect_output_fields(self):
        self.field_list = [("nifti", "intensity")]

    def _count_grids(self):
        self.num_grids = 1

    def _parse_index(self):
        self.grid_left_edge[0][:] = self.ds.domain_left_edge[:]
        self.grid_right_edge[0][:] = self.ds.domain_right_edge[:]
        self.grid_dimensions[0][:] = self.ds.domain_dimensions[:]
        self.grid_particle_count[0][0] = 0
        self.grid_levels[0][0] = 0
        self.max_level = 0

    def _populate_grid_objects(self):
        # only a single grid, no need to loop
        g = self.grid(0, self, self.grid_levels.flat[0], self.grid_dimensions[0])
        g._prepare_grid()
        g._setup_dx()
        self.grids = np.array([g], dtype="object")



class NiftiDataset(Dataset):
    _index_class = NiftiHierarchy
    _field_info_class = NiftiFieldInfo

    def __init__(
        self,
        filename,
        dataset_type="nifti",
        storage_filename=None,
        units_override=None,
        unit_system="cgs",
        default_species_fields=None,
    ):
        self.fluid_types += ("nifti",)
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
        img = nib.load(self.filename)
        self.parameters = {}

        # This needs to set up the following items.  Note that these are all
        # assumed to be in code units; domain_left_edge and domain_right_edge
        # will be converted to YTArray automatically at a later time.
        # This includes the cosmological parameters.
        #
        
        self.parameters['header'] = img.header
        img_data_shape = img.get_fdata().shape
        self.domain_left_edge = np.array((0,0,0))
        self.domain_right_edge = np.array(img_data_shape)
        self.dimensionality = 3
        self.domain_dimensions = img_data_shape
        self.periodicity = (False, False, False)
        self.current_time = 0
        #
        # We also set up cosmological information.  Set these to zero if
        # non-cosmological.
        #
        self.cosmological_simulation = 0
        self.current_redshift = 0
        self.omega_lambda = 0
        self.omega_matter = 0
        self.hubble_constant = 0

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        if isinstance(nib.load, NotAModule):
            return False
        try:
            nib.load(filename)
        except: (nib.ImageFileError)

        return False
