"""
ExodusII data structures



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.funcs import *
from yt.config import ytcfg
from yt.data_objects.data_containers import \
    YTFieldData, \
    YTDataContainer, \
    YTSelectionContainer
from yt.data_objects.grid_patch import \
    AMRGridPatch
from yt.geometry.grid_geometry_handler import \
    GridIndex
from yt.data_objects.static_output import \
    Dataset
from yt.geometry.unstructured_mesh_handler import \
           UnstructuredIndex
from yt.data_objects.unstructured_mesh import \
           SemiStructuredMesh, \
           UnstructuredMesh
from .fields import ExodusIIFieldInfo
from moose_parser.data import Data

class ExodusIIHandler(object):
    def __init__(self, left_edges, right_edges):
        self.left_edges = np.array(left_edges)
        self.right_edges = np.array(right_edges)
            
    def get_fields(self):
        return self.fields.all_fields

    def get_particle_type(self, field) :

        if field in self.particle_types :
            return self.particle_types[field]
        else :
            return False

class ExodusIIGrid(AMRGridPatch):
    _id_offset = 0
    def __init__(self, id, index, level, start, dimensions):
        AMRGridPatch.__init__(self, id, filename=index.index_filename,
                              index=index)
        self.Parent = []
        self.Children = []
        self.Level = level
        self.start_index = start.copy()
        self.stop_index = self.start_index + dimensions
        self.ActiveDimensions = dimensions.copy()

    def __repr__(self):
        return "ExodusIIGrid_%04i (%s)" % (self.id, self.ActiveDimensions)

class ExodusIIHierarchy(GridIndex):
    grid = ExodusIIGrid

    def __init__(self, ds, dataset_type='exodusii'):
        self.dataset_type = dataset_type
        # for now, the index file is the dataset!
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        GridIndex.__init__(self, ds, dataset_type)

    def _detect_output_fields(self):
        # This needs to set a self.field_list that contains all the available,
        # on-disk fields.
        # NOTE: Each should be a tuple, where the first element is the on-disk
        # fluid type or particle type.  Convention suggests that the on-disk
        # fluid type is usually the dataset_type and the on-disk particle type
        # (for a single population of particles) is "io".
        pass

    def _count_grids(self):
        # This needs to set self.num_grids
        pass

    def _parse_index(self):
        # This needs to fill the following arrays, where N is self.num_grids:
        #   self.grid_left_edge         (N, 3) <= float64
        #   self.grid_right_edge        (N, 3) <= float64
        #   self.grid_dimensions        (N, 3) <= int
        #   self.grid_particle_count    (N, 1) <= int
        #   self.grid_levels            (N, 1) <= int
        #   self.grids                  (N, 1) <= grid objects
        #
        pass

    def _populate_grid_objects(self):
        # For each grid, this must call:
        #   grid._prepare_grid()
        #   grid._setup_dx()
        # This must also set:
        #   grid.Children <= list of child grids
        #   grid.Parent   <= parent grid
        # This is handled by the frontend because often the children must be
        # identified.
        pass

class ExodusIIDataset(Dataset):
    _index_class = ExodusIIHierarchy
    _field_info_class = ExodusIIFieldInfo
    _dataset_type = 'exodusii'

    def __init__(self, exodusii_handler, storage_filename=None, geometry="cartesian"):
        self.geometry = geometry
        self.exodusii_handler = exodusii_handler
        name = "InMemoryParameterFile_%s" % (uuid.uuid4().hex)
        from yt.data_objects.static_output import _cached_datasets
        _cached_datasets[name] = self
        Dataset.__init__(self, name, self._dataset_type)

    def _parse_parameter_file(self):
        # This needs to set up the following items.  Note that these are all
        # assumed to be in code units; domain_left_edge and domain_right_edge
        # will be updated to be in code units at a later time.  This includes
        # the cosmological parameters.
        #
        #   self.unique_identifier
        #   self.parameters             <= full of code-specific items of use
        #   self.domain_left_edge       <= array of float64
        #   self.domain_right_edge      <= array of float64
        #   self.dimensionality         <= int
        #   self.domain_dimensions      <= array of int64
        #   self.periodicity            <= three-element tuple of booleans
        #   self.current_time           <= simulation time in code units
        #
        # We also set up cosmological information.  Set these to zero if
        # non-cosmological.
        #
        #   self.cosmological_simulation    <= int, 0 or 1
        #   self.current_redshift           <= float
        #   self.omega_lambda               <= float
        #   self.omega_matter               <= float
        #   self.hubble_constant            <= float
        pass

    def _set_units(self):
        self.field_units = self.exodusii_handler.field_units

    def _set_code_unit_attributes(self):
        base_units = self.exodusii_handler.code_units
        attrs = ('length_unit', 'mass_unit', 'time_unit', 'velocity_unit', 'magnetic_unit')
        cgs_units = ('cm', 'g', 's', 'cm/s', 'gauss')
        for unit, attr, cgs_unit in zip(base_units, attrs, cgs_units):
            if isinstance(unit, string_types):
                uq = self.quan(1.0, unit)
            elif isinstance(unit, numeric_type):
                uq = self.quan(unit, cgs_unit)
            elif isinstance(unit, YTQuantity):
                uq = unit
            elif isinstance(unit, tuple):
                uq = self.quan(unit[0], unit[1])
            else:
                raise RuntimeError("%s (%s) is invalid." % (attr, unit))
            setattr(self, attr, uq)


    @classmethod
    def _is_valid(self, *args, **kwargs):
        # This accepts a filename or a set of arguments and returns True or
        # False depending on if the file is of the type requested.
        return False

class ExodusIIDictFieldHandler(dict):
    _additional_fields = ()

    @property
    def all_fields(self):
        self_fields = chain.from_iterable(s.keys() for s in self.values())
        self_fields = list(set(self_fields))
        fields = list(self._additional_fields) + self_fields
        fields = list(set(fields))
        return fields    

class ExodusIIUnstructuredMesh(UnstructuredMesh):
    _index_offset = 1

    def __init__(self, *args, **kwargs):
        super(ExodusIIUnstructuredMesh, self).__init__(*args, **kwargs)

class ExodusIIUnstructuredIndex(UnstructuredIndex):

    def __init__(self, ds, dataset_type = None):
        self.exodusii_handler = ds.exodusii_handler
        super(ExodusIIUnstructuredIndex, self).__init__(ds, dataset_type)

    def _initialize_mesh(self):
        coords = ensure_list(self.stream_handler.fields.pop("coordinates"))
        connec = ensure_list(self.stream_handler.fields.pop("connectivity"))
        self.meshes = [StreamUnstructuredMesh(
          i, self.index_filename, c1, c2, self)
          for i, (c1, c2) in enumerate(zip(connec, coords))]

    def _setup_data_io(self):
        if self.stream_handler.io is not None:
            self.io = self.stream_handler.io
        else:
            self.io = io_registry[self.dataset_type](self.ds)

    def _detect_output_fields(self):
        self.field_list = list(set(self.stream_handler.get_fields()))

class ExodusIIUnstructuredMeshDataset(ExodusIIDataset):
    _index_class = ExodusIIUnstructuredIndex
    _field_info_class = ExodusIIFieldInfo
    _dataset_type = "exodusii_unstructured"

def load_unstructured_mesh(data, connectivity, coordinates,
                         length_unit = None, bbox=None, sim_time=0.0,
                         mass_unit = None, time_unit = None,
                         velocity_unit = None, magnetic_unit = None,
                         periodicity=(False, False, False),
                         geometry = "cartesian"):
    
    handler = ExodusIIHandler()
    eds = ExodusIIUnstructuredMeshDataset(handler, geometry = geometry)

    return eds
    
