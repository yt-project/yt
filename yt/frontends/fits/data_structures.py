"""
FITS-specific data structures
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

try:
    import astropy.io.fits as pyfits
    import astropy.wcs as pywcs
except ImportError:
    pass

import stat
import numpy as np
import weakref

from yt.config import ytcfg
from yt.funcs import *
from yt.data_objects.grid_patch import \
    AMRGridPatch
from yt.geometry.grid_geometry_handler import \
    GridIndex
from yt.geometry.geometry_handler import \
    YTDataChunk
from yt.data_objects.static_output import \
    Dataset
from yt.utilities.definitions import \
    mpc_conversion, sec_conversion
from yt.utilities.io_handler import \
    io_registry
from yt.utilities.physical_constants import cm_per_mpc
from .fields import FITSFieldInfo
from yt.utilities.decompose import \
    decompose_array, get_psize

angle_units = ["deg","arcsec","arcmin","mas"]
all_units = angle_units + mpc_conversion.keys()

class FITSGrid(AMRGridPatch):
    _id_offset = 0
    def __init__(self, id, hierarchy, level):
        AMRGridPatch.__init__(self, id, filename = hierarchy.hierarchy_filename,
                              hierarchy = hierarchy)
        self.Parent = None
        self.Children = []
        self.Level = 0

    def __repr__(self):
        return "FITSGrid_%04i (%s)" % (self.id, self.ActiveDimensions)
    
class FITSHierarchy(GridIndex):

    grid = FITSGrid
    
    def __init__(self,pf,data_style='fits'):
        self.data_style = data_style
        self.field_indexes = {}
        self.parameter_file = weakref.proxy(pf)
        # for now, the hierarchy file is the parameter file!
        self.hierarchy_filename = self.parameter_file.parameter_filename
        self.directory = os.path.dirname(self.hierarchy_filename)
        self._handle = pf._handle
        self.float_type = np.float64
        GridIndex.__init__(self,pf,data_style)

    def _initialize_data_storage(self):
        pass

    def _detect_output_fields(self):
        self.field_list = []
        for h in self._handle[self.parameter_file.first_image:]:
            if h.is_image:
                self.field_list.append(("fits", h.name.lower()))
                        
    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        GridIndex._setup_classes(self, dd)
        self.object_types.sort()

    def _count_grids(self):
        self.num_grids = self.pf.nprocs
                
    def _parse_hierarchy(self):
        f = self._handle # shortcut
        pf = self.parameter_file # shortcut

        if pf.nprocs > 1:
            bbox = np.array([[le,re] for le, re in zip(pf.domain_left_edge,
                                                       pf.domain_right_edge)])
            psize = get_psize(np.array(pf.domain_dimensions), pf.nprocs)
            temp_arr = np.zeros(pf.domain_dimensions)
            gle, gre, temp_arr = decompose_array(temp_arr, psize, bbox)
            self.grid_left_edge = self.pf.arr(gle, "code_length")
            self.grid_right_edge = self.pf.arr(gre, "code_length")
            self.grid_dimensions = np.array([grid.shape for grid in temp_arr], dtype="int32")
            del temp_arr
        else:
            self.grid_left_edge[0,:] = pf.domain_left_edge
            self.grid_right_edge[0,:] = pf.domain_right_edge
            self.grid_dimensions[0] = pf.domain_dimensions
        
        self.grid_levels.flat[:] = 0
        self.grids = np.empty(self.num_grids, dtype='object')
        for i in xrange(self.num_grids):
            self.grids[i] = self.grid(i, self, self.grid_levels[i,0])
        
    def _populate_grid_objects(self):
        for i in xrange(self.num_grids):
            self.grids[i]._prepare_grid()
            self.grids[i]._setup_dx()
        self.max_level = 0

    def _setup_derived_fields(self):
        super(FITSHierarchy, self)._setup_derived_fields()
        [self.parameter_file.conversion_factors[field] 
         for field in self.field_list]
        for field in self.field_list:
            if field not in self.derived_field_list:
                self.derived_field_list.append(field)

        for field in self.derived_field_list:
            f = self.parameter_file.field_info[field]
            if f._function.func_name == "_TranslationFunc":
                # Translating an already-converted field
                self.parameter_file.conversion_factors[field] = 1.0 
                
    def _setup_data_io(self):
        self.io = io_registry[self.data_style](self.parameter_file)

class FITSDataset(Dataset):
    _index_class = FITSHierarchy
    _field_info_class = FITSFieldInfo
    _data_style = "fits"
    _handle = None
    
    def __init__(self, filename, data_style='fits',
                 primary_header = None,
                 sky_conversion = None,
                 storage_filename = None,
                 mask_nans = True,
                 nprocs=1):
        self.fluid_types += ("fits",)
        self.mask_nans = mask_nans
        self.nprocs = nprocs
        if isinstance(filename, pyfits.HDUList):
            self._handle = filename
            fname = filename.filename()
        else:
            self._handle = pyfits.open(filename)
            fname = filename
        for i, h in enumerate(self._handle):
            if h.is_image and h.data is not None:
                self.first_image = i
                break
        
        if primary_header is None:
            self.primary_header = self._handle[self.first_image].header
        else:
            self.primary_header = primary_header
        self.shape = self._handle[self.first_image].shape

        self.wcs = pywcs.WCS(header=self.primary_header)

        self.file_unit = None
        for i, unit in enumerate(self.wcs.wcs.cunit):
            if unit in all_units:
                self.file_unit = unit.name
                idx = i
                break
        self.new_unit = None
        self.pixel_scale = 1.0
        if self.file_unit in angle_units:
            if sky_conversion is not None:
                self.new_unit = sky_conversion[1]
                self.pixel_scale = np.abs(self.wcs.wcs.cdelt[idx])*sky_conversion[0]
        elif self.file_unit in mpc_conversion:
            self.new_unit = self.file_unit
            self.pixel_scale = self.wcs.wcs.cdelt[idx]

        Dataset.__init__(self, fname, data_style)
        self.storage_filename = storage_filename
            
        self.refine_by = 2

    def _set_code_unit_attributes(self):
        """
        Generates the conversion to various physical _units based on the parameter file
        """
        if self.new_unit is not None:
            length_factor = self.pixel_scale
            length_unit = self.new_unit
        else:
            mylog.warning("No length conversion provided. Assuming 1 = 1 cm.")
            length_factor = 1.0
            length_unit = "cm"
        self.length_unit = self.quan(length_factor,length_unit).in_cgs()
        self.mass_unit = self.quan(1.0, "g")
        self.time_unit = self.quan(1.0, "s")
        self.velocity_unit = self.quan(1.0, "cm/s")        
        DW = self.domain_right_edge-self.domain_left_edge
        self.unit_registry.modify("unitary", DW[:self.dimensionality].max())

    def _parse_parameter_file(self):
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])
        for k, v in self.primary_header.items():
            self.parameters[k] = v

        # Determine dimensionality

        self.dimensionality = self.primary_header["naxis"]
        self.geometry = "cartesian"

        dims = self._handle[self.first_image].shape[::-1]
        self.domain_dimensions = np.array(dims)
        if self.dimensionality == 2:
            self.domain_dimensions = np.append(self.domain_dimensions,
                                               [int(1)])
            
        self.domain_left_edge = np.array([0.5]*3)
        self.domain_right_edge = np.array([float(dim)+0.5 for dim in self.domain_dimensions])

        if self.dimensionality == 2:
            self.domain_left_edge[-1] = 0.5
            self.domain_right_edge[-1] = 1.5
            
        # Get the simulation time
        try:
            self.current_time = self.parameters["time"]
        except:
            mylog.warning("Cannot find time")
            self.current_time = 0.0
            pass
        
        # For now we'll ignore these
        self.periodicity = (False,)*3
        self.current_redshift = self.omega_lambda = self.omega_matter = \
            self.hubble_constant = self.cosmological_simulation = 0.0

    def __del__(self):
        self._handle.close()

    @classmethod
    def _is_valid(self, *args, **kwargs):
        try:
            if isinstance(args[0], pyfits.HDUList):
                for h in args[0]:
                    if h.is_image and h.data is not None:
                        return True
        except:
            pass
        try:
            fileh = pyfits.open(args[0])
            for h in fileh:
                if h.is_image and h.data is not None:
                    fileh.close()
                    return True
            fileh.close()
        except:
            pass
        return False


