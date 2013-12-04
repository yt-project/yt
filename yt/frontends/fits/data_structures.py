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
    GridGeometryHandler
from yt.geometry.geometry_handler import \
    YTDataChunk
from yt.data_objects.static_output import \
    StaticOutput
from yt.utilities.definitions import \
    mpc_conversion, sec_conversion
from yt.utilities.io_handler import \
    io_registry
from yt.utilities.physical_constants import cm_per_mpc
from .fields import FITSFieldInfo, add_fits_field, KnownFITSFields
from yt.data_objects.field_info_container import FieldInfoContainer, NullFunc, \
     ValidateDataField, TranslationFunc

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
    
class FITSHierarchy(GridGeometryHandler):

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
        GridGeometryHandler.__init__(self,pf,data_style)

    def _initialize_data_storage(self):
        pass

    def _detect_fields(self):
        self.field_list = []
        for h in self._handle[self.parameter_file.first_image:]:
            if h.is_image:
                self.field_list.append(h.name.lower())
                        
    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        GridGeometryHandler._setup_classes(self, dd)
        self.object_types.sort()

    def _count_grids(self):
        self.num_grids = 1
                
    def _parse_hierarchy(self):
        f = self._handle # shortcut
        pf = self.parameter_file # shortcut
        
        # Initialize to the domain left / domain right
        self.grid_left_edge[0,:] = pf.domain_left_edge
        self.grid_right_edge[0,:] = pf.domain_right_edge
        self.grid_dimensions[0] = pf.domain_dimensions
        
        # This will become redundant, as _prepare_grid will reset it to its
        # current value.  Note that FLASH uses 1-based indexing for refinement
        # levels, but we do not, so we reduce the level by 1.
        self.grid_levels.flat[:] = 0
        self.grids = np.empty(self.num_grids, dtype='object')
        for i in xrange(self.num_grids):
            self.grids[i] = self.grid(i, self, self.grid_levels[i,0])
        
    def _populate_grid_objects(self):
        self.grids[0]._prepare_grid()
        self.grids[0]._setup_dx()
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

class FITSStaticOutput(StaticOutput):
    _hierarchy_class = FITSHierarchy
    _fieldinfo_fallback = FITSFieldInfo
    _fieldinfo_known = KnownFITSFields
    _handle = None
    
    def __init__(self, filename, data_style='fits',
                 primary_header = None,
                 sky_conversion = None,
                 storage_filename = None,
                 conversion_override = None):

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
        if conversion_override is None: conversion_override = {}
        self._conversion_override = conversion_override

        self.wcs = pywcs.WCS(self.primary_header)

        if self.wcs.wcs.cunit[0].name in ["deg","arcsec","arcmin","mas"]:
            self.sky_wcs = self.wcs.deepcopy()
            if sky_conversion is None:
                self._set_minimalist_wcs()
            else:
                dims = np.array(self.shape)
                ndims = len(self.shape)
                new_unit = sky_conversion[1]
                new_deltx = np.abs(self.wcs.wcs.cdelt[0])*sky_conversion[0]
                new_delty = np.abs(self.wcs.wcs.cdelt[1])*sky_conversion[0]
                self.wcs.wcs.cdelt = [new_deltx, new_delty]
                self.wcs.wcs.crpix = 0.5*(dims+1)
                self.wcs.wcs.crval = [0.0]*2
                self.wcs.wcs.cunit = [new_unit]*2
                self.wcs.wcs.ctype = ["LINEAR"]*2

        if not all(key in self.primary_header for key in
                   ["CRPIX1","CRVAL1","CDELT1","CUNIT1"]):
            self._set_minimalist_wcs()

        StaticOutput.__init__(self, fname, data_style)
        self.storage_filename = storage_filename
            
        self.refine_by = 2
        self._set_units()

    def _set_minimalist_wcs(self):
        mylog.warning("Could not determine WCS information. Using pixel units.")
        dims = np.array(self.shape)
        ndims = len(dims)
        self.wcs.wcs.crpix = 0.5*(dims+1)
        self.wcs.wcs.cdelt = [1.]*ndims
        self.wcs.wcs.crval = 0.5*(dims+1)
        self.wcs.wcs.cunit = ["pixel"]*ndims
        self.wcs.wcs.ctype = ["LINEAR"]*ndims

    def _set_units(self):
        """
        Generates the conversion to various physical _units based on the parameter file
        """
        self.units = {}
        self.time_units = {}
        if len(self.parameters) == 0:
            self._parse_parameter_file()
        self.conversion_factors = defaultdict(lambda: 1.0)
        file_unit = self.wcs.wcs.cunit[0].name.lower()
        if file_unit in mpc_conversion:
            self._setup_getunits_units()
        else:
            self._setup_nounits_units()
        self.parameters["Time"] = self.conversion_factors["Time"]
        self.time_units['1'] = 1
        self.units['1'] = 1.0
        self.units['unitary'] = 1.0 / (self.domain_right_edge - self.domain_left_edge).max()
        for unit in sec_conversion.keys():
            self.time_units[unit] = self.conversion_factors["Time"] / sec_conversion[unit]
        for p, v in self._conversion_override.items():
            self.conversion_factors[p] = v

    def _setup_comoving_units(self):
        pass

    def _setup_getunits_units(self):
        file_unit = self.wcs.wcs.cunit[0].name.lower()
        for unit in mpc_conversion.keys():
            self.units[unit] = mpc_conversion[unit]/mpc_conversion[file_unit]
        self.conversion_factors["Time"] = 1.0
                                            
    def _setup_nounits_units(self):
        for unit in mpc_conversion.keys():
            self.units[unit] = mpc_conversion[unit] / mpc_conversion["cm"]
        self.conversion_factors["Time"] = 1.0

    def _parse_parameter_file(self):
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])
        for k, v in self.primary_header.items():
            self.parameters[k] = v

        # Determine dimensionality

        self.dimensionality = self.primary_header["naxis"]
        self.geometry = "cartesian"

        self.domain_dimensions = np.array(self._handle[self.first_image].shape)
        if self.dimensionality == 2:
            self.domain_dimensions = np.append(self.domain_dimensions,
                                               [int(1)])
        ND = self.dimensionality
        
        le = [0.5]*ND
        re = [float(dim)+0.5 for dim in self.domain_dimensions]
        if ND == 2:
            xe, ye = self.wcs.wcs_pix2world([le[0],re[0]],
                                            [le[1],re[1]], 1)
            self.domain_left_edge = np.array([xe[0], ye[0], 0.0])
            self.domain_right_edge = np.array([xe[1], ye[1], 1.0]) 
        elif ND == 3:
            xe, ye, ze = world_edges = self.wcs.wcs_pix2world([le[0],re[0]],
                                                              [le[1],re[1]],
                                                              [le[2],re[2]], 1)
            self.domain_left_edge = np.array([xe[0], ye[0], ze[0]])
            self.domain_right_edge = np.array([xe[1], ye[1], ze[1]])

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


