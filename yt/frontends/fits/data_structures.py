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

import stat
import types
import numpy as np
import numpy.core.defchararray as np_char
import weakref
import warnings

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
from .fields import FITSFieldInfo, FITSXYVFieldInfo
from yt.utilities.decompose import \
    decompose_array, get_psize, decompose_array_nocopy

class astropy_imports:
    _pyfits = None
    @property
    def pyfits(self):
        if self._pyfits is None:
            try:
                import astropy.io.fits as pyfits
                self.log
            except ImportError:
                pyfits = None
            self._pyfits = pyfits
        return self._pyfits

    _pywcs = None
    @property
    def pywcs(self):
        if self._pywcs is None:
            try:
                import astropy.wcs as pywcs
                self.log
            except ImportError:
                pywcs = None
            self._pywcs = pywcs
        return self._pywcs

    _log = None
    @property
    def log(self):
        if self._log is None:
            try:
                from astropy import log
                if log.exception_logging_enabled():
                    log.disable_exception_logging()
            except ImportError:
                log = None
            self._log = log
        return self._log

ap = astropy_imports()

angle_units = ["deg","arcsec","arcmin","mas"]
all_units = angle_units + mpc_conversion.keys()

known_units = {"k":"K",
               "jy":"Jy"}

def fits_file_validator(ds, *args, **kwargs):
    ext = args[0].rsplit(".", 1)[-1]
    if ext.upper() == "GZ":
        # We don't know for sure that there will be > 1
        ext = args[0].rsplit(".", 1)[0].rsplit(".", 1)[-1]
    if ext.upper() not in ("FITS", "FTS"):
        return False
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=UserWarning, append=True)
            fileh = ap.pyfits.open(args[0])
        if ds._check_axes(fileh): return True
    except:
        pass
    return False

class FITSGrid(AMRGridPatch):
    _id_offset = 0
    def __init__(self, id, index, level):
        AMRGridPatch.__init__(self, id, filename = index.index_filename,
                              index = index)
        self.Parent = None
        self.Children = []
        self.Level = 0

    def __repr__(self):
        return "FITSGrid_%04i (%s)" % (self.id, self.ActiveDimensions)

class FITSHierarchy(GridIndex):

    grid = FITSGrid

    def __init__(self,pf,dataset_type='fits'):
        self.dataset_type = dataset_type
        self.field_indexes = {}
        self.parameter_file = weakref.proxy(pf)
        # for now, the index file is the parameter file!
        self.index_filename = self.parameter_file.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        self._handle = pf._handle
        self.float_type = np.float64
        GridIndex.__init__(self,pf,dataset_type)

    def _initialize_data_storage(self):
        pass

    def _detect_image_units(self, fname, header):
        try:
            field_units = header["bunit"].lower().strip(" ")
            # FITS units always return upper-case, so we need to get
            # the right case by comparing against known units
            for name in known_units:
                if field_units.find(name) > -1:
                    field_units = field_units.replace(name, known_units[name])
            self.parameter_file.field_units[fname] = field_units
        except:
            pass

    def _detect_output_fields(self):
        self.field_list = []
        self._field_map = {}
        for h in self._handle[self.parameter_file.first_image:]:
            if h.header["naxis"] >= 2:
                if self.parameter_file.four_dims:
                    for idx in range(h.header["naxis4"]):
                        fname = h.name.lower()+"_%d" % (idx)
                        self._field_map[fname] = idx
                        self.field_list.append((self.dataset_type, fname))
                        self._detect_image_units(fname, h.header)
                else:
                    fname = h.name.lower()
                    self.field_list.append((self.dataset_type, fname))
                    self._detect_image_units(fname, h.header)
        line_db = self.parameter_file.line_database
        for k, v in line_db.iteritems():
            print "Adding line: ", k, v
            self.field_list.append((self.dataset_type, k))
            self._field_map[k] = 0
            self.parameter_file.field_units[k] = "Jy/beam"

    def _count_grids(self):
        self.num_grids = self.pf.nprocs

    def _parse_index(self):
        f = self._handle # shortcut
        pf = self.parameter_file # shortcut

        if pf.nprocs > 1:
            bbox = np.array([[le,re] for le, re in zip(pf.domain_left_edge,
                                                       pf.domain_right_edge)])
            psize = get_psize(np.array(pf.domain_dimensions), pf.nprocs)
            gle, gre, shapes = decompose_array_nocopy(pf.domain_dimensions, psize, bbox)
            self.grid_left_edge = self.pf.arr(gle, "code_length")
            self.grid_right_edge = self.pf.arr(gre, "code_length")
            self.grid_dimensions = np.array([shape for shape in shapes], dtype="int32")
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
        self.io = io_registry[self.dataset_type](self.parameter_file)

    def _chunk_io(self, dobj, cache = True, local_only = False):
        # local_only is only useful for inline datasets and requires
        # implementation by subclasses.
        gfiles = defaultdict(list)
        gobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for g in gobjs:
            gfiles[g.id].append(g)
        for fn in sorted(gfiles):
            gs = gfiles[fn]
            yield YTDataChunk(dobj, "io", gs, self._count_selection(dobj, gs),
                              cache = cache)

class FITSDataset(Dataset):
    _index_class = FITSHierarchy
    _field_info_class = FITSFieldInfo
    _dataset_type = "fits"
    _handle = None

    def __init__(self, filename, dataset_type='fits',
                 nprocs = None,
                 storage_filename = None,
                 mask_nans = False,
                 folded_axis=None,
                 folded_width=None,
                 line_database=None
                 ):
        self.folded_axis = folded_axis
        self.folded_width = folded_width
        if line_database is None:
            line_database = {}
        self.line_database = line_database
        self.fluid_types += ("fits",)
        self.mask_nans = mask_nans
        self.nprocs = nprocs
        self._handle = ap.pyfits.open(filename, memmap=True, do_not_scale_image_data=True)
        for i, h in enumerate(self._handle):
            if h.header["naxis"] >= 2:
                self.first_image = i
                break

        self.primary_header = self._handle[self.first_image].header
        self.shape = self._handle[self.first_image].shape
        self.wcs = ap.pywcs.WCS(header=self.primary_header)

        self.file_unit = None
        for i, unit in enumerate(self.wcs.wcs.cunit):
            if unit in all_units:
                self.file_unit = unit.name
                idx = i
                break
        self.new_unit = None
        self.pixel_scale = 1.0
        if self.file_unit in mpc_conversion:
            self.new_unit = self.file_unit
            self.pixel_scale = self.wcs.wcs.cdelt[idx]

        self.refine_by = 2
        self.four_dims = False

        Dataset.__init__(self, filename, dataset_type)
        self.storage_filename = storage_filename

    def _set_code_unit_attributes(self):
        """
        Generates the conversion to various physical _units based on the parameter file
        """
        if self.new_unit is not None:
            length_factor = self.pixel_scale
            length_unit = str(self.new_unit)
        else:
            mylog.warning("No length conversion provided. Assuming 1 = 1 cm.")
            length_factor = 1.0
            length_unit = "cm"
        self.length_unit = self.quan(length_factor,length_unit)
        self.mass_unit = self.quan(1.0, "g")
        self.time_unit = self.quan(1.0, "s")
        self.velocity_unit = self.quan(1.0, "cm/s")

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

        if self.nprocs is None:
            self.nprocs = np.around(np.prod(self.domain_dimensions) /
                                    32**self.dimensionality).astype("int")
            self.nprocs = min(self.nprocs, 2500)

    def __del__(self):
        self._handle.close()

    @classmethod
    def _check_axes(cls, handle):
        for h in handle:
            if h.header["naxis"] >= 2:
                axes_names = [h.header["CTYPE%d" % (ax)] for ax in xrange(1,4)]
                a = np_char.startswith(axes_names, "RA")
                b = np_char.startswith(axes_names, "DEC")
                c = np_char.startswith(axes_names, "VEL")
                d = np_char.startswith(axes_names, "FREQ")
                e = np_char.startswith(axes_names, "ENER")
                if (a+b+c+d+e).sum() != 3:
                    handle.close()
                    return True
        return False

    @classmethod
    def _is_valid(cls, *args, **kwargs):
        return fits_file_validator(cls, *args, **kwargs)

class FITSXYVDataset(FITSDataset):
    _dataset_type = "xyv_fits"
    _field_info_class = FITSXYVFieldInfo

    def __init__(self, filename,
                 dataset_type='xyv_fits',
                 nprocs=None,
                 storage_filename=None,
                 mask_nans=False,
                 folded_axis=None,
                 folded_width=None,
                 line_database=None
                 ):

        self.fluid_types += ("xyv_fits",)

        super(FITSXYVDataset, self).__init__(filename, dataset_type=dataset_type,
                                             nprocs=nprocs,
                                             storage_filename=storage_filename,
                                             mask_nans=mask_nans,
                                             folded_axis=folded_axis,
                                             folded_width=folded_width,
                                             line_database=line_database)
        self.axes_names = [self.primary_header["CTYPE%d" % (ax)] for ax in xrange(1,4)]
        self.ra_axis = np.where(np_char.startswith(self.axes_names, "RA"))[0][0]
        self.dec_axis = np.where(np_char.startswith(self.axes_names, "DEC"))[0][0]
        self.vel_axis = np_char.startswith(self.axes_names, "VEL")
        self.vel_axis += np_char.startswith(self.axes_names, "FREQ")
        self.vel_axis += np_char.startswith(self.axes_names, "ENER")
        self.vel_axis = np.where(self.vel_axis)[0][0]

        self.wcs_2d = ap.pywcs.WCS(naxis=2)
        self.wcs_2d.wcs.crpix = self.wcs.wcs.crpix[[self.ra_axis, self.dec_axis]]
        self.wcs_2d.wcs.cdelt = self.wcs.wcs.cdelt[[self.ra_axis, self.dec_axis]]
        self.wcs_2d.wcs.crval = self.wcs.wcs.crval[[self.ra_axis, self.dec_axis]]
        self.wcs_2d.wcs.cunit = [str(self.wcs.wcs.cunit[self.ra_axis]),
                                 str(self.wcs.wcs.cunit[self.dec_axis])]
        self.wcs_2d.wcs.ctype = [self.wcs.wcs.ctype[self.ra_axis],
                                 self.wcs.wcs.ctype[self.dec_axis]]

        self.wcs_1d = ap.pywcs.WCS(naxis=1)
        self.wcs_1d.wcs.crpix = [self.wcs.wcs.crpix[self.vel_axis]]
        self.wcs_1d.wcs.cdelt = [self.wcs.wcs.cdelt[self.vel_axis]]
        self.wcs_1d.wcs.crval = [self.wcs.wcs.crval[self.vel_axis]]
        self.wcs_1d.wcs.cunit = [str(self.wcs.wcs.cunit[self.vel_axis])]
        self.wcs_1d.wcs.ctype = [self.wcs.wcs.ctype[self.vel_axis]]

    def _parse_parameter_file(self):

        super(FITSXYVDataset, self)._parse_parameter_file()

        if self.dimensionality == 4:
            self.dimensionality = 3
            self.four_dims = True
            self.domain_dimensions = self.domain_dimensions[:3]
            self.domain_left_edge = self.domain_left_edge[:3]
            self.domain_right_edge = self.domain_right_edge[:3]

        if self.folded_axis is not None:
            ax = self.folded_axis
            ratio = self.folded_width/self.domain_dimensions[ax]
            self.domain_dimensions[ax] = int(self.folded_width)
            self.domain_left_edge[ax] = -self.folded_width/2.
            self.domain_right_edge[ax] = self.folded_width/2.

        if self.nprocs is None:
            self.nprocs = np.around(np.prod(self.domain_dimensions) /
                                    32**self.dimensionality).astype("int")
            self.nprocs = min(self.nprocs, 2500)

    @classmethod
    def _check_axes(cls, handle):
        for h in handle:
            if h.header["naxis"] >= 3:
                axes_names = [h.header["CTYPE%d" % (ax)] for ax in xrange(1,4)]
                a = np_char.startswith(axes_names, "RA")
                b = np_char.startswith(axes_names, "DEC")
                c = np_char.startswith(axes_names, "VEL")
                d = np_char.startswith(axes_names, "FREQ")
                e = np_char.startswith(axes_names, "ENER")
                if (a+b+c+d+e).sum() == 3:
                    handle.close()
                    return True
        return False

    @classmethod
    def _is_valid(cls, *args, **kwargs):
        return fits_file_validator(cls, *args, **kwargs)
