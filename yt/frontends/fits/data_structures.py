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
import re

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
    mpc_conversion
from yt.utilities.io_handler import \
    io_registry
from .fields import FITSFieldInfo
from yt.utilities.decompose import \
    decompose_array, get_psize
from yt.units.unit_lookup_table import default_unit_symbol_lut
from yt.utilities.on_demand_imports import ap

known_units = dict([(unit.lower(),unit) for unit in default_unit_symbol_lut])
axes_prefixes = ["RA","DEC","V","ENER","FREQ"]

delimiters = ["*", "/", "-", "^"]
delimiters += [str(i) for i in xrange(10)]
regex_pattern = '|'.join(re.escape(_) for _ in delimiters)

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

    def _determine_image_units(self, fname, header):
        try:
            field_units = header["bunit"].lower().strip(" ")
            # FITS units always return upper-case, so we need to get
            # the right case by comparing against known units. This
            # only really works for common units.
            units = re.split(regex_pattern, field_units)
            for unit in units:
                if unit in known_units:
                    field_units = field_units.replace(unit, known_units[unit])
            self.parameter_file.field_units[fname] = field_units
        except:
            self.parameter_file.field_units[fname] = "dimensionless"

    def _ensure_same_dims(self, hdu):
        ds = self.parameter_file
        conditions = [hdu.header["naxis"] != ds.primary_header["naxis"]]
        for i in xrange(ds.naxis):
            nax = "naxis%d" % (i+1)
            conditions.append(hdu.header[nax] != ds.primary_header[nax])
        if np.any(conditions):
            return False
        else:
            return True

    def _detect_output_fields(self):
        ds = self.parameter_file
        self.field_list = []
        self._axis_map = {}
        self._file_map = {}
        self._ext_map = {}
        self._scale_map = {}
        # We create a field from each slice on the 4th axis
        if self.parameter_file.naxis == 4:
            naxis4 = self.parameter_file.primary_header["naxis4"]
        else:
            naxis4 = 1
        for i, fits_file in enumerate(self.parameter_file._fits_files):
            for j, hdu in enumerate(fits_file):
                if self._ensure_same_dims(hdu):
                    try:
                        fname = hdu.header["btype"].lower()
                    except:
                        fname = hdu.name.lower()
                    for k in xrange(naxis4):
                        if naxis4 > 1:
                            fname += "_%s_%d" % (hdu.header["CTYPE4"], k+1)
                        if self.pf.num_files > 1:
                            try:
                                fname += "_%5.3fGHz" % (hdu.header["restfreq"]/1.0e9)
                            except:
                                fname += "_%5.3fGHz" % (hdu.header["restfrq"]/1.0e9)
                            else:
                                fname += "_field_%d" % (i)
                        self._axis_map[fname] = k
                        self._file_map[fname] = fits_file
                        self._ext_map[fname] = j
                        self._scale_map[fname] = [0.0,1.0]
                        if "bzero" in hdu.header:
                            self._scale_map[fname][0] = hdu.header["bzero"]
                        if "bscale" in hdu.header:
                            self._scale_map[fname][1] = hdu.header["bscale"]
                        self.field_list.append((self.dataset_type, fname))
                        mylog.info("Adding field %s to the list of fields." % (fname))
                        self._determine_image_units(fname, hdu.header)
                else:
                    mylog.warning("Image block %s does not have " % (hdu.name.lower()) +
                                  "the same dimensions as the primary and will not be " +
                                  "available as a field.")


        # For line fields, we still read the primary field. Not sure how to extend this
        # For now, we pick off the first field from the field list.
        line_db = self.parameter_file.line_database
        primary_fname = self.field_list[0][1]
        for k, v in line_db.iteritems():
            mylog.info("Adding line field: %s at offset %i" % (k, v))
            self.field_list.append((self.dataset_type, k))
            self._ext_map[k] = self._ext_map[primary_fname]
            self._axis_map[k] = self._axis_map[primary_fname]
            self._file_map[k] = self._file_map[primary_fname]
            self.parameter_file.field_units[k] = self.parameter_file.field_units[primary_fname]

    def _count_grids(self):
        self.num_grids = self.pf.nprocs

    def _parse_index(self):
        f = self._handle # shortcut
        pf = self.parameter_file # shortcut

        # If nprocs > 1, decompose the domain into virtual grids
        if pf.nprocs > 1:
            bbox = np.array([[le,re] for le, re in zip(pf.domain_left_edge,
                                                       pf.domain_right_edge)])
            psize = get_psize(np.array(pf.domain_dimensions), pf.nprocs)
            gle, gre, shapes, slices = decompose_array(pf.domain_dimensions, psize, bbox)
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

    def __init__(self, filename,
                 dataset_type = 'fits',
                 slave_files = [],
                 nprocs = None,
                 storage_filename = None,
                 nan_mask = None,
                 folded_axis = None,
                 folded_width = None,
                 line_database = None,
                 suppress_astropy_warnings = True):
        self.folded_axis = folded_axis
        self.folded_width = folded_width
        self._unfolded_domain_dimensions = None
        if line_database is None:
            line_database = {}
        self.line_database = line_database

        if suppress_astropy_warnings:
            warnings.filterwarnings('ignore', module="astropy", append=True)
        self.filenames = [filename] + slave_files
        self.num_files = len(self.filenames)
        self.fluid_types += ("fits",)
        if nan_mask is None:
            self.nan_mask = {}
        elif isinstance(nan_mask, float):
            self.nan_mask = {"all":nan_mask}
        elif isinstance(nan_mask, dict):
            self.nan_mask = nan_mask
        self.nprocs = nprocs
        self._handle = ap.pyfits.open(self.filenames[0],
                                      memmap=True,
                                      do_not_scale_image_data=True,
                                      ignore_blank=True)
        self._fits_files = [self._handle]
        if self.num_files > 1:
            for fits_file in slave_files:
                self._fits_files.append(ap.pyfits.open(fits_file,
                                                       memmap=True,
                                                       do_not_scale_image_data=True,
                                                       ignore_blank=True))

        self.first_image = 0 # Assumed for now
        self.primary_header = self._handle[self.first_image].header
        self.wcs = ap.pywcs.WCS(header=self.primary_header)
        self.axis_names = {}
        self.naxis = self.primary_header["naxis"]
        for i, ax in enumerate("xyz"[:self.naxis]):
            self.axis_names[self.primary_header["ctype%d" % (i+1)]] = ax
        self.refine_by = 2

        Dataset.__init__(self, filename, dataset_type)
        self.storage_filename = storage_filename

    def _set_code_unit_attributes(self):
        """
        Generates the conversion to various physical _units based on the parameter file
        """
        file_unit = None
        for i, unit in enumerate(self.wcs.wcs.cunit):
            if unit in mpc_conversion.keys():
                file_unit = unit.name
                idx = i
                break
        if file_unit is None:
            self.no_cgs_equiv_length = True
            mylog.warning("No length conversion provided. Assuming 1 = 1 cm.")
            length_factor = 1.0
            length_unit = "cm"
        else:
            length_factor = self.wcs.wcs.cdelt[idx]
            length_unit = str(file_unit)
        self.length_unit = self.quan(length_factor,length_unit)
        self.mass_unit = self.quan(1.0, "g")
        self.time_unit = self.quan(1.0, "s")
        self.velocity_unit = self.quan(1.0, "cm/s")

    def _parse_parameter_file(self):
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])

        # Determine dimensionality

        self.dimensionality = self.naxis
        self.geometry = "cartesian"

        # Sometimes a FITS file has a 4D datacube, in which case
        # we take the 4th axis and assume it consists of different fields.
        if self.dimensionality == 4: self.dimensionality = 3

        dims = [self.primary_header["naxis%d" % (i+1)] for i in xrange(self.naxis)]
        self.domain_dimensions = np.array(dims)[:self.dimensionality]
        if self.dimensionality == 2:
            self.domain_dimensions = np.append(self.domain_dimensions,
                                               [int(1)])

        self.domain_left_edge = np.array([0.5]*3)
        self.domain_right_edge = np.array([float(dim)+0.5 for dim in self.domain_dimensions])

        if self.folded_axis is not None:
            self.domain_left_edge[self.folded_axis] = -self.folded_width/2.
            self.domain_right_edge[self.folded_axis] = self.folded_width/2.
            self._unfolded_domain_dimensions = self.domain_dimensions.copy()
            self.domain_dimensions[self.folded_axis] = int(self.folded_width)

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

        # If nprocs is None, do some automatic decomposition of the domain
        if self.nprocs is None:
            self.nprocs = np.around(np.prod(self.domain_dimensions) /
                                    32**self.dimensionality).astype("int")
            self.nprocs = min(self.nprocs, 512)

        # Check to see if this data is in (RA,Dec,?) format
        self.ppv_data = False
        x = np.zeros((self.dimensionality), dtype="bool")
        for ap in axes_prefixes:
            x += np_char.startswith(self.axis_names.keys()[:self.dimensionality], ap)
        if x.sum() == self.dimensionality: self._setup_ppv()

    def _setup_ppv(self):

        self.ppv_data = True

        end = min(self.dimensionality+1,4)
        ctypes = np.array([self.primary_header["CTYPE%d" % (i)] for i in xrange(1,end)])
        self.ra_axis = np.where(np_char.startswith(ctypes, "RA"))[0][0]
        self.dec_axis = np.where(np_char.startswith(ctypes, "DEC"))[0][0]

        if self.wcs.naxis > 2:

            self.vel_axis = np_char.startswith(ctypes, "V")
            self.vel_axis += np_char.startswith(ctypes, "FREQ")
            self.vel_axis += np_char.startswith(ctypes, "ENER")
            self.vel_axis = np.where(self.vel_axis)[0][0]
            self.vel_name = ctypes[self.vel_axis].lower()

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

        else:

            self.wcs_2d = self.wcs
            self.wcs_1d = None
            self.vel_axis = 2
            self.vel_name = "z"

    def __del__(self):
        for file in self._fits_files:
            file.close()
        self._handle.close()

    @classmethod
    def _is_valid(cls, *args, **kwargs):
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
            valid = fileh[0].header["naxis"] >= 2
            fileh.close()
            return valid
        except:
            pass
        return False
