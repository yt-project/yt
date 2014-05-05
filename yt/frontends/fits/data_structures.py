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
from yt.utilities.io_handler import \
    io_registry
from .fields import FITSFieldInfo
from yt.utilities.decompose import \
    decompose_array, get_psize
from yt.units.unit_lookup_table import \
    default_unit_symbol_lut, \
    prefixable_units, \
    unit_prefixes
from yt.units import dimensions
from yt.units.yt_array import YTQuantity
from yt.utilities.on_demand_imports import _astropy

lon_prefixes = ["X","RA","GLON"]
lat_prefixes = ["Y","DEC","GLAT"]
vel_prefixes = ["V","ENER","FREQ","WAV"]
delimiters = ["*", "/", "-", "^"]
delimiters += [str(i) for i in xrange(10)]
regex_pattern = '|'.join(re.escape(_) for _ in delimiters)

field_from_unit = {"Jy":"intensity",
                   "K":"temperature"}

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

    def _guess_name_from_units(self, units):
        for k,v in field_from_unit.items():
            if k in units:
                mylog.warning("Guessing this is a %s field based on its units of %s." % (v,k))
                return v
        return None

    def _determine_image_units(self, header, known_units):
        try:
            field_units = header["bunit"].lower().strip(" ").replace(" ", "")
            # FITS units always return upper-case, so we need to get
            # the right case by comparing against known units. This
            # only really works for common units.
            units = set(re.split(regex_pattern, field_units))
            units.remove('')
            n = int(0)
            for unit in units:
                if unit in known_units:
                    field_units = field_units.replace(unit, known_units[unit])
                    n += 1
            if n != len(units): field_units = "dimensionless"
            return field_units
        except KeyError:
            return "dimensionless"

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
        if ds.events_data:
            for k,v in ds.events_info.items():
                fname = "event_"+k
                mylog.info("Adding field %s to the list of fields." % (fname))
                self.field_list.append(("io",fname))
                if k in ["x","y"]:
                    unit = "code_length"
                else:
                    unit = v
                self.parameter_file.field_units[("io",fname)] = unit
            return
        self._axis_map = {}
        self._file_map = {}
        self._ext_map = {}
        self._scale_map = {}
        # Since FITS header keywords are case-insensitive, we only pick a subset of
        # prefixes, ones that we expect to end up in headers.
        known_units = dict([(unit.lower(),unit) for unit in self.pf.unit_registry.lut])
        for unit in known_units.values():
            if unit in prefixable_units:
                for p in ["n","u","m","c","k"]:
                    known_units[(p+unit).lower()] = p+unit
        # We create a field from each slice on the 4th axis
        if self.parameter_file.naxis == 4:
            naxis4 = self.parameter_file.primary_header["naxis4"]
        else:
            naxis4 = 1
        for i, fits_file in enumerate(self.parameter_file._fits_files):
            for j, hdu in enumerate(fits_file):
                if self._ensure_same_dims(hdu):
                    units = self._determine_image_units(hdu.header, known_units)
                    try:
                        # Grab field name from btype
                        fname = hdu.header["btype"].lower()
                    except KeyError:
                        # Try to guess the name from the units
                        fname = self._guess_name_from_units(units)
                        # When all else fails
                        if fname is None: fname = "image_%d" % (j)
                    if self.pf.num_files > 1 and fname.startswith("image"):
                        fname += "_file_%d" % (i)
                    for k in xrange(naxis4):
                        if naxis4 > 1:
                            fname += "_%s_%d" % (hdu.header["CTYPE4"], k+1)
                        if fname in self.field_list:
                            mylog.error("You have two fields with the same name. Change one of " +
                                        "the names in the BTYPE header keyword to distinguish " +
                                        "them.")
                        self._axis_map[fname] = k
                        self._file_map[fname] = fits_file
                        self._ext_map[fname] = j
                        self._scale_map[fname] = [0.0,1.0]
                        if "bzero" in hdu.header:
                            self._scale_map[fname][0] = hdu.header["bzero"]
                        if "bscale" in hdu.header:
                            self._scale_map[fname][1] = hdu.header["bscale"]
                        self.field_list.append(("fits", fname))
                        self.parameter_file.field_units[fname] = units
                        mylog.info("Adding field %s to the list of fields." % (fname))
                        if units == "dimensionless":
                            mylog.warning("Could not determine dimensions for field %s, " % (fname) +
                                          "setting to dimensionless.")
                else:
                    mylog.warning("Image block %s does not have " % (hdu.name.lower()) +
                                  "the same dimensions as the primary and will not be " +
                                  "available as a field.")

        # For line fields, we still read the primary field. Not sure how to extend this
        # For now, we pick off the first field from the field list.
        line_db = self.parameter_file.line_database
        primary_fname = self.field_list[0][1]
        for k, v in line_db.iteritems():
            mylog.info("Adding line field: %s at frequency %g GHz" % (k, v))
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
            dims = np.array(pf.domain_dimensions)
            # If we are creating a dataset of lines, only decompose along the position axes
            if len(pf.line_database) > 0:
                dims[pf.vel_axis] = 1
            psize = get_psize(dims, pf.nprocs)
            gle, gre, shapes, slices = decompose_array(dims, psize, bbox)
            self.grid_left_edge = self.pf.arr(gle, "code_length")
            self.grid_right_edge = self.pf.arr(gre, "code_length")
            self.grid_dimensions = np.array([shape for shape in shapes], dtype="int32")
            # If we are creating a dataset of lines, only decompose along the position axes
            if len(pf.line_database) > 0:
                self.grid_left_edge[:,pf.vel_axis] = pf.domain_left_edge[pf.vel_axis]
                self.grid_right_edge[:,pf.vel_axis] = pf.domain_right_edge[pf.vel_axis]
                self.grid_dimensions[:,pf.vel_axis] = pf.domain_dimensions[pf.vel_axis]

        else:
            self.grid_left_edge[0,:] = pf.domain_left_edge
            self.grid_right_edge[0,:] = pf.domain_right_edge
            self.grid_dimensions[0] = pf.domain_dimensions

        if self.pf.events_data:
            try:
                self.grid_particle_count[:] = pf.primary_header["naxis2"]
            except KeyError:
                self.grid_particle_count[:] = 0.0
            self._particle_indices = np.zeros(self.num_grids + 1, dtype='int64')
            self._particle_indices[1] = self.grid_particle_count.squeeze()

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
                 auxiliary_files = [],
                 nprocs = None,
                 storage_filename = None,
                 nan_mask = None,
                 line_database = None,
                 line_width = None,
                 suppress_astropy_warnings = True,
                 parameters = None):

        if parameters is None:
            parameters = {}
        self.specified_parameters = parameters

        if line_width is not None:
            self.line_width = YTQuantity(line_width[0], line_width[1])
            self.line_units = line_width[1]
        else:
            self.line_width = None

        self.line_database = {}
        if line_database is not None:
            for k in line_database:
                self.line_database[k] = YTQuantity(line_database[k], self.line_units)

        if suppress_astropy_warnings:
            warnings.filterwarnings('ignore', module="astropy", append=True)
        auxiliary_files = ensure_list(auxiliary_files)
        self.filenames = [filename] + auxiliary_files
        self.num_files = len(self.filenames)
        self.fluid_types += ("fits",)
        if nan_mask is None:
            self.nan_mask = {}
        elif isinstance(nan_mask, float):
            self.nan_mask = {"all":nan_mask}
        elif isinstance(nan_mask, dict):
            self.nan_mask = nan_mask
        self.nprocs = nprocs
        self._handle = _astropy.pyfits.open(self.filenames[0],
                                      memmap=True,
                                      do_not_scale_image_data=True,
                                      ignore_blank=True)
        self._fits_files = [self._handle]
        if self.num_files > 1:
            for fits_file in auxiliary_files:
                if os.path.exists(fits_file):
                    fn = fits_file
                else:
                    fn = os.path.join(ytcfg.get("yt","test_data_dir"),fits_file)
                f = _astropy.pyfits.open(fn, memmap=True,
                                   do_not_scale_image_data=True,
                                   ignore_blank=True)
                self._fits_files.append(f)

        if len(self._handle) > 1 and self._handle[1].name == "EVENTS":
            self.events_data = True
            self.first_image = 1
            self.primary_header = self._handle[self.first_image].header
            self.naxis = 2
            self.wcs = _astropy.pywcs.WCS(naxis=2)
            self.events_info = {}
            for k,v in self.primary_header.items():
                if k.startswith("TTYP"):
                    if v.lower() in ["x","y"]:
                        num = k.strip("TTYPE")
                        self.events_info[v.lower()] = (self.primary_header["TLMIN"+num],
                                                       self.primary_header["TLMAX"+num],
                                                       self.primary_header["TCTYP"+num],
                                                       self.primary_header["TCRVL"+num],
                                                       self.primary_header["TCDLT"+num],
                                                       self.primary_header["TCRPX"+num])
                    elif v.lower() in ["energy","time"]:
                        num = k.strip("TTYPE")
                        unit = self.primary_header["TUNIT"+num].lower()
                        if unit.endswith("ev"): unit = unit.replace("ev","eV")
                        self.events_info[v.lower()] = unit
            self.axis_names = [self.events_info[ax][2] for ax in ["x","y"]]
            self.reblock = 1
            if "reblock" in self.specified_parameters:
                self.reblock = self.specified_parameters["reblock"]
            self.wcs.wcs.cdelt = [self.events_info["x"][4]*self.reblock,
                                  self.events_info["y"][4]*self.reblock]
            self.wcs.wcs.crpix = [(self.events_info["x"][5]-0.5)/self.reblock+0.5,
                                  (self.events_info["y"][5]-0.5)/self.reblock+0.5]
            self.wcs.wcs.ctype = [self.events_info["x"][2],self.events_info["y"][2]]
            self.wcs.wcs.cunit = ["deg","deg"]
            self.wcs.wcs.crval = [self.events_info["x"][3],self.events_info["y"][3]]
            self.dims = [(self.events_info["x"][1]-self.events_info["x"][0])/self.reblock,
                         (self.events_info["y"][1]-self.events_info["y"][0])/self.reblock]
        else:
            self.events_data = False
            self.first_image = 0
            self.primary_header = self._handle[self.first_image].header
            self.wcs = _astropy.pywcs.WCS(header=self.primary_header)
            self.naxis = self.primary_header["naxis"]
            self.axis_names = [self.primary_header["ctype%d" % (i+1)]
                               for i in xrange(self.naxis)]
            self.dims = [self.primary_header["naxis%d" % (i+1)]
                         for i in xrange(self.naxis)]

        self.refine_by = 2

        Dataset.__init__(self, filename, dataset_type)
        self.storage_filename = storage_filename

    def _set_code_unit_attributes(self):
        """
        Generates the conversion to various physical _units based on the parameter file
        """
        default_length_units = [u for u,v in default_unit_symbol_lut.items()
                                if str(v[-1]) == "(length)"]
        more_length_units = []
        for unit in default_length_units:
            if unit in prefixable_units:
                more_length_units += [prefix+unit for prefix in unit_prefixes]
        default_length_units += more_length_units
        file_units = []
        cunits = [self.wcs.wcs.cunit[i] for i in xrange(self.dimensionality)]
        for i, unit in enumerate(cunits):
            if unit in default_length_units:
                file_units.append(unit.name)
        if len(set(file_units)) == 1:
            length_factor = self.wcs.wcs.cdelt[0]
            length_unit = str(file_units[0])
            mylog.info("Found length units of %s." % (length_unit))
        else:
            self.no_cgs_equiv_length = True
            mylog.warning("No length conversion provided. Assuming 1 = 1 cm.")
            length_factor = 1.0
            length_unit = "cm"
        self.length_unit = self.quan(length_factor,length_unit)
        self.mass_unit = self.quan(1.0, "g")
        self.time_unit = self.quan(1.0, "s")
        self.velocity_unit = self.quan(1.0, "cm/s")
        if "beam_size" in self.specified_parameters:
            beam_size = self.quan(beam_size[0], beam_size[1]).in_cgs().value
        else:
            beam_size = 1.0
        self.unit_registry.add("beam",beam_size,dimensions=dimensions.solid_angle)
        if self.ppv_data:
            units = self.wcs_2d.wcs.cunit[0]
            if units == "deg": units = "degree"
            if units == "rad": units = "radian"
            pixel_area = np.prod(np.abs(self.wcs_2d.wcs.cdelt))
            pixel_area = self.quan(pixel_area, "%s**2" % (units)).in_cgs()
            pixel_dims = pixel_area.units.dimensions
            self.unit_registry.add("pixel",float(pixel_area.value),dimensions=pixel_dims)

    def _parse_parameter_file(self):

        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])

        # Determine dimensionality

        self.dimensionality = self.naxis
        self.geometry = "cartesian"

        # Sometimes a FITS file has a 4D datacube, in which case
        # we take the 4th axis and assume it consists of different fields.
        if self.dimensionality == 4: self.dimensionality = 3

        self.domain_dimensions = np.array(self.dims)[:self.dimensionality]
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

        # If this is a 2D events file, no need to decompose
        if self.events_data: self.nprocs = 1

        # If nprocs is None, do some automatic decomposition of the domain
        if self.nprocs is None:
            self.nprocs = np.around(np.prod(self.domain_dimensions) /
                                    32**self.dimensionality).astype("int")
            self.nprocs = max(min(self.nprocs, 512), 1)

        self.reversed = False

        # Check to see if this data is in some kind of (Lat,Lon,Vel) format
        self.ppv_data = False
        x = 0
        for p in lon_prefixes+lat_prefixes+vel_prefixes:
            y = np_char.startswith(self.axis_names[:self.dimensionality], p)
            x += y.sum()
        if x == self.dimensionality: self._setup_ppv()

    def _setup_ppv(self):

        self.ppv_data = True
        self.geometry = "ppv"

        end = min(self.dimensionality+1,4)
        if self.events_data:
            ctypes = self.axis_names
        else:
            ctypes = np.array([self.primary_header["CTYPE%d" % (i)]
                               for i in xrange(1,end)])

        log_str = "Detected these axes: "+"%s "*len(ctypes)
        mylog.info(log_str % tuple([ctype for ctype in ctypes]))

        self.lat_axis = np.zeros((end-1), dtype="bool")
        for p in lat_prefixes:
            self.lat_axis += np_char.startswith(ctypes, p)
        self.lat_axis = np.where(self.lat_axis)[0][0]
        self.lat_name = ctypes[self.lat_axis].split("-")[0].lower()

        self.lon_axis = np.zeros((end-1), dtype="bool")
        for p in lon_prefixes:
            self.lon_axis += np_char.startswith(ctypes, p)
        self.lon_axis = np.where(self.lon_axis)[0][0]
        self.lon_name = ctypes[self.lon_axis].split("-")[0].lower()

        if self.wcs.naxis > 2:

            self.vel_axis = np.zeros((end-1), dtype="bool")
            for p in vel_prefixes:
                self.vel_axis += np_char.startswith(ctypes, p)
            self.vel_axis = np.where(self.vel_axis)[0][0]
            self.vel_name = ctypes[self.vel_axis].split("-")[0].lower()

            self.wcs_2d = _astropy.pywcs.WCS(naxis=2)
            self.wcs_2d.wcs.crpix = self.wcs.wcs.crpix[[self.lon_axis, self.lat_axis]]
            self.wcs_2d.wcs.cdelt = self.wcs.wcs.cdelt[[self.lon_axis, self.lat_axis]]
            self.wcs_2d.wcs.crval = self.wcs.wcs.crval[[self.lon_axis, self.lat_axis]]
            self.wcs_2d.wcs.cunit = [str(self.wcs.wcs.cunit[self.lon_axis]),
                                     str(self.wcs.wcs.cunit[self.lat_axis])]
            self.wcs_2d.wcs.ctype = [self.wcs.wcs.ctype[self.lon_axis],
                                     self.wcs.wcs.ctype[self.lat_axis]]

            x0 = self.wcs.wcs.crpix[self.vel_axis]
            dz = self.wcs.wcs.cdelt[self.vel_axis]
            z0 = self.wcs.wcs.crval[self.vel_axis]
            self.vel_unit = str(self.wcs.wcs.cunit[self.vel_axis])

            if dz < 0.0:
                self.reversed = True
                le = self.dims[self.vel_axis]+0.5
                re = 0.5
            else:
                le = 0.5
                re = self.dims[self.vel_axis]+0.5
            self.domain_left_edge[self.vel_axis] = (le-x0)*dz + z0
            self.domain_right_edge[self.vel_axis] = (re-x0)*dz + z0
            if self.reversed: dz *= -1

            if self.line_width is not None:
                self.line_width = self.line_width.in_units(self.vel_unit)
                self.freq_begin = self.domain_left_edge[self.vel_axis]
                nz = np.rint(self.line_width.value/dz).astype("int")
                self.line_width = dz*nz
                self.domain_left_edge[self.vel_axis] = -self.line_width/2.
                self.domain_right_edge[self.vel_axis] = self.line_width/2.
                self.domain_dimensions[self.vel_axis] = nz

        else:

            self.wcs_2d = self.wcs
            self.vel_axis = 2
            self.vel_name = "z"
            self.vel_unit = "code length"

    def __del__(self):
        for file in self._fits_files:
            file.close()
            del file
        self._handle.close()
        del self._handle

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
                fileh = _astropy.pyfits.open(args[0])
            valid = fileh[0].header["naxis"] >= 2
            if len(fileh) > 1 and fileh[1].name == "EVENTS":
                valid = fileh[1].header["naxis"] >= 2
            fileh.close()
            return valid
        except:
            pass
        return False
