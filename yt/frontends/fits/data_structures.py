import os
import sys
import time
import uuid
import warnings
import weakref
from collections import defaultdict
from typing import Type

import numpy as np
import numpy.core.defchararray as np_char
from more_itertools import always_iterable

from yt.config import ytcfg
from yt.data_objects.index_subobjects.grid_patch import AMRGridPatch
from yt.data_objects.static_output import Dataset
from yt.fields.field_info_container import FieldInfoContainer
from yt.funcs import mylog, setdefaultattr
from yt.geometry.geometry_handler import YTDataChunk
from yt.geometry.grid_geometry_handler import GridIndex
from yt.units import dimensions
from yt.units.unit_lookup_table import (  # type: ignore
    default_unit_symbol_lut,
    unit_prefixes,
)
from yt.units.unit_object import UnitParseError  # type: ignore
from yt.units.yt_array import YTQuantity
from yt.utilities.decompose import decompose_array, get_psize
from yt.utilities.file_handler import FITSFileHandler
from yt.utilities.io_handler import io_registry
from yt.utilities.on_demand_imports import NotAModule, _astropy

from .fields import FITSFieldInfo, WCSFITSFieldInfo, YTFITSFieldInfo

if sys.version_info >= (3, 8):
    from functools import cached_property
else:
    from yt._maintenance.backports import cached_property
lon_prefixes = ["X", "RA", "GLON", "LINEAR"]
lat_prefixes = ["Y", "DEC", "GLAT", "LINEAR"]

spec_names = {"V": "Velocity", "F": "Frequency", "E": "Energy", "W": "Wavelength"}

space_prefixes = list(set(lon_prefixes + lat_prefixes))
unique_sky_prefixes = set(space_prefixes)
unique_sky_prefixes.difference_update({"X", "Y", "LINEAR"})
sky_prefixes = list(unique_sky_prefixes)
spec_prefixes = list(spec_names.keys())


class FITSGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, id, index, level):
        AMRGridPatch.__init__(self, id, filename=index.index_filename, index=index)
        self.Parent = None
        self.Children = []
        self.Level = 0


class FITSHierarchy(GridIndex):

    grid = FITSGrid

    def __init__(self, ds, dataset_type="fits"):
        self.dataset_type = dataset_type
        self.field_indexes = {}
        self.dataset = weakref.proxy(ds)
        # for now, the index file is the dataset
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        self._handle = ds._handle
        self.float_type = np.float64
        GridIndex.__init__(self, ds, dataset_type)

    def _initialize_data_storage(self):
        pass

    def _guess_name_from_units(self, units):
        field_from_unit = {"Jy": "intensity", "K": "temperature"}
        for k, v in field_from_unit.items():
            if k in units:
                mylog.warning(
                    "Guessing this is a %s field based on its units of %s.", v, k
                )
                return v
        return None

    def _determine_image_units(self, bunit):
        try:
            try:
                # First let AstroPy attempt to figure the unit out
                u = 1.0 * _astropy.units.Unit(bunit, format="fits")
                u = YTQuantity.from_astropy(u).units
            except ValueError:
                try:
                    # Let yt try it by itself
                    u = self.ds.quan(1.0, bunit).units
                except UnitParseError:
                    return "dimensionless"
            return str(u)
        except KeyError:
            return "dimensionless"

    def _ensure_same_dims(self, hdu):
        ds = self.dataset
        conditions = [hdu.header["naxis"] != ds.primary_header["naxis"]]
        for i in range(ds.naxis):
            nax = "naxis%d" % (i + 1)
            conditions.append(hdu.header[nax] != ds.primary_header[nax])
        if np.any(conditions):
            return False
        else:
            return True

    def _detect_output_fields(self):
        self.field_list = []
        self._axis_map = {}
        self._file_map = {}
        self._ext_map = {}
        self._scale_map = {}
        dup_field_index = {}
        # Since FITS header keywords are case-insensitive, we only pick a subset of
        # prefixes, ones that we expect to end up in headers.
        known_units = {unit.lower(): unit for unit in self.ds.unit_registry.lut}
        for unit in list(known_units.values()):
            if unit in self.ds.unit_registry.prefixable_units:
                for p in ["n", "u", "m", "c", "k"]:
                    known_units[(p + unit).lower()] = p + unit
        # We create a field from each slice on the 4th axis
        if self.dataset.naxis == 4:
            naxis4 = self.dataset.primary_header["naxis4"]
        else:
            naxis4 = 1
        for i, fits_file in enumerate(self.dataset._handle._fits_files):
            for j, hdu in enumerate(fits_file):
                if (
                    isinstance(hdu, _astropy.pyfits.BinTableHDU)
                    or hdu.header["naxis"] == 0
                ):
                    continue
                if self._ensure_same_dims(hdu):
                    units = self._determine_image_units(hdu.header["bunit"])
                    try:
                        # Grab field name from btype
                        fname = hdu.header["btype"]
                    except KeyError:
                        # Try to guess the name from the units
                        fname = self._guess_name_from_units(units)
                        # When all else fails
                        if fname is None:
                            fname = "image_%d" % (j)
                    if self.ds.num_files > 1 and fname.startswith("image"):
                        fname += "_file_%d" % (i)
                    if ("fits", fname) in self.field_list:
                        if fname in dup_field_index:
                            dup_field_index[fname] += 1
                        else:
                            dup_field_index[fname] = 1
                        mylog.warning(
                            "This field has the same name as a previously loaded "
                            "field. Changing the name from %s to %s_%d. To avoid "
                            "this, change one of the BTYPE header keywords.",
                            fname,
                            fname,
                            dup_field_index[fname],
                        )
                        fname += "_%d" % (dup_field_index[fname])
                    for k in range(naxis4):
                        if naxis4 > 1:
                            fname += "_%s_%d" % (hdu.header["CTYPE4"], k + 1)
                        self._axis_map[fname] = k
                        self._file_map[fname] = fits_file
                        self._ext_map[fname] = j
                        self._scale_map[fname] = [0.0, 1.0]
                        if "bzero" in hdu.header:
                            self._scale_map[fname][0] = hdu.header["bzero"]
                        if "bscale" in hdu.header:
                            self._scale_map[fname][1] = hdu.header["bscale"]
                        self.field_list.append(("fits", fname))
                        self.dataset.field_units[fname] = units
                        mylog.info("Adding field %s to the list of fields.", fname)
                        if units == "dimensionless":
                            mylog.warning(
                                "Could not determine dimensions for field %s, "
                                "setting to dimensionless.",
                                fname,
                            )
                else:
                    mylog.warning(
                        "Image block %s does not have the same dimensions "
                        "as the primary and will not be available as a field.",
                        hdu.name.lower(),
                    )

    def _count_grids(self):
        self.num_grids = self.ds.parameters["nprocs"]

    def _parse_index(self):
        ds = self.dataset

        # If nprocs > 1, decompose the domain into virtual grids
        if self.num_grids > 1:
            self._domain_decomp()
        else:
            self.grid_left_edge[0, :] = ds.domain_left_edge
            self.grid_right_edge[0, :] = ds.domain_right_edge
            self.grid_dimensions[0] = ds.domain_dimensions

        self.grid_levels.flat[:] = 0
        self.grids = np.empty(self.num_grids, dtype="object")
        for i in range(self.num_grids):
            self.grids[i] = self.grid(i, self, self.grid_levels[i, 0])

    def _domain_decomp(self):
        bbox = np.array(
            [self.ds.domain_left_edge, self.ds.domain_right_edge]
        ).transpose()
        dims = self.ds.domain_dimensions
        psize = get_psize(dims, self.num_grids)
        gle, gre, shapes, slices = decompose_array(dims, psize, bbox)
        self.grid_left_edge = self.ds.arr(gle, "code_length")
        self.grid_right_edge = self.ds.arr(gre, "code_length")
        self.grid_dimensions = np.array(shapes, dtype="int32")

    def _populate_grid_objects(self):
        for i in range(self.num_grids):
            self.grids[i]._prepare_grid()
            self.grids[i]._setup_dx()
        self.max_level = 0

    def _setup_derived_fields(self):
        super()._setup_derived_fields()
        [self.dataset.conversion_factors[field] for field in self.field_list]
        for field in self.field_list:
            if field not in self.derived_field_list:
                self.derived_field_list.append(field)

        for field in self.derived_field_list:
            f = self.dataset.field_info[field]
            if f.is_alias:
                # Translating an already-converted field
                self.dataset.conversion_factors[field] = 1.0

    def _setup_data_io(self):
        self.io = io_registry[self.dataset_type](self.dataset)

    def _chunk_io(self, dobj, cache=True, local_only=False):
        # local_only is only useful for inline datasets and requires
        # implementation by subclasses.
        gfiles = defaultdict(list)
        gobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for g in gobjs:
            gfiles[g.id].append(g)
        for fn in sorted(gfiles):
            gs = gfiles[fn]
            yield YTDataChunk(
                dobj, "io", gs, self._count_selection(dobj, gs), cache=cache
            )


def find_primary_header(fileh):
    # Sometimes the primary hdu doesn't have an image
    if len(fileh) > 1 and fileh[0].header["naxis"] == 0:
        first_image = 1
    else:
        first_image = 0
    header = fileh[first_image].header
    return header, first_image


def check_fits_valid(filename):
    ext = filename.rsplit(".", 1)[-1]
    if ext.upper() in ("GZ", "FZ"):
        # We don't know for sure that there will be > 1
        ext = filename.rsplit(".", 1)[0].rsplit(".", 1)[-1]
    if ext.upper() not in ("FITS", "FTS"):
        return None
    elif isinstance(_astropy.pyfits, NotAModule):
        raise RuntimeError(
            "This appears to be a FITS file, but AstroPy is not installed."
        )
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning, append=True)
            fileh = _astropy.pyfits.open(filename)
        header, _ = find_primary_header(fileh)
        if header["naxis"] >= 2:
            return fileh
        else:
            fileh.close()
    except Exception:
        pass
    return None


def check_sky_coords(filename, ndim):
    fileh = check_fits_valid(filename)
    if fileh is not None:
        try:
            if len(fileh) > 1 and fileh[1].name == "EVENTS" and ndim == 2:
                fileh.close()
                return True
            else:
                header, _ = find_primary_header(fileh)
                if header["naxis"] < ndim:
                    return False
                axis_names = [
                    header.get("ctype%d" % (i + 1), "") for i in range(header["naxis"])
                ]
                if len(axis_names) == 3 and axis_names.count("LINEAR") == 2:
                    return any(a[0] in spec_prefixes for a in axis_names)
                x = find_axes(axis_names, sky_prefixes + spec_prefixes)
                fileh.close()
                return x >= ndim
        except Exception:
            pass
    return False


class FITSDataset(Dataset):
    _index_class = FITSHierarchy
    _field_info_class: Type[FieldInfoContainer] = FITSFieldInfo
    _dataset_type = "fits"
    _handle = None

    def __init__(
        self,
        filename,
        dataset_type="fits",
        auxiliary_files=None,
        nprocs=None,
        storage_filename=None,
        nan_mask=None,
        suppress_astropy_warnings=True,
        parameters=None,
        units_override=None,
        unit_system="cgs",
    ):

        if parameters is None:
            parameters = {}
        parameters["nprocs"] = nprocs
        self.specified_parameters = parameters

        if suppress_astropy_warnings:
            warnings.filterwarnings("ignore", module="astropy", append=True)

        self.filenames = [filename] + list(always_iterable(auxiliary_files))
        self.num_files = len(self.filenames)
        self.fluid_types += ("fits",)
        if nan_mask is None:
            self.nan_mask = {}
        elif isinstance(nan_mask, float):
            self.nan_mask = {"all": nan_mask}
        elif isinstance(nan_mask, dict):
            self.nan_mask = nan_mask
        self._handle = FITSFileHandler(self.filenames[0])
        if isinstance(
            self.filenames[0], _astropy.pyfits.hdu.image._ImageBaseHDU
        ) or isinstance(self.filenames[0], _astropy.pyfits.HDUList):
            fn = f"InMemoryFITSFile_{uuid.uuid4().hex}"
        else:
            fn = self.filenames[0]
        self._handle._fits_files.append(self._handle)
        if self.num_files > 1:
            for fits_file in auxiliary_files:
                if isinstance(fits_file, _astropy.pyfits.hdu.image._ImageBaseHDU):
                    f = _astropy.pyfits.HDUList([fits_file])
                elif isinstance(fits_file, _astropy.pyfits.HDUList):
                    f = fits_file
                else:
                    if os.path.exists(fits_file):
                        fn = fits_file
                    else:
                        fn = os.path.join(ytcfg.get("yt", "test_data_dir"), fits_file)
                    f = _astropy.pyfits.open(
                        fn, memmap=True, do_not_scale_image_data=True, ignore_blank=True
                    )
                self._handle._fits_files.append(f)

        self.refine_by = 2

        Dataset.__init__(
            self,
            fn,
            dataset_type,
            units_override=units_override,
            unit_system=unit_system,
        )
        self.storage_filename = storage_filename

    def _set_code_unit_attributes(self):
        """
        Generates the conversion to various physical _units based on the
        parameter file
        """
        if getattr(self, "length_unit", None) is None:
            default_length_units = [
                u for u, v in default_unit_symbol_lut.items() if str(v[1]) == "(length)"
            ]
            more_length_units = []
            for unit in default_length_units:
                if unit in self.unit_registry.prefixable_units:
                    more_length_units += [prefix + unit for prefix in unit_prefixes]
            default_length_units += more_length_units
            file_units = []
            cunits = [self.wcs.wcs.cunit[i] for i in range(self.dimensionality)]
            for unit in (_.to_string() for _ in cunits):
                if unit in default_length_units:
                    file_units.append(unit)
            if len(set(file_units)) == 1:
                length_factor = self.wcs.wcs.cdelt[0]
                length_unit = str(file_units[0])
                mylog.info("Found length units of %s.", length_unit)
            else:
                self.no_cgs_equiv_length = True
                mylog.warning("No length conversion provided. Assuming 1 = 1 cm.")
                length_factor = 1.0
                length_unit = "cm"
            setdefaultattr(self, "length_unit", self.quan(length_factor, length_unit))
        for unit, cgs in [("time", "s"), ("mass", "g")]:
            # We set these to cgs for now, but they may have been overridden
            if getattr(self, unit + "_unit", None) is not None:
                continue
            mylog.warning("Assuming 1.0 = 1.0 %s", cgs)
            setdefaultattr(self, f"{unit}_unit", self.quan(1.0, cgs))
        self.magnetic_unit = np.sqrt(
            4 * np.pi * self.mass_unit / (self.time_unit**2 * self.length_unit)
        )
        self.magnetic_unit.convert_to_units("gauss")
        self.velocity_unit = self.length_unit / self.time_unit

    @cached_property
    def unique_identifier(self) -> str:
        if self.parameter_filename.startswith("InMemory"):
            return str(time.time())
        else:
            return super().unique_identifier

    def _parse_parameter_file(self):

        self._determine_structure()
        self._determine_axes()

        # Determine dimensionality

        self.dimensionality = self.naxis
        self.geometry = "cartesian"

        # Sometimes a FITS file has a 4D datacube, in which case
        # we take the 4th axis and assume it consists of different fields.
        if self.dimensionality == 4:
            self.dimensionality = 3

        self._determine_wcs()

        self.current_time = 0.0

        self.domain_dimensions = np.array(self.dims)[: self.dimensionality]
        if self.dimensionality == 2:
            self.domain_dimensions = np.append(self.domain_dimensions, [int(1)])
        self._determine_bbox()

        # Get the simulation time
        try:
            self.current_time = self.parameters["time"]
        except Exception:
            mylog.warning("Cannot find time")
            self.current_time = 0.0
            pass

        # For now we'll ignore these
        self._periodicity = (False,) * 3
        self.current_redshift = 0.0
        self.omega_lambda = 0.0
        self.omega_matter = 0.0
        self.hubble_constant = 0.0
        self.cosmological_simulation = 0

        self._determine_nprocs()

        # Now we can set up some of our parameters for convenience.
        for k, v in self.primary_header.items():
            self.parameters[k] = v
        # Remove potential default keys
        self.parameters.pop("", None)

    def _determine_nprocs(self):
        # If nprocs is None, do some automatic decomposition of the domain
        if self.specified_parameters["nprocs"] is None:
            nprocs = np.around(
                np.prod(self.domain_dimensions) / 32**self.dimensionality
            ).astype("int")
            self.parameters["nprocs"] = max(min(nprocs, 512), 1)
        else:
            self.parameters["nprocs"] = self.specified_parameters["nprocs"]

    def _determine_structure(self):
        self.primary_header, self.first_image = find_primary_header(self._handle)
        self.naxis = self.primary_header["naxis"]
        self.axis_names = [
            self.primary_header.get("ctype%d" % (i + 1), "LINEAR")
            for i in range(self.naxis)
        ]
        self.dims = [
            self.primary_header["naxis%d" % (i + 1)] for i in range(self.naxis)
        ]

    def _determine_wcs(self):
        wcs = _astropy.pywcs.WCS(header=self.primary_header)
        if self.naxis == 4:
            self.wcs = _astropy.pywcs.WCS(naxis=3)
            self.wcs.wcs.crpix = wcs.wcs.crpix[:3]
            self.wcs.wcs.cdelt = wcs.wcs.cdelt[:3]
            self.wcs.wcs.crval = wcs.wcs.crval[:3]
            self.wcs.wcs.cunit = [str(unit) for unit in wcs.wcs.cunit][:3]
            self.wcs.wcs.ctype = [type for type in wcs.wcs.ctype][:3]
        else:
            self.wcs = wcs

    def _determine_bbox(self):
        domain_left_edge = np.array([0.5] * 3)
        domain_right_edge = np.array(
            [float(dim) + 0.5 for dim in self.domain_dimensions]
        )

        if self.dimensionality == 2:
            domain_left_edge[-1] = 0.5
            domain_right_edge[-1] = 1.5

        self.domain_left_edge = domain_left_edge
        self.domain_right_edge = domain_right_edge

    def _determine_axes(self):
        self.lat_axis = 1
        self.lon_axis = 0
        self.lat_name = "Y"
        self.lon_name = "X"

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        try:
            fileh = check_fits_valid(filename)
        except Exception:
            return False

        if fileh is None:
            return False
        else:
            fileh.close()
            return True

    @classmethod
    def _guess_candidates(cls, base, directories, files):
        candidates = []
        for fn, fnl in ((_, _.lower()) for _ in files):
            if (
                fnl.endswith(".fits")
                or fnl.endswith(".fits.gz")
                or fnl.endswith(".fits.fz")
            ):
                candidates.append(fn)
        # FITS files don't preclude subdirectories
        return candidates, True

    def close(self):
        self._handle.close()


def find_axes(axis_names, prefixes):
    x = 0
    for p in prefixes:
        y = np_char.startswith(axis_names, p)
        x += np.any(y)
    return x


class YTFITSDataset(FITSDataset):
    _field_info_class = YTFITSFieldInfo

    def _parse_parameter_file(self):
        super()._parse_parameter_file()
        # Get the current time
        if "time" in self.primary_header:
            self.current_time = self.primary_header["time"]

    def _set_code_unit_attributes(self):
        """
        Generates the conversion to various physical _units based on the parameter file
        """
        for unit, cgs in [
            ("length", "cm"),
            ("time", "s"),
            ("mass", "g"),
            ("velocity", "cm/s"),
            ("magnetic", "gauss"),
        ]:
            if unit == "magnetic":
                short_unit = "bfunit"
            else:
                short_unit = f"{unit[0]}unit"
            if short_unit in self.primary_header:
                # units should now be in header
                u = self.quan(
                    self.primary_header[short_unit],
                    self.primary_header.comments[short_unit].strip("[]"),
                )
                mylog.info("Found %s units of %s.", unit, u)
            else:
                if unit == "length":
                    # Falling back to old way of getting units for length
                    # in old files
                    u = self.quan(1.0, str(self.wcs.wcs.cunit[0]))
                    mylog.info("Found %s units of %s.", unit, u)
                else:
                    # Give up otherwise
                    u = self.quan(1.0, cgs)
                    mylog.warning(
                        "No unit for %s found. Assuming 1.0 code_%s = 1.0 %s",
                        unit,
                        unit,
                        cgs,
                    )
            setdefaultattr(self, f"{unit}_unit", u)

    def _determine_bbox(self):
        dx = np.zeros(3)
        dx[: self.dimensionality] = self.wcs.wcs.cdelt
        domain_left_edge = np.zeros(3)
        domain_left_edge[: self.dimensionality] = self.wcs.wcs.crval - dx[
            : self.dimensionality
        ] * (self.wcs.wcs.crpix - 0.5)
        domain_right_edge = domain_left_edge + dx * self.domain_dimensions

        if self.dimensionality == 2:
            domain_left_edge[-1] = 0.0
            domain_right_edge[-1] = dx[0]

        self.domain_left_edge = domain_left_edge
        self.domain_right_edge = domain_right_edge

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        try:
            fileh = check_fits_valid(filename)
        except Exception:
            return False

        if fileh is None:
            return False
        else:
            if "WCSNAME" in fileh[0].header:
                isyt = fileh[0].header["WCSNAME"].strip() == "yt"
            else:
                isyt = False
            fileh.close()
            return isyt


class SkyDataFITSDataset(FITSDataset):
    _field_info_class = WCSFITSFieldInfo

    def _determine_wcs(self):
        super()._determine_wcs()
        end = min(self.dimensionality + 1, 4)
        self.ctypes = np.array(
            [self.primary_header["CTYPE%d" % (i)] for i in range(1, end)]
        )
        self.wcs_2d = self.wcs

    def _parse_parameter_file(self):
        super()._parse_parameter_file()

        end = min(self.dimensionality + 1, 4)

        self.geometry = "spectral_cube"

        log_str = "Detected these axes: " + "%s " * len(self.ctypes)
        mylog.info(log_str, *self.ctypes)

        self.lat_axis = np.zeros((end - 1), dtype="bool")
        for p in lat_prefixes:
            self.lat_axis += np_char.startswith(self.ctypes, p)
        self.lat_axis = np.where(self.lat_axis)[0][0]
        self.lat_name = self.ctypes[self.lat_axis].split("-")[0].lower()

        self.lon_axis = np.zeros((end - 1), dtype="bool")
        for p in lon_prefixes:
            self.lon_axis += np_char.startswith(self.ctypes, p)
        self.lon_axis = np.where(self.lon_axis)[0][0]
        self.lon_name = self.ctypes[self.lon_axis].split("-")[0].lower()

        if self.lat_axis == self.lon_axis and self.lat_name == self.lon_name:
            self.lat_axis = 1
            self.lon_axis = 0
            self.lat_name = "Y"
            self.lon_name = "X"

        self.spec_axis = 2
        self.spec_name = "z"
        self.spec_unit = ""

    def _set_code_unit_attributes(self):
        super()._set_code_unit_attributes()
        units = self.wcs_2d.wcs.cunit[0]
        if units == "deg":
            units = "degree"
        if units == "rad":
            units = "radian"
        pixel_area = np.prod(np.abs(self.wcs_2d.wcs.cdelt))
        pixel_area = self.quan(pixel_area, f"{units}**2").in_cgs()
        pixel_dims = pixel_area.units.dimensions
        self.unit_registry.add("pixel", float(pixel_area.value), dimensions=pixel_dims)
        if "beam_size" in self.specified_parameters:
            beam_size = self.specified_parameters["beam_size"]
            beam_size = self.quan(beam_size[0], beam_size[1]).in_cgs().value
            self.unit_registry.add("beam", beam_size, dimensions=dimensions.solid_angle)

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        try:
            return check_sky_coords(filename, ndim=2)
        except Exception:
            return False


class SpectralCubeFITSHierarchy(FITSHierarchy):
    def _domain_decomp(self):
        dz = self.ds.quan(1.0, "code_length") * self.ds.spectral_factor
        self.grid_dimensions[:, 2] = np.around(
            float(self.ds.domain_dimensions[2]) / self.num_grids
        ).astype("int")
        self.grid_dimensions[-1, 2] += self.ds.domain_dimensions[2] % self.num_grids
        self.grid_left_edge[0, 2] = self.ds.domain_left_edge[2]
        self.grid_left_edge[1:, 2] = (
            self.ds.domain_left_edge[2] + np.cumsum(self.grid_dimensions[:-1, 2]) * dz
        )
        self.grid_right_edge[:, 2] = (
            self.grid_left_edge[:, 2] + self.grid_dimensions[:, 2] * dz
        )
        self.grid_left_edge[:, :2] = self.ds.domain_left_edge[:2]
        self.grid_right_edge[:, :2] = self.ds.domain_right_edge[:2]
        self.grid_dimensions[:, :2] = self.ds.domain_dimensions[:2]


class SpectralCubeFITSDataset(SkyDataFITSDataset):
    _index_class = SpectralCubeFITSHierarchy

    def __init__(
        self,
        filename,
        auxiliary_files=None,
        nprocs=None,
        storage_filename=None,
        nan_mask=None,
        spectral_factor=1.0,
        suppress_astropy_warnings=True,
        parameters=None,
        units_override=None,
        unit_system="cgs",
    ):
        if auxiliary_files is None:
            auxiliary_files = []
        self.spectral_factor = spectral_factor
        super().__init__(
            filename,
            nprocs=nprocs,
            auxiliary_files=auxiliary_files,
            storage_filename=storage_filename,
            suppress_astropy_warnings=suppress_astropy_warnings,
            nan_mask=nan_mask,
            parameters=parameters,
            units_override=units_override,
            unit_system=unit_system,
        )

    def _parse_parameter_file(self):
        super()._parse_parameter_file()

        self.geometry = "spectral_cube"

        end = min(self.dimensionality + 1, 4)

        self.spec_axis = np.zeros(end - 1, dtype="bool")
        for p in spec_names.keys():
            self.spec_axis += np_char.startswith(self.ctypes, p)
        self.spec_axis = np.where(self.spec_axis)[0][0]
        self.spec_name = spec_names[self.ctypes[self.spec_axis].split("-")[0][0]]

        # Extract a subimage from a WCS object
        self.wcs_2d = self.wcs.sub(["longitude", "latitude"])

        self._p0 = self.wcs.wcs.crpix[self.spec_axis]
        self._dz = self.wcs.wcs.cdelt[self.spec_axis]
        self._z0 = self.wcs.wcs.crval[self.spec_axis]
        self.spec_unit = str(self.wcs.wcs.cunit[self.spec_axis])

        if self.spectral_factor == "auto":
            self.spectral_factor = float(
                max(self.domain_dimensions[[self.lon_axis, self.lat_axis]])
            )
            self.spectral_factor /= self.domain_dimensions[self.spec_axis]
            mylog.info("Setting the spectral factor to %f", self.spectral_factor)
        Dz = (
            self.domain_right_edge[self.spec_axis]
            - self.domain_left_edge[self.spec_axis]
        )
        dre = self.domain_right_edge
        dre[self.spec_axis] = (
            self.domain_left_edge[self.spec_axis] + self.spectral_factor * Dz
        )
        self.domain_right_edge = dre
        self._dz /= self.spectral_factor
        self._p0 = (self._p0 - 0.5) * self.spectral_factor + 0.5

    def _determine_nprocs(self):
        # If nprocs is None, do some automatic decomposition of the domain
        if self.specified_parameters["nprocs"] is None:
            nprocs = np.around(self.domain_dimensions[2] / 8).astype("int")
            self.parameters["nprocs"] = max(min(nprocs, 512), 1)
        else:
            self.parameters["nprocs"] = self.specified_parameters["nprocs"]

    def spec2pixel(self, spec_value):
        sv = self.arr(spec_value).in_units(self.spec_unit)
        return self.arr((sv.v - self._z0) / self._dz + self._p0, "code_length")

    def pixel2spec(self, pixel_value):
        pv = self.arr(pixel_value, "code_length")
        return self.arr((pv.v - self._p0) * self._dz + self._z0, self.spec_unit)

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        try:
            return check_sky_coords(filename, ndim=3)
        except Exception:
            return False


class EventsFITSHierarchy(FITSHierarchy):
    def _detect_output_fields(self):
        ds = self.dataset
        self.field_list = []
        for k, v in ds.events_info.items():
            fname = "event_" + k
            mylog.info("Adding field %s to the list of fields.", fname)
            self.field_list.append(("io", fname))
            if k in ["x", "y"]:
                field_unit = "code_length"
            else:
                field_unit = v
            self.dataset.field_units[("io", fname)] = field_unit
        return

    def _parse_index(self):
        super()._parse_index()
        try:
            self.grid_particle_count[:] = self.dataset.primary_header["naxis2"]
        except KeyError:
            self.grid_particle_count[:] = 0.0
        self._particle_indices = np.zeros(self.num_grids + 1, dtype="int64")
        self._particle_indices[1] = self.grid_particle_count.squeeze()


class EventsFITSDataset(SkyDataFITSDataset):
    _index_class = EventsFITSHierarchy

    def __init__(
        self,
        filename,
        storage_filename=None,
        suppress_astropy_warnings=True,
        reblock=1,
        parameters=None,
        units_override=None,
        unit_system="cgs",
    ):
        self.reblock = reblock
        super().__init__(
            filename,
            nprocs=1,
            storage_filename=storage_filename,
            parameters=parameters,
            suppress_astropy_warnings=suppress_astropy_warnings,
            units_override=units_override,
            unit_system=unit_system,
        )

    def _determine_structure(self):
        self.first_image = 1
        self.primary_header = self._handle[self.first_image].header
        self.naxis = 2

    def _determine_wcs(self):
        self.wcs = _astropy.pywcs.WCS(naxis=2)
        self.events_info = {}
        for k, v in self.primary_header.items():
            if k.startswith("TTYP"):
                if v.lower() in ["x", "y"]:
                    num = k.replace("TTYPE", "")
                    self.events_info[v.lower()] = (
                        self.primary_header["TLMIN" + num],
                        self.primary_header["TLMAX" + num],
                        self.primary_header["TCTYP" + num],
                        self.primary_header["TCRVL" + num],
                        self.primary_header["TCDLT" + num],
                        self.primary_header["TCRPX" + num],
                    )
                elif v.lower() in ["energy", "time"]:
                    num = k.replace("TTYPE", "")
                    unit = self.primary_header["TUNIT" + num].lower()
                    if unit.endswith("ev"):
                        unit = unit.replace("ev", "eV")
                    self.events_info[v.lower()] = unit
        self.axis_names = [self.events_info[ax][2] for ax in ["x", "y"]]
        self.wcs.wcs.cdelt = [
            self.events_info["x"][4] * self.reblock,
            self.events_info["y"][4] * self.reblock,
        ]
        self.wcs.wcs.crpix = [
            (self.events_info["x"][5] - 0.5) / self.reblock + 0.5,
            (self.events_info["y"][5] - 0.5) / self.reblock + 0.5,
        ]
        self.wcs.wcs.ctype = [self.events_info["x"][2], self.events_info["y"][2]]
        self.wcs.wcs.cunit = ["deg", "deg"]
        self.wcs.wcs.crval = [self.events_info["x"][3], self.events_info["y"][3]]
        self.dims = [
            (self.events_info["x"][1] - self.events_info["x"][0]) / self.reblock,
            (self.events_info["y"][1] - self.events_info["y"][0]) / self.reblock,
        ]
        self.ctypes = self.axis_names
        self.wcs_2d = self.wcs

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        try:
            fileh = check_fits_valid(filename)
        except Exception:
            return False
        if fileh is not None:
            try:
                valid = fileh[1].name == "EVENTS"
                fileh.close()
                return valid
            except Exception:
                pass
        return False
