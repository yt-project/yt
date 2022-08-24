import os
import struct
from typing import Type

import numpy as np

from yt.data_objects.static_output import ParticleFile
from yt.fields.field_info_container import FieldInfoContainer
from yt.frontends.sph.data_structures import SPHDataset, SPHParticleIndex
from yt.funcs import only_on_root
from yt.geometry.geometry_handler import Index
from yt.utilities.chemical_formulas import compute_mu
from yt.utilities.cosmology import Cosmology
from yt.utilities.fortran_utils import read_record
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.on_demand_imports import _h5py as h5py

from .definitions import (
    SNAP_FORMAT_2_OFFSET,
    gadget_field_specs,
    gadget_header_specs,
    gadget_ptype_specs,
)
from .fields import GadgetFieldInfo


def _fix_unit_ordering(unit):
    if isinstance(unit[0], str):
        unit = unit[1], unit[0]
    return unit


def _byte_swap_32(x):
    return struct.unpack(">I", struct.pack("<I", x))[0]


class GadgetBinaryHeader:
    """A convenient interface to Gadget binary header.

    This is a helper class to facilitate the main dataset and IO classes.
    The main usage is through the GadgetDataset._header attribute. It is also
    used stand-alone in GadgetDataset._is_valid method.
    """

    _placeholder_keys = ["unused", "empty"]

    def __init__(self, filename, header_spec):
        self.filename = filename
        if isinstance(header_spec, str):
            header_spec = [header_spec]
        self.spec = [
            GadgetDataset._setup_binary_spec(hs, gadget_header_specs)
            for hs in header_spec
        ]

    @property
    def float_type(self):
        """Determine whether the floats are single or double precision."""
        gformat, endianswap = self.gadget_format
        # We do not currently support double precision snapshot format 2
        if gformat == 2:
            return "f4"
        else:
            hvals = self.value
            np0 = sum(hvals["Npart"])
            with self.open() as f:
                f.seek(self.position_offset)
                # Calculate particle number assuming single precision
                np1 = struct.unpack(endianswap + "I", f.read(4))[0] / (4 * 3)
            if np1 == np0:
                return "f4"
            elif np1 == 2 * np0:
                return "f8"
            else:
                raise RuntimeError(
                    "Gadget snapshot file is likely corrupted! "
                    "Cannot determine float type."
                )

    @property
    def gadget_format(self):
        """Determine Gadget snapshot format and endianness.

        The difference between Gadget snapshot format 1 and 2 can be told from
        the first 4 bytes of the file. For format 1, it's the header size. For
        format 2, it's always 8.
        """
        first_header_size = self.size[0]
        # Read the first 4 bytes assuming little endian int32
        with self.open() as f:
            (rhead,) = struct.unpack("<I", f.read(4))
        # Format 1?
        if rhead == first_header_size:
            return 1, "<"
        elif rhead == _byte_swap_32(first_header_size):
            return 1, ">"
        # Format 2?
        elif rhead == 8:
            return 2, "<"
        elif rhead == _byte_swap_32(8):
            return 2, ">"
        else:
            raise RuntimeError(
                "Gadget snapshot file is likely corrupted! "
                "The first 4 bytes represent %s (as little endian int32). "
                "But we are looking for %s (for format 1) or 8 (for format 2)."
                % (rhead, first_header_size)
            )

    @property
    def position_offset(self):
        """Offset to the position block."""
        n_header = len(self.size)
        offset = 8 * n_header + sum(self.size)
        if self.gadget_format[0] == 2:
            offset += SNAP_FORMAT_2_OFFSET * (n_header + 1)
        return offset

    @property
    def size(self):
        """Header size in bytes."""

        def single_header_size(single_header_spec):
            return sum(
                field[1] * np.dtype(field[2]).itemsize for field in single_header_spec
            )

        return [single_header_size(spec) for spec in self.spec]

    @property
    def value(self):
        """Header values as a dictionary."""
        # The entries in this header are capitalized and named to match Table 4
        # in the GADGET-2 user guide.
        gformat, endianswap = self.gadget_format
        # Read header
        with self.open() as f:
            hvals = {}
            for spec in self.spec:
                if gformat == 2:
                    f.seek(f.tell() + SNAP_FORMAT_2_OFFSET)
                hvals_new = read_record(f, spec, endian=endianswap)
                hvals.update(hvals_new)
        # Remove placeholder keys
        for key in self._placeholder_keys:
            if key in hvals:
                del hvals[key]
        # Convert length 1 list to scalar value
        for i in hvals:
            if len(hvals[i]) == 1:
                hvals[i] = hvals[i][0]
        return hvals

    def open(self):
        """Open snapshot file."""
        for filename in [self.filename, self.filename + ".0"]:
            if os.path.exists(filename):
                return open(filename, "rb")
        raise FileNotFoundError(f"Snapshot file {self.filename} does not exist.")

    def validate(self):
        """Validate data integrity."""
        try:
            self.open().close()
            self.gadget_format
            self.float_type
        except Exception:
            return False
        return True


class GadgetBinaryFile(ParticleFile):
    def __init__(self, ds, io, filename, file_id, range=None):
        header = GadgetBinaryHeader(filename, ds._header.spec)
        self.header = header.value
        self._position_offset = header.position_offset
        with header.open() as f:
            self._file_size = f.seek(0, os.SEEK_END)

        super().__init__(ds, io, filename, file_id, range)

    def _calculate_offsets(self, field_list, pcounts):
        # Note that we ignore pcounts here because it's the global count.  We
        # just want the local count, which we store here.
        self.field_offsets = self.io._calculate_field_offsets(
            field_list,
            self.header["Npart"].copy(),
            self._position_offset,
            self.start,
            self._file_size,
        )


class GadgetBinaryIndex(SPHParticleIndex):
    def __init__(self, ds, dataset_type):
        super().__init__(ds, dataset_type)
        self._initialize_index()

    def _initialize_index(self):
        # Normally this function is called during field detection. We call it
        # here because we need to know which fields exist on-disk so that we can
        # read in the smoothing lengths for SPH data before we construct the
        # Morton bitmaps.
        self._detect_output_fields()
        super()._initialize_index()

    def _initialize_frontend_specific(self):
        super()._initialize_frontend_specific()
        self.io._float_type = self.ds._header.float_type


class GadgetDataset(SPHDataset):
    _index_class: Type[Index] = GadgetBinaryIndex
    _file_class: Type[ParticleFile] = GadgetBinaryFile
    _field_info_class: Type[FieldInfoContainer] = GadgetFieldInfo
    _particle_mass_name = "Mass"
    _particle_coordinates_name = "Coordinates"
    _particle_velocity_name = "Velocities"
    _sph_ptypes = ("Gas",)
    _suffix = ""

    def __init__(
        self,
        filename,
        dataset_type="gadget_binary",
        additional_fields=(),
        unit_base=None,
        index_order=None,
        index_filename=None,
        kdtree_filename=None,
        kernel_name=None,
        bounding_box=None,
        header_spec="default",
        field_spec="default",
        ptype_spec="default",
        long_ids=False,
        units_override=None,
        mean_molecular_weight=None,
        header_offset=0,
        unit_system="cgs",
        use_dark_factor=False,
        w_0=-1.0,
        w_a=0.0,
        default_species_fields=None,
    ):
        if self._instantiated:
            return
        # Check if filename is a directory
        if os.path.isdir(filename):
            # Get the .0 snapshot file. We know there's only 1 and it's valid since we
            # came through _is_valid in load()
            for f in os.listdir(filename):
                fname = os.path.join(filename, f)
                fext = os.path.splitext(fname)[-1]
                if (
                    (".0" in f)
                    and (fext not in {".ewah", ".kdtree"})
                    and os.path.isfile(fname)
                ):
                    filename = os.path.join(filename, f)
                    break
        self._header = GadgetBinaryHeader(filename, header_spec)
        header_size = self._header.size
        if header_size != [256]:
            only_on_root(
                mylog.warning,
                "Non-standard header size is detected! "
                "Gadget-2 standard header is 256 bytes, but yours is %s. "
                "Make sure a non-standard header is actually expected. "
                "Otherwise something is wrong, "
                "and you might want to check how the dataset is loaded. "
                "Further information about header specification can be found in "
                "https://yt-project.org/docs/dev/examining/loading_data.html#header-specification.",
                header_size,
            )
        self._field_spec = self._setup_binary_spec(field_spec, gadget_field_specs)
        self._ptype_spec = self._setup_binary_spec(ptype_spec, gadget_ptype_specs)
        self.storage_filename = None
        if long_ids:
            self._id_dtype = "u8"
        else:
            self._id_dtype = "u4"
        self.long_ids = long_ids
        self.header_offset = header_offset
        if unit_base is not None and "UnitLength_in_cm" in unit_base:
            # We assume this is comoving, because in the absence of comoving
            # integration the redshift will be zero.
            unit_base["cmcm"] = 1.0 / unit_base["UnitLength_in_cm"]
        self._unit_base = unit_base
        if bounding_box is not None:
            # This ensures that we know a bounding box has been applied
            self._domain_override = True
            bbox = np.array(bounding_box, dtype="float64")
            if bbox.shape == (2, 3):
                bbox = bbox.transpose()
            self.domain_left_edge = bbox[:, 0]
            self.domain_right_edge = bbox[:, 1]
        else:
            self.domain_left_edge = self.domain_right_edge = None
        if units_override is not None:
            raise RuntimeError(
                "units_override is not supported for GadgetDataset. "
                + "Use unit_base instead."
            )

        # Set dark energy parameters before cosmology object is created
        self.use_dark_factor = use_dark_factor
        self.w_0 = w_0
        self.w_a = w_a

        super().__init__(
            filename,
            dataset_type=dataset_type,
            unit_system=unit_system,
            index_order=index_order,
            index_filename=index_filename,
            kdtree_filename=kdtree_filename,
            kernel_name=kernel_name,
            default_species_fields=default_species_fields,
        )
        if self.cosmological_simulation:
            self.time_unit.convert_to_units("s/h")
            self.length_unit.convert_to_units("kpccm/h")
            self.mass_unit.convert_to_units("g/h")
        else:
            self.time_unit.convert_to_units("s")
            self.length_unit.convert_to_units("kpc")
            self.mass_unit.convert_to_units("Msun")
        if mean_molecular_weight is None:
            self.mu = compute_mu(self.default_species_fields)
        else:
            self.mu = mean_molecular_weight

    @classmethod
    def _setup_binary_spec(cls, spec, spec_dict):
        if isinstance(spec, str):
            _hs = ()
            for hs in spec.split("+"):
                _hs += spec_dict[hs]
            spec = _hs
        return spec

    def __str__(self):
        return os.path.basename(self.parameter_filename).split(".")[0]

    def _get_hvals(self):
        self.gen_hsmls = False
        return self._header.value

    def _parse_parameter_file(self):

        hvals = self._get_hvals()

        self.dimensionality = 3
        self.refine_by = 2
        self.parameters["HydroMethod"] = "sph"
        # Set standard values

        # We may have an overridden bounding box.
        if self.domain_left_edge is None and hvals["BoxSize"] != 0:
            self.domain_left_edge = np.zeros(3, "float64")
            self.domain_right_edge = np.ones(3, "float64") * hvals["BoxSize"]

        self.domain_dimensions = np.ones(3, "int32")
        self._periodicity = (True, True, True)

        self.cosmological_simulation = 1

        try:
            self.current_redshift = hvals["Redshift"]
        except KeyError:
            # Probably not a cosmological dataset, we should just set
            # z = 0 and let the user know
            self.current_redshift = 0.0
            only_on_root(mylog.info, "Redshift is not set in Header. Assuming z=0.")

        try:
            self.omega_lambda = hvals["OmegaLambda"]
            self.omega_matter = hvals["Omega0"]
            self.hubble_constant = hvals["HubbleParam"]
        except KeyError:
            # If these are not set it is definitely not a cosmological dataset.
            self.omega_lambda = 0.0
            self.omega_matter = 1.0  # Just in case somebody asks for it.
            # Hubble is set below for Omega Lambda = 0.

        # According to the Gadget manual, OmegaLambda will be zero for
        # non-cosmological datasets.  However, it may be the case that
        # individuals are running cosmological simulations *without* Lambda, in
        # which case we may be doing something incorrect here.
        # It may be possible to deduce whether ComovingIntegration is on
        # somehow, but opinions on this vary.
        if self.omega_lambda == 0.0:
            only_on_root(
                mylog.info, "Omega Lambda is 0.0, so we are turning off Cosmology."
            )
            self.hubble_constant = 1.0  # So that scaling comes out correct
            self.cosmological_simulation = 0
            self.current_redshift = 0.0
            # This may not be correct.
            self.current_time = hvals["Time"]
        else:
            # Now we calculate our time based on the cosmology, because in
            # ComovingIntegration hvals["Time"] will in fact be the expansion
            # factor, not the actual integration time, so we re-calculate
            # global time from our Cosmology.
            cosmo = Cosmology(
                hubble_constant=self.hubble_constant,
                omega_matter=self.omega_matter,
                omega_lambda=self.omega_lambda,
            )
            self.current_time = cosmo.lookback_time(self.current_redshift, 1e6)
            only_on_root(
                mylog.info,
                "Calculating time from %0.3e to be %0.3e seconds",
                hvals["Time"],
                self.current_time,
            )
        self.parameters = hvals

        prefix = os.path.abspath(
            os.path.join(
                os.path.dirname(self.parameter_filename),
                os.path.basename(self.parameter_filename).split(".", 1)[0],
            )
        )

        if hvals["NumFiles"] > 1:
            for t in (
                f"{prefix}.%(num)s{self._suffix}",
                f"{prefix}.gad.%(num)s{self._suffix}",
            ):
                if os.path.isfile(t % {"num": 0}):
                    self.filename_template = t
                    break
            else:
                raise RuntimeError("Could not determine correct data file template.")
        else:
            self.filename_template = self.parameter_filename

        self.file_count = hvals["NumFiles"]

    def _set_code_unit_attributes(self):
        # If no units passed in by user, set a sane default (Gadget-2 users
        # guide).
        if self._unit_base is None:
            if self.cosmological_simulation == 1:
                only_on_root(
                    mylog.info, "Assuming length units are in kpc/h (comoving)"
                )
                self._unit_base = dict(length=(1.0, "kpccm/h"))
            else:
                only_on_root(mylog.info, "Assuming length units are in kpc (physical)")
                self._unit_base = dict(length=(1.0, "kpc"))

        # If units passed in by user, decide what to do about
        # co-moving and factors of h
        unit_base = self._unit_base or {}
        if "length" in unit_base:
            length_unit = unit_base["length"]
        elif "UnitLength_in_cm" in unit_base:
            if self.cosmological_simulation == 0:
                length_unit = (unit_base["UnitLength_in_cm"], "cm")
            else:
                length_unit = (unit_base["UnitLength_in_cm"], "cmcm/h")
        else:
            raise RuntimeError
        length_unit = _fix_unit_ordering(length_unit)
        self.length_unit = self.quan(length_unit[0], length_unit[1])

        unit_base = self._unit_base or {}

        if self.cosmological_simulation:
            # see http://www.mpa-garching.mpg.de/gadget/gadget-list/0113.html
            # for why we need to include a factor of square root of the
            # scale factor
            vel_units = "cm/s * sqrt(a)"
        else:
            vel_units = "cm/s"

        if "velocity" in unit_base:
            velocity_unit = unit_base["velocity"]
        elif "UnitVelocity_in_cm_per_s" in unit_base:
            velocity_unit = (unit_base["UnitVelocity_in_cm_per_s"], vel_units)
        else:
            velocity_unit = (1e5, vel_units)
        velocity_unit = _fix_unit_ordering(velocity_unit)
        self.velocity_unit = self.quan(velocity_unit[0], velocity_unit[1])

        # We set hubble_constant = 1.0 for non-cosmology, so this is safe.
        # Default to 1e10 Msun/h if mass is not specified.
        if "mass" in unit_base:
            mass_unit = unit_base["mass"]
        elif "UnitMass_in_g" in unit_base:
            if self.cosmological_simulation == 0:
                mass_unit = (unit_base["UnitMass_in_g"], "g")
            else:
                mass_unit = (unit_base["UnitMass_in_g"], "g/h")
        else:
            # Sane default
            mass_unit = (1e10, "Msun/h")
        mass_unit = _fix_unit_ordering(mass_unit)
        self.mass_unit = self.quan(mass_unit[0], mass_unit[1])
        if self.cosmological_simulation:
            # self.velocity_unit is the unit to rescale on-disk velocities, The
            # actual internal velocity unit is really in comoving units
            # since the time unit is derived from the internal velocity unit, we
            # infer the internal velocity unit here and name it vel_unit
            #
            # see http://www.mpa-garching.mpg.de/gadget/gadget-list/0113.html
            if "velocity" in unit_base:
                vel_unit = unit_base["velocity"]
            elif "UnitVelocity_in_cm_per_s" in unit_base:
                vel_unit = (unit_base["UnitVelocity_in_cm_per_s"], "cmcm/s")
            else:
                vel_unit = (1, "kmcm/s")
            vel_unit = self.quan(*vel_unit)
        else:
            vel_unit = self.velocity_unit
        self.time_unit = self.length_unit / vel_unit

        if "specific_energy" in unit_base:
            specific_energy_unit = unit_base["specific_energy"]
        elif "UnitEnergy_in_cgs" in unit_base and "UnitMass_in_g" in unit_base:
            specific_energy_unit = (
                unit_base["UnitEnergy_in_cgs"] / unit_base["UnitMass_in_g"]
            )
            specific_energy_unit = (specific_energy_unit, "(cm/s)**2")
        else:
            # Sane default
            specific_energy_unit = (1, "(km/s) ** 2")
        specific_energy_unit = _fix_unit_ordering(specific_energy_unit)
        self.specific_energy_unit = self.quan(*specific_energy_unit)

        if "magnetic" in unit_base:
            magnetic_unit = unit_base["magnetic"]
        elif "UnitMagneticField_in_gauss" in unit_base:
            magnetic_unit = (unit_base["UnitMagneticField_in_gauss"], "gauss")
        else:
            # Sane default
            magnetic_unit = (1.0, "gauss")
        magnetic_unit = _fix_unit_ordering(magnetic_unit)
        self.magnetic_unit = self.quan(*magnetic_unit)

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        if "header_spec" in kwargs:
            header_spec = kwargs["header_spec"]
        else:
            header_spec = "default"
        # Check to see if passed filename is a directory. If so, use it to get
        # the .0 snapshot file. Make sure there's only one such file, otherwise
        # there's an ambiguity about which file the user wants. Ignore ewah files
        if os.path.isdir(filename):
            valid_files = []
            for f in os.listdir(filename):
                fname = os.path.join(filename, f)
                if (".0" in f) and (".ewah" not in f) and os.path.isfile(fname):
                    valid_files.append(f)
            if len(valid_files) == 0:
                return False
            elif len(valid_files) > 1:
                return False
            else:
                validated_file = os.path.join(filename, valid_files[0])
        else:
            validated_file = filename
        header = GadgetBinaryHeader(validated_file, header_spec)
        return header.validate()


class GadgetHDF5File(ParticleFile):
    pass


class GadgetHDF5Dataset(GadgetDataset):
    _file_class = GadgetHDF5File
    _index_class = SPHParticleIndex
    _field_info_class: Type[FieldInfoContainer] = GadgetFieldInfo
    _particle_mass_name = "Masses"
    _sph_ptypes = ("PartType0",)
    _suffix = ".hdf5"

    def __init__(
        self,
        filename,
        dataset_type="gadget_hdf5",
        unit_base=None,
        index_order=None,
        index_filename=None,
        kernel_name=None,
        bounding_box=None,
        units_override=None,
        unit_system="cgs",
        default_species_fields=None,
    ):
        self.storage_filename = None
        filename = os.path.abspath(filename)
        if units_override is not None:
            raise RuntimeError(
                "units_override is not supported for GadgetHDF5Dataset. "
                "Use unit_base instead."
            )
        super().__init__(
            filename,
            dataset_type,
            unit_base=unit_base,
            index_order=index_order,
            index_filename=index_filename,
            kernel_name=kernel_name,
            bounding_box=bounding_box,
            unit_system=unit_system,
            default_species_fields=default_species_fields,
        )

    def _get_hvals(self):
        handle = h5py.File(self.parameter_filename, mode="r")
        hvals = {}
        hvals.update((str(k), v) for k, v in handle["/Header"].attrs.items())
        # Compat reasons.
        hvals["NumFiles"] = hvals["NumFilesPerSnapshot"]
        hvals["Massarr"] = hvals["MassTable"]
        sph_ptypes = [ptype for ptype in self._sph_ptypes if ptype in handle]
        if sph_ptypes:
            self.gen_hsmls = "SmoothingLength" not in handle[sph_ptypes[0]]
        else:
            self.gen_hsmls = False
        # Later versions of Gadget and its derivatives have a "Parameters"
        # group in the HDF5 file.
        if "Parameters" in handle:
            hvals.update((str(k), v) for k, v in handle["/Parameters"].attrs.items())
        handle.close()
        return hvals

    def _get_uvals(self):
        handle = h5py.File(self.parameter_filename, mode="r")
        uvals = {}
        uvals.update((str(k), v) for k, v in handle["/Units"].attrs.items())
        handle.close()
        return uvals

    def _set_owls_eagle(self):

        self.dimensionality = 3
        self.refine_by = 2
        self.parameters["HydroMethod"] = "sph"

        self._unit_base = self._get_uvals()
        self._unit_base["cmcm"] = 1.0 / self._unit_base["UnitLength_in_cm"]

        self.current_redshift = self.parameters["Redshift"]
        self.omega_lambda = self.parameters["OmegaLambda"]
        self.omega_matter = self.parameters["Omega0"]
        self.hubble_constant = self.parameters["HubbleParam"]

        if self.domain_left_edge is None and self.parameters["BoxSize"] != 0:
            self.domain_left_edge = np.zeros(3, "float64")
            self.domain_right_edge = np.ones(3, "float64") * self.parameters["BoxSize"]

        self.domain_dimensions = np.ones(3, "int32")

        self.cosmological_simulation = 1
        self._periodicity = (True, True, True)

        prefix = os.path.abspath(
            os.path.join(
                os.path.dirname(self.parameter_filename),
                os.path.basename(self.parameter_filename).split(".", 1)[0],
            )
        )

        suffix = self.parameter_filename.rsplit(".", 1)[-1]
        if self.parameters["NumFiles"] > 1:
            self.filename_template = f"{prefix}.%(num)i.{suffix}"
        else:
            self.filename_template = self.parameter_filename

        self.file_count = self.parameters["NumFilesPerSnapshot"]

    def _set_owls_eagle_units(self):

        # note the contents of the HDF5 Units group are in _unit_base
        # note the velocity stored on disk is sqrt(a) dx/dt
        # physical velocity [cm/s] = a dx/dt
        # = sqrt(a) * velocity_on_disk * UnitVelocity_in_cm_per_s
        self.length_unit = self.quan(self._unit_base["UnitLength_in_cm"], "cmcm/h")
        self.mass_unit = self.quan(self._unit_base["UnitMass_in_g"], "g/h")
        self.velocity_unit = self.quan(
            self._unit_base["UnitVelocity_in_cm_per_s"], "cm/s * sqrt(a)"
        )
        self.time_unit = self.quan(self._unit_base["UnitTime_in_s"], "s/h")

        specific_energy_unit_cgs = (
            self._unit_base["UnitEnergy_in_cgs"] / self._unit_base["UnitMass_in_g"]
        )
        self.specific_energy_unit = self.quan(specific_energy_unit_cgs, "(cm/s)**2")

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        need_groups = ["Header"]
        veto_groups = ["FOF", "Group", "Subhalo"]
        valid = True
        try:
            fh = h5py.File(filename, mode="r")
            valid = all(ng in fh["/"] for ng in need_groups) and not any(
                vg in fh["/"] for vg in veto_groups
            )
            fh.close()
        except Exception:
            valid = False
            pass

        try:
            fh = h5py.File(filename, mode="r")
            valid = fh["Header"].attrs["Code"].decode("utf-8") != "SWIFT"
            fh.close()
        except (OSError, KeyError, ImportError):
            pass

        return valid
