import abc
import glob
import os
from typing import List, Optional, Set, Tuple, Type

from yt.config import ytcfg
from yt.funcs import mylog
from yt.utilities.cython_fortran_utils import FortranFile

from .io import _read_fluid_file_descriptor
from .io_utils import read_offset

FIELD_HANDLERS: Set[Type["FieldFileHandler"]] = set()


def get_field_handlers():
    return FIELD_HANDLERS


def register_field_handler(ph):
    FIELD_HANDLERS.add(ph)


DETECTED_FIELDS = {}  # type: ignore


class HandlerMixin:
    """This contains all the shared methods to handle RAMSES files.

    This is not supposed to be user-facing.
    """

    def setup_handler(self, domain):
        """
        Initialize an instance of the class. This automatically sets
        the full path to the file. This is not intended to be
        overridden in most cases.

        If you need more flexibility, rewrite this function to your
        need in the inherited class.
        """
        self.ds = ds = domain.ds
        self.domain = domain
        self.domain_id = domain.domain_id
        basename = os.path.abspath(ds.root_folder)
        iout = int(os.path.basename(ds.parameter_filename).split(".")[0].split("_")[1])

        if ds.num_groups > 0:
            igroup = ((domain.domain_id - 1) // ds.group_size) + 1
            full_path = os.path.join(
                basename,
                f"group_{igroup:05d}",
                self.fname.format(iout=iout, icpu=domain.domain_id),
            )
        else:
            full_path = os.path.join(
                basename, self.fname.format(iout=iout, icpu=domain.domain_id)
            )

        if os.path.exists(full_path):
            self.fname = full_path
        else:
            raise FileNotFoundError(
                f"Could not find {self._file_type} file (type: {self.ftype}). "
                f"Tried {full_path}"
            )

        if self.file_descriptor is not None:
            if ds.num_groups > 0:
                # The particle file descriptor is *only* in the first group
                self.file_descriptor = os.path.join(
                    basename, "group_00001", self.file_descriptor
                )
            else:
                self.file_descriptor = os.path.join(basename, self.file_descriptor)

    @property
    def exists(self):
        """
        This function should return True if the *file* the instance
        exists. It is called for each file of the type found on the
        disk.

        By default, it just returns whether the file exists. Override
        it for more complex cases.
        """
        return os.path.exists(self.fname)

    @property
    def has_descriptor(self):
        """
        This function should return True if a *file descriptor*
        exists.

        By default, it just returns whether the file exists. Override
        it for more complex cases.
        """
        return os.path.exists(self.file_descriptor)

    @classmethod
    def any_exist(cls, ds):
        """
        This function should return True if the kind of particle
        represented by the class exists in the dataset. It takes as
        argument the class itself -not an instance- and a dataset.

        Arguments
        ---------
        ds : a Ramses Dataset

        Note
        ----
        This function is usually called once at the initialization of
        the RAMSES Dataset structure to determine if the particle type
        (e.g. regular particles) exists.
        """
        if ds.unique_identifier in cls._unique_registry:
            return cls._unique_registry[ds.unique_identifier]

        iout = int(os.path.basename(ds.parameter_filename).split(".")[0].split("_")[1])

        fname = os.path.join(
            os.path.split(ds.parameter_filename)[0], cls.fname.format(iout=iout, icpu=1)
        )
        exists = os.path.exists(fname)
        cls._unique_registry[ds.unique_identifier] = exists

        return exists


class FieldFileHandler(abc.ABC, HandlerMixin):
    """
    Abstract class to handle particles in RAMSES. Each instance
    represents a single file (one domain).

    To add support to a new particle file, inherit from this class and
    implement all functions containing a `NotImplementedError`.

    See `SinkParticleFileHandler` for an example implementation."""

    _file_type = "field"

    # These properties are static properties
    ftype: Optional[str] = None  # The name to give to the field type
    fname: Optional[str] = None  # The name of the file(s)
    attrs: Optional[
        Tuple[Tuple[str, int, str], ...]
    ] = None  # The attributes of the header
    known_fields = None  # A list of tuple containing the field name and its type
    config_field: Optional[str] = None  # Name of the config section (if any)

    file_descriptor: Optional[str] = None  # The name of the file descriptor (if any)

    # These properties are computed dynamically
    field_offsets = None  # Mapping from field to offset in file
    field_types = (
        None  # Mapping from field to the type of the data (float, integer, ...)
    )

    def __init_subclass__(cls, *args, **kwargs):
        """
        Registers subclasses at creation.
        """
        super().__init_subclass__(*args, **kwargs)

        if cls.ftype is not None:
            register_field_handler(cls)

        cls._unique_registry = {}
        return cls

    def __init__(self, domain):
        self.setup_handler(domain)

    @classmethod
    @abc.abstractmethod
    def detect_fields(cls, ds):
        """
        Called once to setup the fields of this type

        It should set the following static variables:

        * parameters: dictionary
            Dictionary containing the variables. The keys should match
            those of `cls.attrs`
        * field_list: list of (ftype, fname)
            The list of the field present in the file
        """
        pass

    @classmethod
    def get_detected_fields(cls, ds):
        """
        Get the detected fields from the registry.
        """
        if ds.unique_identifier in DETECTED_FIELDS:
            d = DETECTED_FIELDS[ds.unique_identifier]
            if cls.ftype in d:
                return d[cls.ftype]

        return None

    @classmethod
    def set_detected_fields(cls, ds, fields):
        """
        Store the detected fields into the registry.
        """
        if ds.unique_identifier not in DETECTED_FIELDS:
            DETECTED_FIELDS[ds.unique_identifier] = {}

        DETECTED_FIELDS[ds.unique_identifier].update({cls.ftype: fields})

    @classmethod
    def purge_detected_fields(cls, ds):
        """
        Purge the registry.

        This should be called on dataset creation to force the field
        detection to be called.
        """
        if ds.unique_identifier in DETECTED_FIELDS:
            DETECTED_FIELDS.pop(ds.unique_identifier)

    @property
    def level_count(self):
        """
        Return the number of cells per level.
        """
        if getattr(self, "_level_count", None) is not None:
            return self._level_count
        self.offset

        return self._level_count

    @property
    def offset(self):
        """
        Compute the offsets of the fields.

        By default, it skips the header (as defined by `cls.attrs`)
        and computes the offset at each level.

        It should be generic enough for most of the cases, but if the
        *structure* of your fluid file is non-canonical, change this.
        """

        if getattr(self, "_offset", None) is not None:
            return self._offset

        nvars = len(self.field_list)
        with FortranFile(self.fname) as fd:
            # Skip headers
            nskip = len(self.attrs)
            fd.skip(nskip)
            min_level = self.domain.ds.min_level

            # The file is as follows:
            # > headers
            # loop over levels
            #   loop over cpu domains
            #     > <ilevel>: current level
            #     > <nocts>: number of octs in level, domain
            #     loop over <nvars> variables (positions, velocities, density, ...)
            #       loop over <2*2*2> cells in each oct
            #          > <data> with shape (nocts, )
            #
            # So there are 8 * nvars records each with length (nocts, )
            # at each (level, cpus)
            offset, level_count = read_offset(
                fd,
                min_level,
                self.domain.domain_id,
                self.parameters["nvar"],
                self.domain.amr_header,
                skip_len=nvars * 8,
            )

        self._offset = offset
        self._level_count = level_count
        return self._offset

    @classmethod
    def load_fields_from_yt_config(cls) -> List[str]:
        if cls.config_field and ytcfg.has_section(cls.config_field):
            cfg = ytcfg.get(cls.config_field, "fields")
            fields = [_.strip() for _ in cfg if _.strip() != ""]
            return fields

        return []


class HydroFieldFileHandler(FieldFileHandler):
    ftype = "ramses"
    fname = "hydro_{iout:05d}.out{icpu:05d}"
    file_descriptor = "hydro_file_descriptor.txt"
    config_field = "ramses-hydro"

    attrs = (
        ("ncpu", 1, "i"),
        ("nvar", 1, "i"),
        ("ndim", 1, "i"),
        ("nlevelmax", 1, "i"),
        ("nboundary", 1, "i"),
        ("gamma", 1, "d"),
    )

    @classmethod
    def detect_fields(cls, ds):
        # Try to get the detected fields
        detected_fields = cls.get_detected_fields(ds)
        if detected_fields:
            return detected_fields

        num = os.path.basename(ds.parameter_filename).split(".")[0].split("_")[1]
        testdomain = 1  # Just pick the first domain file to read
        basepath = os.path.abspath(os.path.dirname(ds.parameter_filename))
        basename = "%s/%%s_%s.out%05i" % (basepath, num, testdomain)
        fname = basename % "hydro"
        fname_desc = os.path.join(basepath, cls.file_descriptor)

        attrs = cls.attrs
        with FortranFile(fname) as fd:
            hvals = fd.read_attrs(attrs)
        cls.parameters = hvals

        # Store some metadata
        ds.gamma = hvals["gamma"]
        nvar = hvals["nvar"]

        ok = False

        if ds._fields_in_file is not None:
            # Case 1: fields are provided by users on construction of dataset
            fields = list(ds._fields_in_file)
            ok = True
        else:
            # Case 2: fields are provided by users in the config
            fields = cls.load_fields_from_yt_config()
            ok = len(fields) > 0

        if not ok and os.path.exists(fname_desc):
            # Case 3: there is a file descriptor
            # Or there is an hydro file descriptor
            mylog.debug("Reading hydro file descriptor.")
            # For now, we can only read double precision fields
            fields = [e[0] for e in _read_fluid_file_descriptor(fname_desc)]

            # We get no fields for old-style hydro file descriptor
            ok = len(fields) > 0

        if not ok:
            # Case 4: attempt autodetection with usual fields
            foldername = os.path.abspath(os.path.dirname(ds.parameter_filename))
            rt_flag = any(glob.glob(os.sep.join([foldername, "info_rt_*.txt"])))
            if rt_flag:  # rt run
                if nvar < 10:
                    mylog.info("Detected RAMSES-RT file WITHOUT IR trapping.")

                    fields = [
                        "Density",
                        "x-velocity",
                        "y-velocity",
                        "z-velocity",
                        "Pressure",
                        "Metallicity",
                        "HII",
                        "HeII",
                        "HeIII",
                    ]
                else:
                    mylog.info("Detected RAMSES-RT file WITH IR trapping.")

                    fields = [
                        "Density",
                        "x-velocity",
                        "y-velocity",
                        "z-velocity",
                        "Pres_IR",
                        "Pressure",
                        "Metallicity",
                        "HII",
                        "HeII",
                        "HeIII",
                    ]
            else:
                if nvar < 5:
                    mylog.debug(
                        "nvar=%s is too small! YT doesn't currently "
                        "support 1D/2D runs in RAMSES %s"
                    )
                    raise ValueError
                # Basic hydro runs
                if nvar == 5:
                    fields = [
                        "Density",
                        "x-velocity",
                        "y-velocity",
                        "z-velocity",
                        "Pressure",
                    ]
                if nvar > 5 and nvar < 11:
                    fields = [
                        "Density",
                        "x-velocity",
                        "y-velocity",
                        "z-velocity",
                        "Pressure",
                        "Metallicity",
                    ]
                # MHD runs - NOTE:
                # THE MHD MODULE WILL SILENTLY ADD 3 TO THE NVAR IN THE MAKEFILE
                if nvar == 11:
                    fields = [
                        "Density",
                        "x-velocity",
                        "y-velocity",
                        "z-velocity",
                        "B_x_left",
                        "B_y_left",
                        "B_z_left",
                        "B_x_right",
                        "B_y_right",
                        "B_z_right",
                        "Pressure",
                    ]
                if nvar > 11:
                    fields = [
                        "Density",
                        "x-velocity",
                        "y-velocity",
                        "z-velocity",
                        "B_x_left",
                        "B_y_left",
                        "B_z_left",
                        "B_x_right",
                        "B_y_right",
                        "B_z_right",
                        "Pressure",
                        "Metallicity",
                    ]
            mylog.debug(
                "No fields specified by user; automatically setting fields array to %s",
                fields,
            )

        # Allow some wiggle room for users to add too many variables
        count_extra = 0
        while len(fields) < nvar:
            fields.append(f"var_{len(fields)}")
            count_extra += 1
        if count_extra > 0:
            mylog.debug("Detected %s extra fluid fields.", count_extra)
        cls.field_list = [(cls.ftype, e) for e in fields]

        cls.set_detected_fields(ds, fields)

        return fields


class GravFieldFileHandler(FieldFileHandler):
    ftype = "gravity"
    fname = "grav_{iout:05d}.out{icpu:05d}"
    config_field = "ramses-grav"

    attrs = (
        ("ncpu", 1, "i"),
        ("nvar", 1, "i"),
        ("nlevelmax", 1, "i"),
        ("nboundary", 1, "i"),
    )

    @classmethod
    def detect_fields(cls, ds):
        ndim = ds.dimensionality
        iout = int(str(ds).split("_")[1])
        basedir = os.path.split(ds.parameter_filename)[0]
        fname = os.path.join(basedir, cls.fname.format(iout=iout, icpu=1))
        with FortranFile(fname) as fd:
            cls.parameters = fd.read_attrs(cls.attrs)

        nvar = cls.parameters["nvar"]
        ndim = ds.dimensionality

        fields = cls.load_fields_from_yt_config()

        if not fields:
            if nvar == ndim + 1:
                fields = ["Potential"] + [f"{k}-acceleration" for k in "xyz"[:ndim]]
            else:
                fields = [f"{k}-acceleration" for k in "xyz"[:ndim]]
            ndetected = len(fields)

            if ndetected != nvar and not ds._warned_extra_fields["gravity"]:
                mylog.info("Detected %s extra gravity fields.", nvar - ndetected)
                ds._warned_extra_fields["gravity"] = True

                for i in range(nvar - ndetected):
                    fields.append(f"var{i}")

        cls.field_list = [(cls.ftype, e) for e in fields]

        return fields


class RTFieldFileHandler(FieldFileHandler):
    ftype = "ramses-rt"
    fname = "rt_{iout:05d}.out{icpu:05d}"
    config_field = "ramses-rt"

    attrs = (
        ("ncpu", 1, "i"),
        ("nvar", 1, "i"),
        ("ndim", 1, "i"),
        ("nlevelmax", 1, "i"),
        ("nboundary", 1, "i"),
        ("gamma", 1, "d"),
    )

    @classmethod
    def detect_fields(cls, ds):
        # Try to get the detected fields
        detected_fields = cls.get_detected_fields(ds)
        if detected_fields:
            return detected_fields

        fname = ds.parameter_filename.replace("info_", "info_rt_")

        rheader = {}

        def read_rhs(cast):
            line = f.readline()
            p, v = line.split("=")
            rheader[p.strip()] = cast(v)

        with open(fname) as f:
            # Read nRTvar, nions, ngroups, iions
            for _ in range(4):
                read_rhs(int)

            # Try to read rtprecision.
            # Either it is present or the line is simply blank, so
            # we try to parse the line as an int, and if it fails,
            # we simply ignore it.
            try:
                read_rhs(int)
                f.readline()
            except ValueError:
                pass

            # Read X and Y fractions
            for _ in range(2):
                read_rhs(float)
            f.readline()

            # Reat unit_np, unit_pfd
            for _ in range(2):
                read_rhs(float)

            # Read rt_c_frac
            # Note: when using variable speed of light, this line will contain multiple
            # values corresponding the the velocity at each level
            read_rhs(lambda line: [float(e) for e in line.split()])
            f.readline()

            # Read n star, t2star, g_star
            for _ in range(3):
                read_rhs(float)

            # Touchy part, we have to read the photon group properties
            mylog.debug("Not reading photon group properties")

            cls.rt_parameters = rheader

        ngroups = rheader["nGroups"]

        iout = int(str(ds).split("_")[1])
        basedir = os.path.split(ds.parameter_filename)[0]
        fname = os.path.join(basedir, cls.fname.format(iout=iout, icpu=1))
        with FortranFile(fname) as fd:
            cls.parameters = fd.read_attrs(cls.attrs)

        fields = cls.load_fields_from_yt_config()

        if not fields:
            fields = []
            tmp = [
                "Photon_density_%s",
                "Photon_flux_x_%s",
                "Photon_flux_y_%s",
                "Photon_flux_z_%s",
            ]
            for ng in range(ngroups):
                fields.extend([t % (ng + 1) for t in tmp])

        cls.field_list = [(cls.ftype, e) for e in fields]

        cls.set_detected_fields(ds, fields)
        return fields

    @classmethod
    def get_rt_parameters(cls, ds):
        if cls.rt_parameters:
            return cls.rt_parameters

        # Call detect fields to get the rt_parameters
        cls.detect_fields(ds)
        return cls.rt_parameters
