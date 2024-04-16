import abc
import os
import struct
from itertools import chain, count
from typing import TYPE_CHECKING, Any, Callable, Optional

import numpy as np

from yt._typing import FieldKey
from yt.config import ytcfg
from yt.funcs import mylog
from yt.utilities.cython_fortran_utils import FortranFile

from .field_handlers import HandlerMixin
from .io import (
    _ramses_particle_binary_file_handler,
    _ramses_particle_csv_file_handler,
    _read_part_binary_file_descriptor,
    _read_part_csv_file_descriptor,
)

if TYPE_CHECKING:
    from yt.frontends.ramses.data_structures import RAMSESDomainSubset

PARTICLE_HANDLERS: set[type["ParticleFileHandler"]] = set()


def get_particle_handlers():
    return PARTICLE_HANDLERS


def register_particle_handler(ph):
    PARTICLE_HANDLERS.add(ph)


class ParticleFileHandler(abc.ABC, HandlerMixin):
    """
    Abstract class to handle particles in RAMSES. Each instance
    represents a single file (one domain).

    To add support to a new particle file, inherit from this class and
    implement all functions containing a `NotImplementedError`.

    See `SinkParticleFileHandler` for an example implementation."""

    _file_type = "particle"

    ## These properties are static
    # The name to give to the particle type
    ptype: str

    # The name of the file(s).
    fname: str

    # The name of the file descriptor (if any)
    file_descriptor: Optional[str] = None

    # The attributes of the header
    attrs: tuple[tuple[str, int, str], ...]

    # A list of tuple containing the field name and its type
    known_fields: list[FieldKey]

    # The function to employ to read the file
    reader: Callable[
        ["ParticleFileHandler", "RAMSESDomainSubset", list[tuple[str, str]], int],
        dict[tuple[str, str], np.ndarray],
    ]

    # Name of the config section (if any)
    config_field: Optional[str] = None

    ## These properties are computed dynamically
    # Mapping from field to offset in file
    _field_offsets: dict[tuple[str, str], int]

    # Mapping from field to the type of the data (float, integer, ...)
    _field_types: dict[tuple[str, str], str]

    # Number of particle in the domain
    _local_particle_count: int

    # The header of the file
    _header: dict[str, Any]

    def __init_subclass__(cls, *args, **kwargs):
        """
        Registers subclasses at creation.
        """
        super().__init_subclass__(*args, **kwargs)

        if cls.ptype is not None:
            register_particle_handler(cls)

        cls._unique_registry = {}
        return cls

    def __init__(self, domain):
        self.setup_handler(domain)

        # Attempt to read the list of fields from the config file
        if self.config_field and ytcfg.has_section(self.config_field):
            cfg = ytcfg.get(self.config_field, "fields")
            known_fields = []
            for c in (_.strip() for _ in cfg.split("\n") if _.strip() != ""):
                field, field_type = (_.strip() for _ in c.split(","))
                known_fields.append((field, field_type))
            self.known_fields = known_fields

    @abc.abstractmethod
    def read_header(self):
        """
        This function is called once per file. It should:

        * read the header of the file and store any relevant information
        * detect the fields in the file
        * compute the offsets (location in the file) of each field

        It is in charge of setting `self.field_offsets` and `self.field_types`.

        * `_field_offsets`: dictionary: tuple -> integer
            A dictionary that maps `(type, field_name)` to their
            location in the file (integer)
        * `_field_types`: dictionary: tuple -> character
            A dictionary that maps `(type, field_name)` to their type
            (character), following Python's struct convention.
        """
        pass

    @property
    def field_offsets(self) -> dict[tuple[str, str], int]:
        if hasattr(self, "_field_offsets"):
            return self._field_offsets
        self.read_header()
        return self._field_offsets

    @property
    def field_types(self) -> dict[tuple[str, str], str]:
        if hasattr(self, "_field_types"):
            return self._field_types
        self.read_header()
        return self._field_types

    @property
    def local_particle_count(self) -> int:
        if hasattr(self, "_local_particle_count"):
            return self._local_particle_count
        self.read_header()
        return self._local_particle_count

    @property
    def header(self) -> dict[str, Any]:
        if hasattr(self, "_header"):
            return self._header
        self.read_header()
        return self._header


_default_dtypes: dict[int, str] = {
    1: "c",  # char
    2: "h",  # short,
    4: "f",  # float
    8: "d",  # double
}


class DefaultParticleFileHandler(ParticleFileHandler):
    ptype = "io"
    fname = "part_{iout:05d}.out{icpu:05d}"
    file_descriptor = "part_file_descriptor.txt"
    config_field = "ramses-particles"
    reader = _ramses_particle_binary_file_handler

    attrs = (
        ("ncpu", 1, "i"),
        ("ndim", 1, "i"),
        ("npart", 1, "i"),
        ("localseed", -1, "i"),
        ("nstar_tot", 1, "i"),
        ("mstar_tot", 1, "d"),
        ("mstar_lost", 1, "d"),
        ("nsink", 1, "i"),
    )

    known_fields = [
        ("particle_position_x", "d"),
        ("particle_position_y", "d"),
        ("particle_position_z", "d"),
        ("particle_velocity_x", "d"),
        ("particle_velocity_y", "d"),
        ("particle_velocity_z", "d"),
        ("particle_mass", "d"),
        ("particle_identity", "i"),
        ("particle_refinement_level", "i"),
    ]

    def read_header(self):
        if not self.exists:
            self._field_offsets = {}
            self._field_types = {}
            self._local_particle_count = 0
            self._header = {}
            return

        flen = os.path.getsize(self.fname)
        with FortranFile(self.fname) as fd:
            hvals = dict(fd.read_attrs(self.attrs))
            particle_field_pos = fd.tell()

        self._header = hvals
        self._local_particle_count = hvals["npart"]
        extra_particle_fields = self.ds._extra_particle_fields

        if self.has_descriptor:
            particle_fields = _read_part_binary_file_descriptor(self.file_descriptor)
        else:
            particle_fields = list(self.known_fields)

            if extra_particle_fields is not None:
                particle_fields += extra_particle_fields

        if (
            hvals["nstar_tot"] > 0
            and extra_particle_fields is not None
            and ("particle_birth_time", "d") not in particle_fields
        ):
            particle_fields += [
                ("particle_birth_time", "d"),
                ("particle_metallicity", "d"),
            ]

        def build_iterator():
            return chain(
                particle_fields,
                ((f"particle_extra_field_{i}", "d") for i in count(1)),
            )

        field_offsets = {}
        _pfields = {}
        ptype = self.ptype
        blockLen = struct.calcsize("i") * 2

        particle_fields_iterator = build_iterator()
        ipos = particle_field_pos
        while ipos < flen:
            field, vtype = next(particle_fields_iterator)
            field_offsets[ptype, field] = ipos
            _pfields[ptype, field] = vtype
            ipos += blockLen + struct.calcsize(vtype) * hvals["npart"]

        if ipos != flen:
            particle_fields_iterator = build_iterator()
            with FortranFile(self.fname) as fd:
                fd.seek(particle_field_pos)
                ipos = particle_field_pos
                while ipos < flen:
                    field, vtype = next(particle_fields_iterator)
                    old_pos = fd.tell()
                    field_offsets[ptype, field] = old_pos
                    fd.skip(1)
                    ipos = fd.tell()

                    record_len = ipos - old_pos - blockLen
                    exp_len = struct.calcsize(vtype) * hvals["npart"]

                    if record_len != exp_len:
                        # Guess record vtype from length
                        nbytes = record_len // hvals["npart"]
                        vtype = _default_dtypes[nbytes]

                        mylog.warning(
                            "Supposed that `%s` has type %s given record size",
                            field,
                            np.dtype(vtype),
                        )

                    _pfields[ptype, field] = vtype

        if field.startswith("particle_extra_field_"):
            iextra = int(field.split("_")[-1])
        else:
            iextra = 0
        if iextra > 0 and not self.ds._warned_extra_fields["io"]:
            mylog.warning(
                "Detected %s extra particle fields assuming kind "
                "`double`. Consider using the `extra_particle_fields` "
                "keyword argument if you have unexpected behavior.",
                iextra,
            )
            self.ds._warned_extra_fields["io"] = True

        self._field_offsets = field_offsets
        self._field_types = _pfields


class SinkParticleFileHandler(ParticleFileHandler):
    """Handle sink files"""

    ptype = "sink"
    fname = "sink_{iout:05d}.out{icpu:05d}"
    file_descriptor = "sink_file_descriptor.txt"
    config_field = "ramses-sink-particles"
    reader = _ramses_particle_binary_file_handler

    attrs = (("nsink", 1, "i"), ("nindsink", 1, "i"))

    known_fields = [
        ("particle_identifier", "i"),
        ("particle_mass", "d"),
        ("particle_position_x", "d"),
        ("particle_position_y", "d"),
        ("particle_position_z", "d"),
        ("particle_velocity_x", "d"),
        ("particle_velocity_y", "d"),
        ("particle_velocity_z", "d"),
        ("particle_birth_time", "d"),
        ("BH_real_accretion", "d"),
        ("BH_bondi_accretion", "d"),
        ("BH_eddington_accretion", "d"),
        ("BH_esave", "d"),
        ("gas_spin_x", "d"),
        ("gas_spin_y", "d"),
        ("gas_spin_z", "d"),
        ("BH_spin_x", "d"),
        ("BH_spin_y", "d"),
        ("BH_spin_z", "d"),
        ("BH_spin", "d"),
        ("BH_efficiency", "d"),
    ]

    def read_header(self):
        if not self.exists:
            self._field_offsets = {}
            self._field_types = {}
            self._local_particle_count = 0
            self._header = {}
            return
        fd = FortranFile(self.fname)
        flen = os.path.getsize(self.fname)
        hvals = {}
        # Read the header of the file
        attrs = self.attrs

        hvals.update(fd.read_attrs(attrs))
        self._header = hvals

        # This is somehow a trick here: we only want one domain to
        # be read, as ramses writes all the sinks in all the
        # domains. Here, we set the local_particle_count to 0 except
        # for the first domain to be red.
        if getattr(self.ds, "_sink_file_flag", False):
            self._local_particle_count = 0
        else:
            self.ds._sink_file_flag = True
            self._local_particle_count = hvals["nsink"]

        # Read the fields + add the sink properties
        if self.has_descriptor:
            fields = _read_part_binary_file_descriptor(self.file_descriptor)
        else:
            fields = list(self.known_fields)

        # Note: this follows RAMSES convention.
        for i in range(self.ds.dimensionality * 2 + 1):
            for ilvl in range(self.ds.max_level + 1):
                fields.append((f"particle_prop_{ilvl}_{i}", "d"))

        field_offsets = {}
        _pfields = {}

        # Fill the fields, offsets and types
        self.fields = []
        for field, vtype in fields:
            self.fields.append(field)
            if fd.tell() >= flen:
                break
            field_offsets[self.ptype, field] = fd.tell()
            _pfields[self.ptype, field] = vtype
            fd.skip(1)
        self._field_offsets = field_offsets
        self._field_types = _pfields
        fd.close()


class SinkParticleFileHandlerCsv(ParticleFileHandler):
    """Handle sink files from a csv file, the format from the sink particle in ramses"""

    ptype = "sink_csv"
    fname = "sink_{iout:05d}.csv"
    file_descriptor = None
    config_field = "ramses-sink-particles"
    reader = _ramses_particle_csv_file_handler
    attrs = (("nsink", 1, "i"), ("nindsink", 1, "i"))

    def read_header(self):
        if not self.exists:
            self._field_offsets = {}
            self._field_types = {}
            self._local_particle_count = 0
            self._header = {}
            return
        field_offsets = {}
        _pfields = {}

        fields, self._local_particle_count = _read_part_csv_file_descriptor(self.fname)

        for ind, field in enumerate(fields):
            field_offsets[self.ptype, field] = ind
            _pfields[self.ptype, field] = "d"

        self._field_offsets = field_offsets
        self._field_types = _pfields
