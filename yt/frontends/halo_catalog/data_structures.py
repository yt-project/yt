import glob
import sys
import weakref
from collections import defaultdict
from functools import partial

import numpy as np

from yt.data_objects.selection_objects.data_selection_objects import (
    YTSelectionContainer,
)
from yt.data_objects.static_output import (
    ParticleDataset,
    ParticleFile,
    validate_index_order,
)
from yt.frontends.ytdata.data_structures import SavedDataset
from yt.funcs import parse_h5_attr
from yt.geometry.particle_geometry_handler import ParticleIndex
from yt.utilities.on_demand_imports import _h5py as h5py

from .fields import YTHaloCatalogFieldInfo, YTHaloCatalogHaloFieldInfo

if sys.version_info >= (3, 8):
    from functools import cached_property
else:
    from yt._maintenance.backports import cached_property


class HaloCatalogFile(ParticleFile):
    """
    Base class for data files of halo catalog datasets.

    This is mainly here to correct for periodicity when
    reading particle positions.
    """

    def __init__(self, ds, io, filename, file_id, frange):
        super().__init__(ds, io, filename, file_id, frange)

    def _read_particle_positions(self, ptype, f=None):
        raise NotImplementedError

    def _get_particle_positions(self, ptype, f=None):
        pcount = self.total_particles[ptype]
        if pcount == 0:
            return None

        # Correct for periodicity.
        dle = self.ds.domain_left_edge.to("code_length").v
        dw = self.ds.domain_width.to("code_length").v
        pos = self._read_particle_positions(ptype, f=f)
        si, ei = self.start, self.end
        if None not in (si, ei):
            pos = pos[si:ei]

        np.subtract(pos, dle, out=pos)
        np.mod(pos, dw, out=pos)
        np.add(pos, dle, out=pos)

        return pos


class YTHaloCatalogFile(HaloCatalogFile):
    """
    Data file class for the YTHaloCatalogDataset.
    """

    def __init__(self, ds, io, filename, file_id, frange):
        with h5py.File(filename, mode="r") as f:
            self.header = {field: parse_h5_attr(f, field) for field in f.attrs.keys()}
            pids = f.get("particles/ids")
            self.total_ids = 0 if pids is None else pids.size
            self.group_length_sum = self.total_ids
        super().__init__(ds, io, filename, file_id, frange)

    def _read_particle_positions(self, ptype, f=None):
        """
        Read all particle positions in this file.
        """

        if f is None:
            close = True
            f = h5py.File(self.filename, mode="r")
        else:
            close = False

        pcount = self.header["num_halos"]
        pos = np.empty((pcount, 3), dtype="float64")
        for i, ax in enumerate("xyz"):
            pos[:, i] = f[f"particle_position_{ax}"][()]

        if close:
            f.close()

        return pos


class YTHaloCatalogDataset(SavedDataset):
    """
    Dataset class for halo catalogs made with yt.

    This covers yt FoF/HoP halo finders and the halo analysis
    in yt_astro_analysis.
    """

    _index_class = ParticleIndex
    _file_class = YTHaloCatalogFile
    _field_info_class = YTHaloCatalogFieldInfo
    _suffix = ".h5"
    _con_attrs = (
        "cosmological_simulation",
        "current_time",
        "current_redshift",
        "hubble_constant",
        "omega_matter",
        "omega_lambda",
        "domain_left_edge",
        "domain_right_edge",
    )

    def __init__(
        self,
        filename,
        dataset_type="ythalocatalog",
        index_order=None,
        units_override=None,
        unit_system="cgs",
    ):
        self.index_order = validate_index_order(index_order)
        super().__init__(
            filename,
            dataset_type,
            units_override=units_override,
            unit_system=unit_system,
        )

    def add_field(self, *args, **kwargs):
        super().add_field(*args, **kwargs)
        self._halos_ds.add_field(*args, **kwargs)

    @property
    def halos_field_list(self):
        return self._halos_ds.field_list

    @property
    def halos_derived_field_list(self):
        return self._halos_ds.derived_field_list

    @cached_property
    def _halos_ds(self):
        return YTHaloDataset(self)

    def _setup_classes(self):
        super()._setup_classes()
        self.halo = partial(YTHaloCatalogHaloContainer, ds=self._halos_ds)
        self.halo.__doc__ = YTHaloCatalogHaloContainer.__doc__

    def _parse_parameter_file(self):
        self.refine_by = 2
        self.dimensionality = 3
        self.domain_dimensions = np.ones(self.dimensionality, "int32")
        self._periodicity = (True, True, True)
        prefix = ".".join(self.parameter_filename.rsplit(".", 2)[:-2])
        self.filename_template = f"{prefix}.%(num)s{self._suffix}"
        self.file_count = len(glob.glob(prefix + "*" + self._suffix))
        self.particle_types = ("halos",)
        self.particle_types_raw = ("halos",)
        super()._parse_parameter_file()

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        if not filename.endswith(".h5"):
            return False
        with h5py.File(filename, mode="r") as f:
            if (
                "data_type" in f.attrs
                and parse_h5_attr(f, "data_type") == "halo_catalog"
            ):
                return True
        return False


class YTHaloParticleIndex(ParticleIndex):
    """
    Particle index for getting halo particles from YTHaloCatalogDatasets.
    """

    def __init__(self, ds, dataset_type):
        self.real_ds = weakref.proxy(ds.real_ds)
        super().__init__(ds, dataset_type)

    def _calculate_particle_index_starts(self):
        """
        Create a dict of halo id offsets for each file.
        """
        particle_count = defaultdict(int)
        offset_count = 0
        for data_file in self.data_files:
            data_file.index_start = {
                ptype: particle_count[ptype] for ptype in data_file.total_particles
            }
            data_file.offset_start = offset_count
            for ptype in data_file.total_particles:
                particle_count[ptype] += data_file.total_particles[ptype]
            offset_count += getattr(data_file, "total_offset", 0)

        self._halo_index_start = {}
        for ptype in self.ds.particle_types_raw:
            d = [data_file.index_start[ptype] for data_file in self.data_files]
            self._halo_index_start.update({ptype: np.array(d)})

    def _detect_output_fields(self):
        field_list = []
        scalar_field_list = []
        units = {}

        pc = {}
        for ptype in self.ds.particle_types_raw:
            d = [df.total_particles[ptype] for df in self.data_files]
            pc.update({ptype: sum(d)})
        found_fields = {ptype: False for ptype, pnum in pc.items() if pnum > 0}
        has_ids = False

        for data_file in self.data_files:
            fl, sl, idl, _units = self.io._identify_fields(data_file)
            units.update(_units)
            field_list.extend([f for f in fl if f not in field_list])
            scalar_field_list.extend([f for f in sl if f not in scalar_field_list])
            for ptype in found_fields:
                found_fields[ptype] |= data_file.total_particles[ptype]
            has_ids |= len(idl) > 0
            if all(found_fields.values()) and has_ids:
                break

        self.field_list = field_list
        self.scalar_field_list = scalar_field_list
        ds = self.dataset
        ds.scalar_field_list = scalar_field_list
        ds.particle_types = tuple({pt for pt, ds in field_list})
        ds.field_units.update(units)
        ds.particle_types_raw = ds.particle_types

    def _get_halo_file_indices(self, ptype, identifiers):
        """
        Get the index of the data file list where this halo lives.

        Digitize returns i such that bins[i-1] <= x < bins[i], so we subtract
        one because we will open data file i.
        """
        return np.digitize(identifiers, self._halo_index_start[ptype], right=False) - 1

    def _get_halo_scalar_index(self, ptype, identifier):
        i_scalar = self._get_halo_file_indices(ptype, [identifier])[0]
        scalar_index = identifier - self._halo_index_start[ptype][i_scalar]
        return scalar_index

    def _get_halo_values(self, ptype, identifiers, fields, f=None):
        """
        Get field values for halo data containers.
        """

        # if a file is already open, don't open it again
        filename = None if f is None else f.filename

        data = defaultdict(lambda: np.empty(identifiers.size))
        i_scalars = self._get_halo_file_indices(ptype, identifiers)
        for i_scalar in np.unique(i_scalars):
            # mask array to get field data for this halo
            target = i_scalars == i_scalar
            scalar_indices = identifiers - self._halo_index_start[ptype][i_scalar]

            # only open file if it's not already open
            my_f = (
                f
                if self.data_files[i_scalar].filename == filename
                else h5py.File(self.data_files[i_scalar].filename, mode="r")
            )

            for field in fields:
                data[field][target] = self._read_halo_particle_field(
                    my_f, ptype, field, scalar_indices[target]
                )

            if self.data_files[i_scalar].filename != filename:
                my_f.close()

        return data

    def _identify_base_chunk(self, dobj):
        pass

    def _read_halo_particle_field(self, fh, ptype, field, indices):
        return fh[field][indices]

    def _read_particle_fields(self, fields, dobj, chunk=None):
        if not fields:
            return {}, []
        fields_to_read, fields_to_generate = self._split_fields(fields)
        if not fields_to_read:
            return {}, fields_to_generate
        fields_to_return = self.io._read_particle_selection(dobj, fields_to_read)
        return fields_to_return, fields_to_generate

    def _setup_data_io(self):
        super()._setup_data_io()
        if self.real_ds._instantiated_index is None:
            self.real_ds.index
        self.real_ds.index

        # inherit some things from parent index
        self._data_files = self.real_ds.index.data_files
        self._total_particles = self.real_ds.index.total_particles

        self._calculate_particle_index_starts()


class HaloDataset(ParticleDataset):
    """
    Base class for dataset accessing particles from halo catalogs.
    """

    def __init__(self, ds, dataset_type):
        self.real_ds = ds
        for attr in [
            "filename_template",
            "file_count",
            "particle_types_raw",
            "particle_types",
            "_periodicity",
        ]:
            setattr(self, attr, getattr(self.real_ds, attr))

        super().__init__(self.real_ds.parameter_filename, dataset_type)

    def print_key_parameters(self):
        pass

    def _set_derived_attrs(self):
        pass

    def _parse_parameter_file(self):
        for attr in [
            "cosmological_simulation",
            "cosmology",
            "current_redshift",
            "current_time",
            "dimensionality",
            "domain_dimensions",
            "domain_left_edge",
            "domain_right_edge",
            "domain_width",
            "hubble_constant",
            "omega_lambda",
            "omega_matter",
            "unique_identifier",
        ]:
            setattr(self, attr, getattr(self.real_ds, attr))

    def set_code_units(self):
        self._set_code_unit_attributes()
        self.unit_registry = self.real_ds.unit_registry

    def _set_code_unit_attributes(self):
        for unit in ["length", "time", "mass", "velocity", "magnetic", "temperature"]:
            my_unit = f"{unit}_unit"
            setattr(self, my_unit, getattr(self.real_ds, my_unit, None))

    def __str__(self):
        return f"{self.real_ds}"

    def _setup_classes(self):
        self.objects = []


class YTHaloDataset(HaloDataset):
    """
    Dataset used for accessing member particles from YTHaloCatalogDatasets.
    """

    _index_class = YTHaloParticleIndex
    _file_class = YTHaloCatalogFile
    _field_info_class = YTHaloCatalogHaloFieldInfo

    def __init__(self, ds, dataset_type="ythalo"):
        super().__init__(ds, dataset_type)

    def _set_code_unit_attributes(self):
        pass

    @classmethod
    def _is_valid(self, *args, **kwargs):
        # We don't ever want this to be loaded by yt.load.
        return False


class HaloContainer(YTSelectionContainer):
    """
    Base class for data containers providing halo particles.
    """

    _type_name = "halo"
    _con_args = ("ptype", "particle_identifier")
    _skip_add = True
    _spatial = False

    def __init__(self, ptype, particle_identifier, ds=None):
        if ptype not in ds.particle_types_raw:
            raise RuntimeError(
                f'Possible halo types are {ds.particle_types_raw}, supplied "{ptype}".'
            )

        self.ptype = ptype
        self._current_particle_type = ptype
        super().__init__(ds, {})

        self._set_identifiers(particle_identifier)

        # Find the file that has the scalar values for this halo.
        i_scalar = self.index._get_halo_file_indices(ptype, [self.particle_identifier])[
            0
        ]
        self.i_scalar = i_scalar
        self.scalar_data_file = self.index.data_files[i_scalar]

        # Data files containing particles belonging to this halo.
        self.field_data_files = [self.index.data_files[i_scalar]]

        # index within halo arrays that corresponds to this halo
        self.scalar_index = self.index._get_halo_scalar_index(
            ptype, self.particle_identifier
        )

        self._set_io_data()
        self.particle_number = self._get_particle_number()

        # starting and ending indices for each file containing particles
        self._set_field_indices()

    @cached_property
    def mass(self):
        return self[self.ptype, "particle_mass"][0]

    @cached_property
    def radius(self):
        return self[self.ptype, "virial_radius"][0]

    @cached_property
    def position(self):
        return self[self.ptype, "particle_position"][0]

    @cached_property
    def velocity(self):
        return self[self.ptype, "particle_velocity"][0]

    def _set_io_data(self):
        halo_fields = self._get_member_fieldnames()
        my_data = self.index._get_halo_values(
            self.ptype, np.array([self.particle_identifier]), halo_fields
        )
        self._io_data = {field: np.int64(val[0]) for field, val in my_data.items()}

    def __repr__(self):
        return f"{self.ds}_{self.ptype}_{self.particle_identifier:09d}"


class YTHaloCatalogHaloContainer(HaloContainer):
    """
    Data container for accessing particles from a halo.

    Create a data container to get member particles and individual
    values from halos and subhalos. Halo mass, radius, position, and
    velocity are set as attributes. Halo IDs are accessible
    through the field, "member_ids".  Other fields that are one
    value per halo are accessible as normal.  The field list for
    halo objects can be seen in `ds.halos_field_list`.

    Parameters
    ----------
    ptype : string
        The type of halo. Possible options can be found by
        inspecting the value of ds.particle_types_raw.
    particle_identifier : int
        The halo id.

    Examples
    --------

    >>> import yt
    >>> ds = yt.load("tiny_fof_halos/DD0046/DD0046.0.h5")

    >>> halo = ds.halo("halos", 0)
    >>> print(halo.particle_identifier)
    0
    >>> print(halo.mass)
    8724990744704.453 Msun
    >>> print(halo.radius)
    658.8140635766607 kpc
    >>> print(halo.position)
    [0.05496909 0.19451951 0.04056824] code_length
    >>> print(halo.velocity)
    [7034181.07118151 5323471.09102874 3234522.50495914] cm/s
    >>> # particle ids for this halo
    >>> print(halo["member_ids"])
    [ 1248.   129.   128. 31999. 31969. 31933. 31934.   159. 31903. 31841. ...
      2241.  2240.  2239.  2177.  2209.  2207.  2208.] dimensionless

    """

    def _get_member_fieldnames(self):
        return ["particle_number", "particle_index_start"]

    def _get_particle_number(self):
        return self._io_data["particle_number"]

    def _set_field_indices(self):
        self.field_data_start = [self._io_data["particle_index_start"]]
        self.field_data_end = [self.field_data_start[0] + self.particle_number]

    def _set_identifiers(self, particle_identifier):
        self.particle_identifier = particle_identifier
        self.group_identifier = self.particle_identifier
