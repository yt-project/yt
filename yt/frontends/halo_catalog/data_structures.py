from collections import defaultdict
from functools import partial
import glob
from yt.utilities.on_demand_imports import _h5py as h5py
import numpy as np
import os
import weakref

from .fields import \
    HaloCatalogFieldInfo, \
    HaloCatalogHaloFieldInfo

from yt.data_objects.data_containers import \
    YTSelectionContainer
from yt.data_objects.static_output import \
    ParticleDataset
from yt.frontends.ytdata.data_structures import \
    SavedDataset
from yt.funcs import \
    parse_h5_attr
from yt.geometry.particle_geometry_handler import \
    ParticleIndex
from yt.data_objects.static_output import \
    ParticleFile, \
    validate_index_order

class HaloCatalogFile(ParticleFile):
    def __init__(self, ds, io, filename, file_id, range):
        super(HaloCatalogFile, self).__init__(
            ds, io, filename, file_id, range)

    def _read_particle_positions(self, ptype, f=None):
        raise NotImplementedError

    def _get_particle_positions(self, ptype, f=None):
        pcount = self.total_particles[ptype]
        if pcount == 0:
            return None

        # Correct for periodicity.
        dle = self.ds.domain_left_edge.to('code_length').v
        dw = self.ds.domain_width.to('code_length').v
        pos = self._read_particle_positions(ptype, f=f)
        si, ei = self.start, self.end
        if None not in (si, ei):
            pos = pos[si:ei]

        np.subtract(pos, dle, out=pos)
        np.mod(pos, dw, out=pos)
        np.add(pos, dle, out=pos)

        return pos

class HaloCatalogHDF5File(HaloCatalogFile):
    def __init__(self, ds, io, filename, file_id, range):
        with h5py.File(filename, "r") as f:
            self.header = dict((field, parse_h5_attr(f, field)) \
                               for field in f.attrs.keys())
            pids = f.get('particles/ids')
            self.total_ids = 0 if pids is None else pids.size
            self.group_length_sum = self.total_ids
        super(HaloCatalogHDF5File, self).__init__(
            ds, io, filename, file_id, range)

    def _read_particle_positions(self, ptype, f=None):
        """
        Read all particle positions in this file.
        """

        if f is None:
            close = True
            f = h5py.File(self.filename, "r")
        else:
            close = False

        pcount = self.header["num_halos"]
        pos = np.empty((pcount, 3), dtype="float64")
        for i, ax in enumerate('xyz'):
            pos[:, i] = f["particle_position_%s" % ax][()]

        if close:
            f.close()

        return pos

class HaloCatalogDataset(SavedDataset):
    _index_class = ParticleIndex
    _file_class = HaloCatalogHDF5File
    _field_info_class = HaloCatalogFieldInfo
    _suffix = ".h5"
    _con_attrs = ("cosmological_simulation",
                  "current_time", "current_redshift",
                  "hubble_constant", "omega_matter", "omega_lambda",
                  "domain_left_edge", "domain_right_edge")

    def __init__(self, filename, dataset_type="halocatalog_hdf5",
                 index_order=None, units_override=None, unit_system="cgs"):
        self.index_order = validate_index_order(index_order)
        super(HaloCatalogDataset, self).__init__(filename, dataset_type,
                                                 units_override=units_override,
                                                 unit_system=unit_system)

    def add_field(self, *args, **kwargs):
        super(HaloCatalogDataset, self).add_field(*args, **kwargs)
        self._halos_ds.add_field(*args, **kwargs)

    @property
    def halos_field_list(self):
        return self._halos_ds.field_list

    @property
    def halos_derived_field_list(self):
        return self._halos_ds.derived_field_list

    _instantiated_halo_ds = None
    @property
    def _halos_ds(self):
        if self._instantiated_halo_ds is None:
            self._instantiated_halo_ds = HaloCatalogHaloDataset(self)
        return self._instantiated_halo_ds

    def _setup_classes(self):
        super(HaloCatalogDataset, self)._setup_classes()
        self.halo = partial(HaloCatalogHaloContainer, ds=self._halos_ds)

    def _parse_parameter_file(self):
        self.refine_by = 2
        self.dimensionality = 3
        self.domain_dimensions = np.ones(self.dimensionality, "int32")
        self.periodicity = (True, True, True)
        prefix = ".".join(self.parameter_filename.rsplit(".", 2)[:-2])
        self.filename_template = "%s.%%(num)s%s" % (prefix, self._suffix)
        self.file_count = len(glob.glob(prefix + "*" + self._suffix))
        self.particle_types = ("halos",)
        self.particle_types_raw = ("halos",)
        super(HaloCatalogDataset, self)._parse_parameter_file()

    @classmethod
    def _is_valid(self, *args, **kwargs):
        if not args[0].endswith(".h5"): return False
        with h5py.File(args[0], "r") as f:
            if "data_type" in f.attrs and \
              parse_h5_attr(f, "data_type") == "halo_catalog":
                return True
        return False

class HaloCatalogHaloParticleIndex(ParticleIndex):
    def __init__(self, ds, dataset_type):
        self.real_ds = weakref.proxy(ds.real_ds)
        super(HaloCatalogHaloParticleIndex, self).__init__(ds, dataset_type)

    def _setup_data_io(self):
        super(HaloCatalogHaloParticleIndex, self)._setup_data_io()
        self._setup_filenames()

    def _setup_geometry(self):
        self._setup_data_io()

        if self.real_ds._instantiated_index is None:
            template = self.real_ds.filename_template
            ndoms = self.real_ds.file_count
            cls = self.real_ds._file_class
            self.data_files = \
              [cls(self.dataset, self.io, template % {'num':i}, i, None)
               for i in range(ndoms)]
        else:
            self.data_files = self.real_ds.index.data_files

        self._calculate_particle_index_starts()
        self._calculate_particle_count()
        self._create_halo_id_table()

    def _calculate_particle_count(self):
        """
        Calculate the total number of each type of particle.
        """
        self.particle_count = \
          dict([(ptype, sum([d.total_particles[ptype] for d in self.data_files]))
                 for ptype in self.ds.particle_types_raw])

    def _calculate_particle_index_starts(self):
        """
        Create a dict of halo id offsets for each file.
        """
        particle_count = defaultdict(int)
        for data_file in self.data_files:
            data_file.index_start = dict([(ptype, particle_count[ptype]) for
                                           ptype in data_file.total_particles])
            for ptype in data_file.total_particles:
                particle_count[ptype] += data_file.total_particles[ptype]

        self._halo_index_start = \
          dict([(ptype, np.array([data_file.index_start[ptype]
                                  for data_file in self.data_files]))
                for ptype in self.ds.particle_types_raw])

    def _create_halo_id_table(self):
        """
        Create a list of halo start ids so we know which file
        contains particles for a given halo.  Note, the halo ids
        are distributed over all files and so the ids for a given
        halo are likely stored in a different file than the halo
        itself.
        """

        self._halo_id_number = np.array([data_file.total_ids
                            for data_file in self.data_files])
        self._halo_id_end = self._halo_id_number.cumsum()
        self._halo_id_start = self._halo_id_end - self._halo_id_number

        self._group_length_sum = \
          np.array([data_file.group_length_sum
                    for data_file in self.data_files])

    def _detect_output_fields(self):
        field_list = []
        scalar_field_list = []
        units = {}
        found_fields = \
          dict([(ptype, False)
                for ptype, pnum in self.particle_count.items()
                if pnum > 0])
        has_ids = False

        for data_file in self.data_files:
            fl, sl, idl, _units = self.io._identify_fields(data_file)
            units.update(_units)
            field_list.extend([f for f in fl
                               if f not in field_list])
            scalar_field_list.extend([f for f in sl
                                      if f not in scalar_field_list])
            for ptype in found_fields:
                found_fields[ptype] |= data_file.total_particles[ptype]
            has_ids |= len(idl) > 0
            if all(found_fields.values()) and has_ids: break

        self.field_list = field_list
        self.scalar_field_list = scalar_field_list
        ds = self.dataset
        ds.scalar_field_list = scalar_field_list
        ds.particle_types = tuple(set(pt for pt, ds in field_list))
        ds.field_units.update(units)
        ds.particle_types_raw = ds.particle_types

    def _identify_base_chunk(self, dobj):
        pass

    def _read_particle_fields(self, fields, dobj, chunk = None):
        if len(fields) == 0: return {}, []
        fields_to_read, fields_to_generate = self._split_fields(fields)
        if len(fields_to_read) == 0:
            return {}, fields_to_generate
        fields_to_return = self.io._read_particle_selection(
            dobj, fields_to_read)
        return fields_to_return, fields_to_generate

    def _get_halo_file_indices(self, ptype, identifiers):
        return np.digitize(identifiers,
            self._halo_index_start[ptype], right=False) - 1

    def _get_halo_scalar_index(self, ptype, identifier):
        i_scalar = self._get_halo_file_indices(ptype, [identifier])[0]
        scalar_index = identifier - self._halo_index_start[ptype][i_scalar]
        return scalar_index

    def _get_halo_values(self, ptype, identifiers, fields,
                         f=None):
        """
        Get field values for halos.  IDs are likely to be
        sequential (or at least monotonic), but not necessarily
        all within the same file.

        This does not do much to minimize file i/o, but with
        halos randomly distributed across files, there's not
        much more we can do.
        """

        # if a file is already open, don't open it again
        filename = None if f is None \
          else f.filename

        data = defaultdict(lambda: np.empty(identifiers.size))
        i_scalars = self._get_halo_file_indices(ptype, identifiers)
        for i_scalar in np.unique(i_scalars):
            target = i_scalars == i_scalar
            scalar_indices = identifiers - \
              self._halo_index_start[ptype][i_scalar]

            # only open file if it's not already open
            my_f = f if self.data_files[i_scalar].filename == filename \
              else h5py.File(self.data_files[i_scalar].filename, "r")

            for field in fields:
                data[field][target] = \
                  my_f[os.path.join(ptype, field)][()][scalar_indices[target]]

            if self.data_files[i_scalar].filename != filename: my_f.close()

        return data

class HaloCatalogHaloDataset(ParticleDataset):
    _index_class = HaloCatalogHaloParticleIndex
    _file_class = HaloCatalogHDF5File
    _field_info_class = HaloCatalogHaloFieldInfo

    def __init__(self, ds, dataset_type="halo_catalog_halo_hdf5"):
        self.real_ds = ds
        for attr in ['filename_template', 'file_count',
                     'particle_types_raw', 'particle_types',
                     'periodicity']:
            setattr(self, attr, getattr(self.real_ds, attr))

        super(HaloCatalogHaloDataset, self).__init__(
            self.real_ds.parameter_filename, dataset_type)

    def print_key_parameters(self):
        pass

    def _set_derived_attrs(self):
        pass

    def _parse_parameter_file(self):
        for attr in ["cosmological_simulation", "cosmology",
                     "current_redshift", "current_time",
                     "dimensionality", "domain_dimensions",
                     "domain_left_edge", "domain_right_edge",
                     "domain_width", "hubble_constant",
                     "omega_lambda", "omega_matter",
                     "unique_identifier"]:
            setattr(self, attr, getattr(self.real_ds, attr))

    def set_code_units(self):
        for unit in ["length", "time", "mass",
                     "velocity", "magnetic", "temperature"]:
            my_unit = "%s_unit" % unit
            setattr(self, my_unit, getattr(self.real_ds, my_unit, None))
        self.unit_registry = self.real_ds.unit_registry

    def __repr__(self):
        return "%s" % self.real_ds

    def _setup_classes(self):
        self.objects = []

class HaloCatalogHaloContainer(YTSelectionContainer):
    """
    Create a data container to get member particles and individual
    values from halos and subhalos. Halo mass, position, and
    velocity are set as attributes. Halo IDs are accessible
    through the field, "member_ids".  Other fields that are one
    value per halo are accessible as normal.  The field list for
    halo objects can be seen in `ds.halos_field_list`.

    Parameters
    ----------
    ptype : string
        The type of halo, either "Group" for the main halo or
        "Subhalo" for subhalos.
    particle_identifier : int or tuple of ints
        The halo or subhalo id.  If requesting a subhalo, the id
        can also be given as a tuple of the main halo id and
        subgroup id, such as (1, 4) for subgroup 4 of halo 1.

    Attributes
    ----------
    particle_identifier : int
        The id of the halo or subhalo.
    group_identifier : int
        For subhalos, the id of the enclosing halo.
    subgroup_identifier : int
        For subhalos, the relative id of the subhalo within
        the enclosing halo.
    particle_number : int
        Number of particles in the halo.
    mass : float
        Halo mass.
    position : array of floats
        Halo position.
    velocity : array of floats
        Halo velocity.

    Note
    ----
    Relevant Fields:

     * particle_number - number of particles
     * subhalo_number - number of subhalos
     * group_identifier - id of parent group for subhalos

    Examples
    --------

    >>> import yt
    >>> ds = yt.load("gadget_halos/data/groups_298/fof_subhalo_tab_298.0.hdf5")
    >>>
    >>> halo = ds.halo("Group", 0)
    >>> print(halo.mass)
    13256.5517578 code_mass
    >>> print(halo.position)
    [ 16.18603706   6.95965052  12.52694607] code_length
    >>> print(halo.velocity)
    [ 6943694.22793569  -762788.90647454  -794749.63819757] cm/s
    >>> print(halo["Group_R_Crit200"])
    [ 0.79668683] code_length
    >>>
    >>> # particle ids for this halo
    >>> print(halo["member_ids"])
    [  723631.   690744.   854212. ...,   608589.   905551.  1147449.] dimensionless
    >>>
    >>> # get the first subhalo of this halo
    >>> subhalo = ds.halo("Subhalo", (0, 0))
    >>> print(subhalo["member_ids"])
    [  723631.   690744.   854212. ...,   808362.   956359.  1248821.] dimensionless

    """

    _type_name = "halo"
    _con_args = ("ptype", "particle_identifier")
    _spatial = False

    def __init__(self, ptype, particle_identifier, ds=None):
        if ptype not in ds.particle_types_raw:
            raise RuntimeError("Possible halo types are %s, supplied \"%s\"." %
                               (ds.particle_types_raw, ptype))

        self.ptype = ptype
        self._current_particle_type = ptype
        super(HaloCatalogHaloContainer, self).__init__(ds, {})

        if ptype == "Subhalo" and isinstance(particle_identifier, tuple):
            self.group_identifier, self.subgroup_identifier = \
              particle_identifier
            my_data = self.index._get_halo_values(
                "Group", np.array([self.group_identifier]),
                ["GroupFirstSub"])
            self.particle_identifier = \
              np.int64(my_data["GroupFirstSub"][0] + self.subgroup_identifier)
        else:
            self.particle_identifier = particle_identifier

        if self.particle_identifier >= self.index.particle_count[ptype]:
            raise RuntimeError("%s %d requested, but only %d %s objects exist." %
                               (ptype, particle_identifier,
                                self.index.particle_count[ptype], ptype))

        # Find the file that has the scalar values for this halo.
        i_scalar = self.index._get_halo_file_indices(
            ptype, [self.particle_identifier])[0]
        self.scalar_data_file = self.index.data_files[i_scalar]

        # index within halo arrays that corresponds to this halo
        self.scalar_index = self.index._get_halo_scalar_index(
            ptype, self.particle_identifier)

        halo_fields = ["%sLen" % ptype]
        if ptype == "Subhalo": halo_fields.append("SubhaloGrNr")
        my_data = self.index._get_halo_values(
            ptype, np.array([self.particle_identifier]),
            halo_fields)
        self.particle_number = np.int64(my_data["%sLen" % ptype][0])

        if ptype == "Group":
            self.group_identifier = self.particle_identifier
            id_offset = 0
            # index of file that has scalar values for the group
            g_scalar = i_scalar
            group_index = self.scalar_index

        # If a subhalo, find the index of the parent.
        elif ptype == "Subhalo":
            self.group_identifier = np.int64(my_data["SubhaloGrNr"][0])

            # Find the file that has the scalar values for the parent group.
            g_scalar = self.index._get_halo_file_indices(
                "Group", [self.group_identifier])[0]

            # index within halo arrays that corresponds to the paent group
            group_index = self.index._get_halo_scalar_index(
                "Group", self.group_identifier)

            my_data = self.index._get_halo_values(
                "Group", np.array([self.group_identifier]),
                ["GroupNsubs", "GroupFirstSub"])
            self.subgroup_identifier = self.particle_identifier - \
              np.int64(my_data["GroupFirstSub"][0])
            parent_subhalos = my_data["GroupNsubs"][0]

            mylog.debug("Subhalo %d is subgroup %s of %d in group %d." % \
                        (self.particle_identifier, self.subgroup_identifier,
                        parent_subhalos, self.group_identifier))

            # ids of the sibling subhalos that come before this one
            if self.subgroup_identifier > 0:
                sub_ids = np.arange(
                    self.particle_identifier - self.subgroup_identifier,
                    self.particle_identifier)
                my_data = self.index._get_halo_values(
                    "Subhalo", sub_ids, ["SubhaloLen"])
                id_offset = my_data["SubhaloLen"].sum(dtype=np.int64)
            else:
                id_offset = 0

        # Calculate the starting index for the member particles.
        # First, add up all the particles in the earlier files.
        all_id_start = self.index._group_length_sum[:g_scalar].sum(dtype=np.int64)

        # Now add the halos in this file that come before.
        with h5py.File(self.index.data_files[g_scalar].filename, "r") as f:
            all_id_start += f["Group"]["GroupLen"][:group_index].sum(dtype=np.int64)

        # Add the subhalo offset.
        all_id_start += id_offset

        # indices of first and last files containing member particles
        i_start = np.digitize([all_id_start],
                              self.index._halo_id_start,
                              right=False)[0] - 1
        i_end = np.digitize([all_id_start+self.particle_number],
                            self.index._halo_id_end,
                            right=True)[0]
        self.field_data_files = self.index.data_files[i_start:i_end+1]

        # starting and ending indices for each file containing particles
        self.field_data_start = \
          (all_id_start -
           self.index._halo_id_start[i_start:i_end+1]).clip(min=0)
        self.field_data_start = self.field_data_start.astype(np.int64)
        self.field_data_end = \
          (all_id_start + self.particle_number -
           self.index._halo_id_start[i_start:i_end+1]).clip(
               max=self.index._halo_id_number[i_start:i_end+1])
        self.field_data_end = self.field_data_end.astype(np.int64)

        for attr in ["mass", "position", "velocity"]:
            setattr(self, attr, self[self.ptype, "particle_%s" % attr][0])

    def __repr__(self):
        return "%s_%s_%09d" % \
          (self.ds, self.ptype, self.particle_identifier)
