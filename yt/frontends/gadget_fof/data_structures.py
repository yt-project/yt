from collections import defaultdict
from functools import partial
from yt.utilities.on_demand_imports import _h5py as h5py
import numpy as np
import stat
import os
import weakref

from yt.data_objects.data_containers import \
    YTSelectionContainer
from yt.data_objects.static_output import \
    ParticleDataset
from yt.frontends.gadget.data_structures import \
    _fix_unit_ordering
from yt.frontends.gadget_fof.fields import \
    GadgetFOFFieldInfo, \
    GadgetFOFHaloFieldInfo
from yt.frontends.halo_catalog.data_structures import \
    HaloCatalogFile, \
    HaloCatalogParticleIndex, \
    HaloContainer, \
    HaloDatasetParticleIndex, \
    HaloDataset
from yt.funcs import \
    only_on_root, \
    setdefaultattr
from yt.geometry.particle_geometry_handler import \
    ParticleIndex
from yt.utilities.cosmology import \
    Cosmology
from yt.utilities.logger import ytLogger as \
    mylog

class GadgetFOFParticleIndex(HaloCatalogParticleIndex):
    def _calculate_file_offset_map(self):
        # After the FOF is performed, a load-balancing step redistributes halos
        # and then writes more fields.  Here, for each file, we create a list of
        # files which contain the rest of the redistributed particles.
        ifof = np.array([data_file.total_particles["Group"]
                         for data_file in self.data_files])
        isub = np.array([data_file.total_offset
                         for data_file in self.data_files])
        subend = isub.cumsum()
        fofend = ifof.cumsum()
        istart = np.digitize(fofend - ifof, subend - isub) - 1
        iend = np.clip(np.digitize(fofend, subend), 0, ifof.size - 2)
        for i, data_file in enumerate(self.data_files):
            data_file.offset_files = self.data_files[istart[i]: iend[i] + 1]

    def _detect_output_fields(self):
        field_list = []
        units = {}
        found_fields = \
          dict([(ptype, False)
                for ptype, pnum in self.particle_count.items()
                if pnum > 0])

        for data_file in self.data_files:
            fl, _units = self.io._identify_fields(data_file)
            units.update(_units)
            field_list.extend([f for f in fl if f not in field_list])
            for ptype in found_fields:
                found_fields[ptype] |= data_file.total_particles[ptype]
            if all(found_fields.values()): break

        self.field_list = field_list
        ds = self.dataset
        ds.particle_types = tuple(set(pt for pt, ds in field_list))
        ds.field_units.update(units)
        ds.particle_types_raw = ds.particle_types

    def _setup_filenames(self):
        if not hasattr(self, "data_files"):
            template = self.ds.filename_template
            ndoms = self.ds.file_count
            cls = self.ds._file_class
            self.data_files = [
                cls(self.ds, self.io, template % {'num':i}, i, frange=None)
                for i in range(ndoms)
            ]
        if not hasattr(self, "total_particles"):
            self.total_particles = sum(
                sum(d.total_particles.values()) for d in self.data_files
            )

    def _setup_data_io(self):
        super(GadgetFOFParticleIndex, self)._setup_data_io()
        self._setup_filenames()
        self._calculate_particle_count()
        self._calculate_particle_index_starts()
        self._calculate_file_offset_map()

class GadgetFOFHDF5File(HaloCatalogFile):
    def __init__(self, ds, io, filename, file_id, frange):
        with h5py.File(filename, mode="r") as f:
            self.header = \
              dict((str(field), val)
                   for field, val in f["Header"].attrs.items())
            self.group_length_sum = f["Group/GroupLen"][()].sum() \
              if "Group/GroupLen" in f else 0
            self.group_subs_sum = f["Group/GroupNsubs"][()].sum() \
              if "Group/GroupNsubs" in f else 0
        self.total_ids = self.header["Nids_ThisFile"]
        self.total_offset = 0
        super(GadgetFOFHDF5File, self).__init__(
            ds, io, filename, file_id, frange)

    def _read_particle_positions(self, ptype, f=None):
        """
        Read all particle positions in this file.
        """

        if f is None:
            close = True
            f = h5py.File(self.filename, mode="r")
        else:
            close = False

        pos = f[ptype]["%sPos" % ptype][()].astype("float64")

        if close:
            f.close()

        return pos

class GadgetFOFDataset(ParticleDataset):
    _index_class = GadgetFOFParticleIndex
    _file_class = GadgetFOFHDF5File
    _field_info_class = GadgetFOFFieldInfo

    def __init__(self, filename, dataset_type="gadget_fof_hdf5",
                 index_order=None, index_filename=None,
                 unit_base=None, units_override=None, unit_system="cgs"):
        if unit_base is not None and "UnitLength_in_cm" in unit_base:
            # We assume this is comoving, because in the absence of comoving
            # integration the redshift will be zero.
            unit_base['cmcm'] = 1.0 / unit_base["UnitLength_in_cm"]
        self._unit_base = unit_base
        if units_override is not None:
            raise RuntimeError("units_override is not supported for GadgetFOFDataset. "+
                               "Use unit_base instead.")
        super(GadgetFOFDataset, self).__init__(
            filename, dataset_type, units_override=units_override,
            index_order=index_order, index_filename=index_filename,
            unit_system=unit_system)

    def add_field(self, *args, **kwargs):
        super(GadgetFOFDataset, self).add_field(*args, **kwargs)
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
            self._instantiated_halo_ds = GadgetFOFHaloDataset(self)
        return self._instantiated_halo_ds

    def _setup_classes(self):
        super(GadgetFOFDataset, self)._setup_classes()
        self.halo = partial(GadgetFOFHaloContainer, ds=self._halos_ds)

    def _parse_parameter_file(self):
        with h5py.File(self.parameter_filename, mode="r") as f:
            self.parameters = \
              dict((str(field), val)
                   for field, val in f["Header"].attrs.items())

        self.dimensionality = 3
        self.refine_by = 2
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])

        # Set standard values
        self.domain_left_edge = np.zeros(3, "float64")
        self.domain_right_edge = np.ones(3, "float64") * \
          self.parameters["BoxSize"]
        self.domain_dimensions = np.ones(3, "int32")
        self.cosmological_simulation = 1
        self.periodicity = (True, True, True)
        self.current_redshift = self.parameters["Redshift"]
        self.omega_lambda = self.parameters["OmegaLambda"]
        self.omega_matter = self.parameters["Omega0"]
        self.hubble_constant = self.parameters["HubbleParam"]
        cosmology = Cosmology(hubble_constant=self.hubble_constant,
                              omega_matter=self.omega_matter,
                              omega_lambda=self.omega_lambda)
        self.current_time = cosmology.t_from_z(self.current_redshift)

        prefix = os.path.abspath(
            os.path.join(os.path.dirname(self.parameter_filename), 
                         os.path.basename(self.parameter_filename).split(".", 1)[0]))
        suffix = self.parameter_filename.rsplit(".", 1)[-1]
        self.filename_template = "%s.%%(num)i.%s" % (prefix, suffix)
        self.file_count = self.parameters["NumFiles"]
        self.particle_types = ("Group", "Subhalo")
        self.particle_types_raw = ("Group", "Subhalo")

    def _set_code_unit_attributes(self):
        # Set a sane default for cosmological simulations.
        if self._unit_base is None and self.cosmological_simulation == 1:
            only_on_root(mylog.info, "Assuming length units are in Mpc/h (comoving)")
            self._unit_base = dict(length = (1.0, "Mpccm/h"))
        # The other same defaults we will use from the standard Gadget
        # defaults.
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
        setdefaultattr(self, 'length_unit',
                       self.quan(length_unit[0], length_unit[1]))
        
        if "velocity" in unit_base:
            velocity_unit = unit_base["velocity"]
        elif "UnitVelocity_in_cm_per_s" in unit_base:
            velocity_unit = (unit_base["UnitVelocity_in_cm_per_s"], "cm/s")
        else:
            if self.cosmological_simulation == 0:
                velocity_unit = (1e5, "cm/s")
            else:
                velocity_unit = (1e5, "cm/s * sqrt(a)")
        velocity_unit = _fix_unit_ordering(velocity_unit)
        setdefaultattr(self, 'velocity_unit',
                       self.quan(velocity_unit[0], velocity_unit[1]))

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
            mass_unit = (1.0, "1e10*Msun/h")
        mass_unit = _fix_unit_ordering(mass_unit)
        setdefaultattr(self, 'mass_unit', self.quan(mass_unit[0], mass_unit[1]))

        if "time" in unit_base:
            time_unit = unit_base["time"]
        elif "UnitTime_in_s" in unit_base:
            time_unit = (unit_base["UnitTime_in_s"], "s")
        else:
            tu = (self.length_unit / self.velocity_unit).to("yr/h")
            time_unit = (tu.d, tu.units)
        setdefaultattr(self, 'time_unit', self.quan(time_unit[0], time_unit[1]))

    def __repr__(self):
        return self.basename.split(".", 1)[0]

    @classmethod
    def _is_valid(self, *args, **kwargs):
        need_groups = ['Group', 'Header', 'Subhalo']
        veto_groups = ['FOF']
        valid = True
        try:
            fh = h5py.File(args[0], mode='r')
            valid = all(ng in fh["/"] for ng in need_groups) and \
              not any(vg in fh["/"] for vg in veto_groups)
            fh.close()
        except:
            valid = False
            pass
        return valid

class GadgetFOFHaloParticleIndex(GadgetFOFParticleIndex, HaloDatasetParticleIndex):
    _detect_output_fields = HaloDatasetParticleIndex._detect_output_fields
    _setup_data_io = GadgetFOFParticleIndex._setup_data_io

    def _setup_data_io(self):
        super(GadgetFOFHaloParticleIndex, self)._setup_data_io()
        self._create_halo_id_table()

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

    def _read_halo_particle_field(self, fh, ptype, field, indices):
        return fh[os.path.join(ptype, field)][indices]

class GadgetFOFHaloDataset(HaloDataset):
    _index_class = GadgetFOFHaloParticleIndex
    _file_class = GadgetFOFHDF5File
    _field_info_class = GadgetFOFHaloFieldInfo

    def __init__(self, ds, dataset_type="gadget_fof_halo_hdf5"):
        super(GadgetFOFHaloDataset, self).__init__(ds, dataset_type)

class GadgetFOFHaloContainer(HaloContainer):
    def _get_member_fieldnames(self):
        halo_fields = ["%sLen" % self.ptype]
        if self.ptype == "Subhalo":
            halo_fields.append("SubhaloGrNr")
        return halo_fields

    def _get_particle_number(self):
        return self._io_data["%sLen" % self.ptype]

    def _set_field_indices(self):
        if self.ptype == "Group":
            self.group_identifier = self.particle_identifier
            id_offset = 0
            # index of file that has scalar values for the group
            g_scalar = self.i_scalar
            group_index = self.scalar_index

        # If a subhalo, find the index of the parent.
        elif self.ptype == "Subhalo":
            self.group_identifier = self._io_data["SubhaloGrNr"]

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
        with h5py.File(self.index.data_files[g_scalar].filename, mode="r") as f:
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

    def _set_identifiers(self, particle_identifier):
        if self.ptype == "Subhalo" and isinstance(particle_identifier, tuple):
            self.group_identifier, self.subgroup_identifier = \
              particle_identifier
            my_data = self.index._get_halo_values(
                "Group", np.array([self.group_identifier]),
                ["GroupFirstSub"])
            self.particle_identifier = \
              np.int64(my_data["GroupFirstSub"][0] + self.subgroup_identifier)
        else:
            self.particle_identifier = particle_identifier
