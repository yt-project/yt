"""
Data structures for Gadget frontend




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import h5py
import numpy as np
import stat
import os
import types

from yt.data_objects.static_output import \
    ParticleFile
from yt.frontends.sph.data_structures import \
    ParticleDataset
from yt.frontends.sph.fields import \
    SPHFieldInfo
from yt.geometry.particle_geometry_handler import \
    ParticleIndex
from yt.utilities.cosmology import \
    Cosmology
from yt.utilities.definitions import \
    sec_conversion
from yt.utilities.fortran_utils import read_record
from yt.utilities.logger import ytLogger as mylog

from .definitions import \
    gadget_header_specs, \
    gadget_field_specs, \
    gadget_ptype_specs

def _fix_unit_ordering(unit):
    if isinstance(unit[0], types.StringTypes):
        unit = unit[1], unit[0]
    return unit

class GadgetBinaryFile(ParticleFile):
    def __init__(self, ds, io, filename, file_id):
        with open(filename, "rb") as f:
            self.header = read_record(f, ds._header_spec)
            self._position_offset = f.tell()
            f.seek(0, os.SEEK_END)
            self._file_size = f.tell()

        super(GadgetBinaryFile, self).__init__(ds, io, filename, file_id)

    def _calculate_offsets(self, field_list):
        self.field_offsets = self.io._calculate_field_offsets(
            field_list, self.total_particles,
            self._position_offset, self._file_size)

class GadgetDataset(ParticleDataset):
    _index_class = ParticleIndex
    _file_class = GadgetBinaryFile
    _field_info_class = SPHFieldInfo
    _particle_mass_name = "Mass"
    _particle_coordinates_name = "Coordinates"
    _particle_velocity_name = "Velocities"
    _suffix = ""

    def __init__(self, filename, dataset_type="gadget_binary",
                 additional_fields=(),
                 unit_base=None, n_ref=64,
                 over_refine_factor=1,
                 bounding_box = None,
                 header_spec = "default",
                 field_spec = "default",
                 ptype_spec = "default",
                 units_override=None):
        if self._instantiated: return
        self._header_spec = self._setup_binary_spec(
            header_spec, gadget_header_specs)
        self._field_spec = self._setup_binary_spec(
            field_spec, gadget_field_specs)
        self._ptype_spec = self._setup_binary_spec(
            ptype_spec, gadget_ptype_specs)
        self.n_ref = n_ref
        self.over_refine_factor = over_refine_factor
        self.storage_filename = None
        if unit_base is not None and "UnitLength_in_cm" in unit_base:
            # We assume this is comoving, because in the absence of comoving
            # integration the redshift will be zero.
            unit_base['cmcm'] = 1.0 / unit_base["UnitLength_in_cm"]
        self._unit_base = unit_base
        if bounding_box is not None:
            bbox = np.array(bounding_box, dtype="float64")
            if bbox.shape == (2, 3):
                bbox = bbox.transpose()
            self.domain_left_edge = bbox[:,0]
            self.domain_right_edge = bbox[:,1]
        else:
            self.domain_left_edge = self.domain_right_edge = None
        if units_override is not None:
            raise RuntimeError("units_override is not supported for GadgetDataset. "+
                               "Use unit_base instead.")
        super(GadgetDataset, self).__init__(filename, dataset_type)

    def _setup_binary_spec(self, spec, spec_dict):
        if isinstance(spec, types.StringTypes):
            _hs = ()
            for hs in spec.split("+"):
                _hs += spec_dict[hs]
            spec = _hs
        return spec

    def __repr__(self):
        return os.path.basename(self.parameter_filename).split(".")[0]

    def _get_hvals(self):
        # The entries in this header are capitalized and named to match Table 4
        # in the GADGET-2 user guide.

        f = open(self.parameter_filename)
        hvals = read_record(f, self._header_spec)
        for i in hvals:
            if len(hvals[i]) == 1:
                hvals[i] = hvals[i][0]
        return hvals

    def _parse_parameter_file(self):

        hvals = self._get_hvals()

        self.dimensionality = 3
        self.refine_by = 2
        self.parameters["HydroMethod"] = "sph"
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])
        # Set standard values

        # We may have an overridden bounding box.
        if self.domain_left_edge is None:
            self.domain_left_edge = np.zeros(3, "float64")
            self.domain_right_edge = np.ones(3, "float64") * hvals["BoxSize"]
        nz = 1 << self.over_refine_factor
        self.domain_dimensions = np.ones(3, "int32") * nz
        self.periodicity = (True, True, True)

        self.cosmological_simulation = 1

        self.current_redshift = hvals["Redshift"]
        self.omega_lambda = hvals["OmegaLambda"]
        self.omega_matter = hvals["Omega0"]
        self.hubble_constant = hvals["HubbleParam"]
        # According to the Gadget manual, OmegaLambda will be zero for
        # non-cosmological datasets.  However, it may be the case that
        # individuals are running cosmological simulations *without* Lambda, in
        # which case we may be doing something incorrect here.
        # It may be possible to deduce whether ComovingIntegration is on
        # somehow, but opinions on this vary.
        if self.omega_lambda == 0.0:
            mylog.info("Omega Lambda is 0.0, so we are turning off Cosmology.")
            self.hubble_constant = 1.0  # So that scaling comes out correct
            self.cosmological_simulation = 0
            self.current_redshift = 0.0
            # This may not be correct.
            self.current_time = hvals["Time"] * sec_conversion["Gyr"]
        else:
            # Now we calculate our time based on the cosmology, because in
            # ComovingIntegration hvals["Time"] will in fact be the expansion
            # factor, not the actual integration time, so we re-calculate
            # global time from our Cosmology.
            cosmo = Cosmology(self.hubble_constant,
                              self.omega_matter, self.omega_lambda)
            self.current_time = cosmo.hubble_time(self.current_redshift)
            mylog.info("Calculating time from %0.3e to be %0.3e seconds",
                       hvals["Time"], self.current_time)
        self.parameters = hvals

        prefix = os.path.abspath(
            os.path.join(os.path.dirname(self.parameter_filename),
                         os.path.basename(self.parameter_filename).split(".", 1)[0]))

        if hvals["NumFiles"] > 1:
            self.filename_template = "%s.%%(num)s%s" % (prefix, self._suffix)
        else:
            self.filename_template = self.parameter_filename

        self.file_count = hvals["NumFiles"]

    def _set_code_unit_attributes(self):
        # If no units passed in by user, set a sane default (Gadget-2 users guide).
        if self._unit_base is None:
            if self.cosmological_simulation == 1:
                mylog.info("Assuming length units are in kpc/h (comoving)")
                self._unit_base = dict(length = (1.0, "kpccm/h"))
            else:
                mylog.info("Assuming length units are in kpc (physical)")
                self._unit_base = dict(length = (1.0, "kpc"))
                
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
        if "velocity" in unit_base:
            velocity_unit = unit_base["velocity"]
        elif "UnitVelocity_in_cm_per_s" in unit_base:
            velocity_unit = (unit_base["UnitVelocity_in_cm_per_s"], "cm/s")
        else:
            velocity_unit = (1e5, "cm/s")
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
            mass_unit = (1.0, "1e10*Msun/h")
        mass_unit = _fix_unit_ordering(mass_unit)
        self.mass_unit = self.quan(mass_unit[0], mass_unit[1])
        self.time_unit = self.length_unit / self.velocity_unit

    @classmethod
    def _is_valid(self, *args, **kwargs):
        # We do not allow load() of these files.
        return False

class GadgetHDF5Dataset(GadgetDataset):
    _file_class = ParticleFile
    _field_info_class = SPHFieldInfo
    _particle_mass_name = "Masses"
    _suffix = ".hdf5"

    def __init__(self, filename, dataset_type="gadget_hdf5", 
                 unit_base = None, n_ref=64,
                 over_refine_factor=1,
                 bounding_box = None,
                 units_override=None):
        self.storage_filename = None
        filename = os.path.abspath(filename)
        if units_override is not None:
            raise RuntimeError("units_override is not supported for GadgetHDF5Dataset. "+
                               "Use unit_base instead.")
        super(GadgetHDF5Dataset, self).__init__(
            filename, dataset_type, unit_base=unit_base, n_ref=n_ref,
            over_refine_factor=over_refine_factor,
            bounding_box = bounding_box)

    def _get_hvals(self):
        handle = h5py.File(self.parameter_filename, mode="r")
        hvals = {}
        hvals.update((str(k), v) for k, v in handle["/Header"].attrs.items())
        # Compat reasons.
        hvals["NumFiles"] = hvals["NumFilesPerSnapshot"]
        hvals["Massarr"] = hvals["MassTable"]
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
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])

        self._unit_base = self._get_uvals()
        self._unit_base['cmcm'] = 1.0 / self._unit_base["UnitLength_in_cm"]

        self.current_redshift = self.parameters["Redshift"]
        self.omega_lambda = self.parameters["OmegaLambda"]
        self.omega_matter = self.parameters["Omega0"]
        self.hubble_constant = self.parameters["HubbleParam"]

        if self.domain_left_edge is None:
            self.domain_left_edge = np.zeros(3, "float64")
            self.domain_right_edge = np.ones(3, "float64") * self.parameters["BoxSize"]

        nz = 1 << self.over_refine_factor
        self.domain_dimensions = np.ones(3, "int32") * nz

        self.cosmological_simulation = 1
        self.periodicity = (True, True, True)

        prefix = os.path.abspath(
            os.path.join(os.path.dirname(self.parameter_filename),
                         os.path.basename(self.parameter_filename).split(".", 1)[0]))

        suffix = self.parameter_filename.rsplit(".", 1)[-1]
        self.filename_template = "%s.%%(num)i.%s" % (prefix, suffix)
        self.file_count = self.parameters["NumFilesPerSnapshot"]

    def _set_owls_eagle_units(self):

        # note the contents of the HDF5 Units group are in _unit_base 
        # note the velocity stored on disk is sqrt(a) dx/dt 
        self.length_unit = self.quan(self._unit_base["UnitLength_in_cm"], 'cmcm/h')
        self.mass_unit = self.quan(self._unit_base["UnitMass_in_g"], 'g/h')
        self.velocity_unit = self.quan(self._unit_base["UnitVelocity_in_cm_per_s"], 'cm/s')
        self.time_unit = self.quan(self._unit_base["UnitTime_in_s"], 's/h')

    @classmethod
    def _is_valid(self, *args, **kwargs):
        try:
            fileh = h5py.File(args[0], mode='r')
            if "Constants" not in fileh["/"].keys() and \
               "Header" in fileh["/"].keys():
                fileh.close()
                return True
            fileh.close()
        except:
            pass
        return False
