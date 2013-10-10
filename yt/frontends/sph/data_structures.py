"""
Data structures for a generic SPH/Gadget frontend.




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import h5py
import numpy as np
import stat
import weakref
import struct
import glob
import os

from yt.utilities.fortran_utils import read_record
from yt.utilities.logger import ytLogger as mylog
from yt.geometry.particle_geometry_handler import \
    ParticleGeometryHandler
from yt.data_objects.static_output import \
    StaticOutput, ParticleFile
from yt.utilities.definitions import \
    mpc_conversion, sec_conversion
from yt.utilities.physical_constants import \
    G, \
    cm_per_kpc, \
    mass_sun_cgs
from yt.utilities.cosmology import Cosmology
from .fields import \
    OWLSFieldInfo, \
    KnownOWLSFields, \
    GadgetFieldInfo, \
    KnownGadgetFields, \
    GadgetHDF5FieldInfo, \
    KnownGadgetHDF5Fields, \
    TipsyFieldInfo, \
    KnownTipsyFields


class GadgetBinaryFile(ParticleFile):
    def __init__(self, pf, io, filename, file_id):
        with open(filename, "rb") as f:
            self.header = read_record(f, pf._header_spec)
            self._position_offset = f.tell()
            f.seek(0, os.SEEK_END)
            self._file_size = f.tell()

        super(GadgetBinaryFile, self).__init__(pf, io, filename, file_id)

    def _calculate_offsets(self, field_list):
        self.field_offsets = self.io._calculate_field_offsets(
            field_list, self.total_particles,
            self._position_offset, self._file_size)


class ParticleStaticOutput(StaticOutput):
    _unit_base = None
    over_refine_factor = 1

    def _set_units(self):
        self.units = {}
        self.time_units = {}
        self.conversion_factors = {}
        self.units['1'] = 1.0
        DW = self.domain_right_edge - self.domain_left_edge
        self.units["unitary"] = 1.0 / DW.max()
        # Check
        base = None
        mpch = {}
        mpch.update(mpc_conversion)
        unit_base = self._unit_base or {}
        if self.cosmological_simulation:
            for unit in mpc_conversion:
                mpch['%sh' % unit] = mpch[unit] * self.hubble_constant
                mpch['%shcm' % unit] = (mpch["%sh" % unit] *
                                        (1 + self.current_redshift))
                mpch['%scm' % unit] = mpch[unit] * (1 + self.current_redshift)
        elif 'cmcm' in unit_base:
            unit_base['cm'] = self.units['cm'] = unit_base.pop("cmcm")
        # ud == unit destination
        # ur == unit registry
        for ud, ur in [(self.units, mpch), (self.time_units, sec_conversion)]:
            for unit in sorted(unit_base):
                if unit in ur:
                    ratio = (ur[unit] / ur['mpc'])
                    base = unit_base[unit] * ratio
                    break
            if base is None:
                continue
            for unit in ur:
                ud[unit] = ur[unit] / base


class GadgetStaticOutput(ParticleStaticOutput):
    _hierarchy_class = ParticleGeometryHandler
    _file_class = GadgetBinaryFile
    _fieldinfo_fallback = GadgetFieldInfo
    _fieldinfo_known = KnownGadgetFields
    _particle_mass_name = "Mass"
    _particle_coordinates_name = "Coordinates"
    _suffix = ""
    _header_spec = (('Npart', 6, 'i'),
                    ('Massarr', 6, 'd'),
                    ('Time', 1, 'd'),
                    ('Redshift', 1, 'd'),
                    ('FlagSfr', 1, 'i'),
                    ('FlagFeedback', 1, 'i'),
                    ('Nall', 6, 'i'),
                    ('FlagCooling', 1, 'i'),
                    ('NumFiles', 1, 'i'),
                    ('BoxSize', 1, 'd'),
                    ('Omega0', 1, 'd'),
                    ('OmegaLambda', 1, 'd'),
                    ('HubbleParam', 1, 'd'),
                    ('FlagAge', 1, 'i'),
                    ('FlagMEtals', 1, 'i'),
                    ('NallHW', 6, 'i'),
                    ('unused', 16, 'i'))

    def __init__(self, filename, data_style="gadget_binary",
                 additional_fields=(),
                 unit_base=None, n_ref=64,
                 over_refine_factor=1,
                 bounding_box = None):
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
            self.domain_left_edge = bbox[:,0]
            self.domain_right_edge = bbox[:,1]
        else:
            self.domain_left_edge = self.domain_right_edge = None
        super(GadgetStaticOutput, self).__init__(filename, data_style)

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
            cosmo = Cosmology(self.hubble_constant * 100.0,
                              self.omega_matter, self.omega_lambda)
            self.current_time = cosmo.UniverseAge(self.current_redshift)
            mylog.info("Calculating time from %0.3e to be %0.3e seconds",
                       hvals["Time"], self.current_time)
        self.parameters = hvals

        prefix = self.parameter_filename.split(".", 1)[0]

        if hvals["NumFiles"] > 1:
            self.filename_template = "%s.%%(num)s%s" % (prefix, self._suffix)
        else:
            self.filename_template = self.parameter_filename

        self.file_count = hvals["NumFiles"]

    def _set_units(self):
        if self._unit_base is None and self.cosmological_simulation == 1:
            mylog.info("Assuming length units are in Mpc/h (comoving)")
            self._unit_base = dict(mpchcm = 1.0)
        super(GadgetStaticOutput, self)._set_units()
        length_unit = self.units['cm']
        unit_base = self._unit_base or {}
        velocity_unit = unit_base.get("velocity", 1e5)
        velocity_unit = unit_base.get("UnitVelocity_in_cm_per_s",
                                      velocity_unit)
        # We set hubble_constant = 1.0 for non-cosmology
        msun10 = mass_sun_cgs * 1e10 / self.hubble_constant
        mass_unit = unit_base.get("g", msun10)
        mass_unit = unit_base.get("UnitMass_in_g", mass_unit)
        self.conversion_factors["velocity"] = velocity_unit
        self.conversion_factors["mass"] = mass_unit
        self.conversion_factors["density"] = mass_unit / length_unit ** 3
        # Currently, setting time_units is disabled.  The current_time is
        # accurately set, but until a time that we can confirm how
        # FormationTime for stars is set I am disabling these.
        #time_unit = length_unit / velocity_unit
        #for u in sec_conversion:
        #    self.time_units[u] = time_unit / sec_conversion[u]

    @classmethod
    def _is_valid(self, *args, **kwargs):
        # We do not allow load() of these files.
        return False


class OWLSStaticOutput(GadgetStaticOutput):
    _hierarchy_class = ParticleGeometryHandler
    _file_class = ParticleFile
    _fieldinfo_fallback = OWLSFieldInfo  # For now we have separate from Gadget
    _fieldinfo_known = KnownOWLSFields
    _particle_mass_name = "Mass"
    _particle_coordinates_name = "Coordinates"
    _header_spec = None  # Override so that there's no confusion

    def __init__(self, filename, data_style="OWLS", n_ref=64,
                 over_refine_factor=1):
        self.storage_filename = None
        filename = os.path.abspath(filename)
        super(OWLSStaticOutput, self).__init__(
            filename, data_style, unit_base=None, n_ref=n_ref,
            over_refine_factor=over_refine_factor)

    def __repr__(self):
        return os.path.basename(self.parameter_filename).split(".")[0]

    def _parse_parameter_file(self):
        handle = h5py.File(self.parameter_filename, mode="r")
        hvals = {}
        hvals.update((str(k), v) for k, v in handle["/Header"].attrs.items())

        self.dimensionality = 3
        self.refine_by = 2
        self.parameters["HydroMethod"] = "sph"
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])

        # Set standard values
        self.current_time = hvals["Time_GYR"] * sec_conversion["Gyr"]
        self.domain_left_edge = np.zeros(3, "float64")
        self.domain_right_edge = np.ones(3, "float64") * hvals["BoxSize"]
        nz = 1 << self.over_refine_factor
        self.domain_dimensions = np.ones(3, "int32") * nz
        self.cosmological_simulation = 1
        self.periodicity = (True, True, True)
        self.current_redshift = hvals["Redshift"]
        self.omega_lambda = hvals["OmegaLambda"]
        self.omega_matter = hvals["Omega0"]
        self.hubble_constant = hvals["HubbleParam"]
        self.parameters = hvals

        prefix = os.path.abspath(self.parameter_filename.split(".", 1)[0])
        suffix = self.parameter_filename.rsplit(".", 1)[-1]
        self.filename_template = "%s.%%(num)i.%s" % (prefix, suffix)
        self.file_count = hvals["NumFilesPerSnapshot"]

        # To avoid having to open files twice
        self._unit_base = {}
        self._unit_base.update(
            (str(k), v) for k, v in handle["/Units"].attrs.items())
        # Comoving cm is given in the Units
        self._unit_base['cmcm'] = 1.0 / self._unit_base["UnitLength_in_cm"]

        handle.close()

    @classmethod
    def _is_valid(self, *args, **kwargs):
        try:
            fileh = h5py.File(args[0], mode='r')
            if "Constants" in fileh["/"].keys() and \
               "Header" in fileh["/"].keys():
                fileh.close()
                return True
            fileh.close()
        except:
            pass
        return False

class GadgetHDF5StaticOutput(GadgetStaticOutput):
    _file_class = ParticleFile
    _fieldinfo_fallback = GadgetHDF5FieldInfo
    _fieldinfo_known = KnownGadgetHDF5Fields
    _suffix = ".hdf5"

    def __init__(self, filename, data_style="gadget_hdf5", 
                 unit_base = None, n_ref=64,
                 over_refine_factor=1,
                 bounding_box = None):
        self.storage_filename = None
        filename = os.path.abspath(filename)
        super(GadgetHDF5StaticOutput, self).__init__(
            filename, data_style, unit_base=unit_base, n_ref=n_ref,
            over_refine_factor=over_refine_factor,
            bounding_box = bounding_box)

    def _get_hvals(self):
        handle = h5py.File(self.parameter_filename, mode="r")
        hvals = {}
        hvals.update((str(k), v) for k, v in handle["/Header"].attrs.items())
        # Compat reasons.
        hvals["NumFiles"] = hvals["NumFilesPerSnapshot"]
        hvals["Massarr"] = hvals["MassTable"]
        return hvals

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

class TipsyFile(ParticleFile):

    def _calculate_offsets(self, field_list):
        self.field_offsets = self.io._calculate_particle_offsets(self)

    def __init__(self, pf, io, filename, file_id):
        # To go above 1 domain, we need to include an indexing step in the
        # IOHandler, rather than simply reading from a single file.
        assert file_id == 0
        super(TipsyFile, self).__init__(pf, io, filename, file_id)
        io._create_dtypes(self)


class TipsyStaticOutput(ParticleStaticOutput):
    _hierarchy_class = ParticleGeometryHandler
    _file_class = TipsyFile
    _fieldinfo_fallback = TipsyFieldInfo
    _fieldinfo_known = KnownTipsyFields
    _particle_mass_name = "Mass"
    _particle_coordinates_name = "Coordinates"
    _header_spec = (('time',    'd'),
                    ('nbodies', 'i'),
                    ('ndim',    'i'),
                    ('nsph',    'i'),
                    ('ndark',   'i'),
                    ('nstar',   'i'),
                    ('dummy',   'i'))

    def __init__(self, filename, data_style="tipsy",
                 endian=">",
                 field_dtypes=None,
                 domain_left_edge=None,
                 domain_right_edge=None,
                 unit_base=None,
                 cosmology_parameters=None,
                 parameter_file=None,
                 n_ref=64, over_refine_factor=1):
        self.n_ref = n_ref
        self.over_refine_factor = over_refine_factor
        self.endian = endian
        self.storage_filename = None
        if domain_left_edge is None:
            domain_left_edge = np.zeros(3, "float64") - 0.5
        if domain_right_edge is None:
            domain_right_edge = np.zeros(3, "float64") + 0.5

        self.domain_left_edge = np.array(domain_left_edge, dtype="float64")
        self.domain_right_edge = np.array(domain_right_edge, dtype="float64")

        # My understanding is that dtypes are set on a field by field basis,
        # not on a (particle type, field) basis
        if field_dtypes is None:
            field_dtypes = {}
        self._field_dtypes = field_dtypes

        self._unit_base = unit_base or {}
        self._cosmology_parameters = cosmology_parameters
        if parameter_file is not None:
            parameter_file = os.path.abspath(parameter_file)
        self._param_file = parameter_file
        filename = os.path.abspath(filename)
        super(TipsyStaticOutput, self).__init__(filename, data_style)

    def __repr__(self):
        return os.path.basename(self.parameter_filename)

    def _parse_parameter_file(self):

        # Parsing the header of the tipsy file, from this we obtain
        # the snapshot time and particle counts.

        f = open(self.parameter_filename, "rb")
        hh = self.endian + "".join(["%s" % (b) for a, b in self._header_spec])
        hvals = dict([(a, c) for (a, b), c in zip(self._header_spec,
                     struct.unpack(hh, f.read(struct.calcsize(hh))))])
        self.parameters.update(hvals)
        self._header_offset = f.tell()

        # These are always true, for now.
        self.dimensionality = 3
        self.refine_by = 2
        self.parameters["HydroMethod"] = "sph"

        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])

        # Read in parameter file, if available.
        if self._param_file is None:
            pfn = glob.glob(os.path.join(self.directory, "*.param"))
            assert len(pfn) < 2, \
                "More than one param file is in the data directory"
            if pfn == []:
                pfn = None
            else:
                pfn = pfn[0]
        else:
            pfn = self._param_file

        if pfn is not None:
            for line in (l.strip() for l in open(pfn)):
                # skip comment lines and blank lines
                l = line.strip()
                if l.startswith('#') or l == '':
                    continue
                # parse parameters according to tipsy parameter type
                param, val = (i.strip() for i in line.split('=', 1))
                if param.startswith('n') or param.startswith('i'):
                    val = long(val)
                elif param.startswith('d'):
                    val = float(val)
                elif param.startswith('b'):
                    val = bool(float(val))
                self.parameters[param] = val

        self.current_time = hvals["time"]
        nz = 1 << self.over_refine_factor
        self.domain_dimensions = np.ones(3, "int32") * nz
        if self.parameters.get('bPeriodic', True):
            self.periodicity = (True, True, True)
        else:
            self.periodicity = (False, False, False)

        if self.parameters.get('bComove', True):
            self.cosmological_simulation = 1
            cosm = self._cosmology_parameters or {}
            dcosm = dict(current_redshift=0.0,
                         omega_lambda=0.0,
                         omega_matter=0.0,
                         hubble_constant=1.0)
            for param in ['current_redshift', 'omega_lambda',
                          'omega_matter', 'hubble_constant']:
                pval = cosm.get(param, dcosm[param])
                setattr(self, param, pval)
        else:
            self.cosmological_simulation = 0.0
            kpc_unit = self.parameters.get('dKpcUnit', 1.0)
            self._unit_base['cm'] = 1.0 / (kpc_unit * cm_per_kpc)

        self.filename_template = self.parameter_filename
        self.file_count = 1

        f.close()

    def _set_units(self):
        super(TipsyStaticOutput, self)._set_units()
        if self.cosmological_simulation:
            DW = (self.domain_right_edge - self.domain_left_edge).max()
            cosmo = Cosmology(self.hubble_constant * 100.0,
                              self.omega_matter, self.omega_lambda)
            length_unit = DW * self.units['cm']  # Get it in proper cm
            density_unit = cosmo.CriticalDensity(self.current_redshift)
            mass_unit = density_unit * length_unit ** 3
        else:
            mass_unit = self.parameters.get('dMsolUnit', 1.0) * mass_sun_cgs
            length_unit = self.parameters.get('dKpcUnit', 1.0) * cm_per_kpc
            density_unit = mass_unit / length_unit ** 3
        time_unit = 1.0 / np.sqrt(G * density_unit)
        velocity_unit = length_unit / time_unit
        self.conversion_factors["velocity"] = velocity_unit
        self.conversion_factors["mass"] = mass_unit
        self.conversion_factors["density"] = density_unit
        for u in sec_conversion:
            self.time_units[u] = time_unit * sec_conversion[u]

    @classmethod
    def _is_valid(self, *args, **kwargs):
        # We do not allow load() of these files.
        return False
