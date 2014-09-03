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
import time
import os
import types

from yt.utilities.fortran_utils import read_record
from yt.utilities.logger import ytLogger as mylog
from yt.geometry.particle_geometry_handler import \
    ParticleIndex
from yt.data_objects.static_output import \
    Dataset, ParticleFile
from yt.utilities.definitions import \
    mpc_conversion, sec_conversion
from yt.utilities.physical_constants import \
    G, \
    cm_per_kpc, \
    mass_sun_cgs
from yt.utilities.cosmology import Cosmology
from .fields import \
    SPHFieldInfo, OWLSFieldInfo, TipsyFieldInfo, EagleNetworkFieldInfo
from .definitions import \
    gadget_header_specs, \
    gadget_field_specs, \
    gadget_ptype_specs
from .io import \
    IOHandlerTipsyBinary

try:
    import requests
    import json
except ImportError:
    requests = None

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


class ParticleDataset(Dataset):
    _unit_base = None
    over_refine_factor = 1
    filter_bbox = False


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
                 ptype_spec = "default"):
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

        prefix = self.parameter_filename.split(".", 1)[0]

        if hvals["NumFiles"] > 1:
            self.filename_template = "%s.%%(num)s%s" % (prefix, self._suffix)
        else:
            self.filename_template = self.parameter_filename

        self.file_count = hvals["NumFiles"]

    def _set_code_unit_attributes(self):
        # Set a sane default for cosmological simulations.
        if self._unit_base is None and self.cosmological_simulation == 1:
            mylog.info("Assuming length units are in Mpc/h (comoving)")
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
                 bounding_box = None):
        self.storage_filename = None
        filename = os.path.abspath(filename)
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

class OWLSDataset(GadgetHDF5Dataset):
    _particle_mass_name = "Mass"
    _field_info_class = OWLSFieldInfo
    _time_readin = "Time_GYR"



    def _parse_parameter_file(self):
        handle = h5py.File(self.parameter_filename, mode="r")
        hvals = {}
        hvals.update((str(k), v) for k, v in handle["/Header"].attrs.items())
        hvals["NumFiles"] = hvals["NumFilesPerSnapshot"]
        hvals["Massarr"] = hvals["MassTable"]

        self.dimensionality = 3
        self.refine_by = 2
        self.parameters["HydroMethod"] = "sph"
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])

        # Set standard values
        self.current_time = hvals[self._time_readin] * sec_conversion["Gyr"]
        if self.domain_left_edge is None:
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
               "Header" in fileh["/"].keys() and \
               "SUBFIND" not in fileh["/"].keys() and not\
               ("ChemistryAbundances" in fileh["PartType0"].keys()
                or "ChemicalAbundances" in fileh["PartType0"].keys()):
                fileh.close()
                return True
            fileh.close()
        except:
            pass
        return False

class EagleNetworkDataset(OWLSDataset):
    _particle_mass_name = "Mass"
    _field_info_class = EagleNetworkFieldInfo
    _time_readin = 'Time'

    @classmethod
    def _is_valid(self, *args, **kwargs):
        try:
            fileh = h5py.File(args[0], mode='r')
            if "Constants" in fileh["/"].keys() and \
               "Header" in fileh["/"].keys() and \
               "SUBFIND" not in fileh["/"].keys() and \
               ("ChemistryAbundances" in fileh["PartType0"].keys()
                or "ChemicalAbundances" in fileh["PartType0"].keys()):
                fileh.close()
                return True
            fileh.close()
        except:
            pass
        return False

class TipsyFile(ParticleFile):

    def _calculate_offsets(self, field_list):
        self.field_offsets = self.io._calculate_particle_offsets(self)

    def __init__(self, ds, io, filename, file_id):
        # To go above 1 domain, we need to include an indexing step in the
        # IOHandler, rather than simply reading from a single file.
        assert file_id == 0
        super(TipsyFile, self).__init__(ds, io, filename, file_id)
        io._create_dtypes(self)
        io._update_domain(self)#Check automatically what the domain size is


class TipsyDataset(ParticleDataset):
    _index_class = ParticleIndex
    _file_class = TipsyFile
    _field_info_class = TipsyFieldInfo
    _particle_mass_name = "Mass"
    _particle_coordinates_name = "Coordinates"
    _header_spec = (('time',    'd'),
                    ('nbodies', 'i'),
                    ('ndim',    'i'),
                    ('nsph',    'i'),
                    ('ndark',   'i'),
                    ('nstar',   'i'),
                    ('dummy',   'i'))

    def __init__(self, filename, dataset_type="tipsy",
                 field_dtypes=None,
                 unit_base=None,
                 parameter_file=None,
                 cosmology_parameters=None,
                 n_ref=64, over_refine_factor=1):
        self.n_ref = n_ref
        self.over_refine_factor = over_refine_factor
        if field_dtypes is None:
            field_dtypes = {}
        success, self.endian = self._validate_header(filename)
        if not success:
            print "SOMETHING HAS GONE WRONG.  NBODIES != SUM PARTICLES."
            print "%s != (%s == %s + %s + %s)" % (
                self.parameters['nbodies'],
                tot,
                self.parameters['nsph'],
                self.parameters['ndark'],
                self.parameters['nstar'])
            print "Often this can be fixed by changing the 'endian' parameter."
            print "This defaults to '>' but may in fact be '<'."
            raise RuntimeError
        self.storage_filename = None

        # My understanding is that dtypes are set on a field by field basis,
        # not on a (particle type, field) basis
        self._field_dtypes = field_dtypes

        self._unit_base = unit_base or {}
        self._cosmology_parameters = cosmology_parameters
        if parameter_file is not None:
            parameter_file = os.path.abspath(parameter_file)
        self._param_file = parameter_file
        filename = os.path.abspath(filename)
        super(TipsyDataset, self).__init__(filename, dataset_type)

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
                val = val.split('#')[0]
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
        periodic = self.parameters.get('bPeriodic', True)
        period = self.parameters.get('dPeriod', None)
        comoving = self.parameters.get('bComove', False)
        self.periodicity = (periodic, periodic, periodic)
        if comoving and period is None:
            period = 1.0
        if periodic and period is not None:
            # If we are periodic, that sets our domain width to either 1 or dPeriod.
            self.domain_left_edge = np.zeros(3, "float64") - 0.5*period
            self.domain_right_edge = np.zeros(3, "float64") + 0.5*period
        else:
            self.domain_left_edge = None
            self.domain_right_edge = None
        if comoving:
            cosm = self._cosmology_parameters or {}
            self.scale_factor = hvals["time"]#In comoving simulations, time stores the scale factor a
            self.cosmological_simulation = 1
            dcosm = dict(current_redshift=(1.0/self.scale_factor)-1.0,
                         omega_lambda=self.parameters.get('dLambda', cosm.get('omega_lambda',0.0)),
                         omega_matter=self.parameters.get('dOmega0', cosm.get('omega_matter',0.0)),
                         hubble_constant=self.parameters.get('dHubble0', cosm.get('hubble_constant',1.0)))
            for param in dcosm.keys():
                pval = dcosm[param]
                setattr(self, param, pval)
        else:
            self.cosmological_simulation = 0.0
            kpc_unit = self.parameters.get('dKpcUnit', 1.0)
            self._unit_base['cm'] = 1.0 / (kpc_unit * cm_per_kpc)

        self.filename_template = self.parameter_filename
        self.file_count = 1

        f.close()

    def _set_derived_attrs(self):
        if self.domain_left_edge is None or self.domain_right_edge is None:
            self.domain_left_edge = np.nan
            self.domain_right_edge = np.nan
            self.index
        super(TipsyDataset, self)._set_derived_attrs()

    def _set_code_unit_attributes(self):
        if self.cosmological_simulation:
            mu = self.parameters.get('dMsolUnit', 1.)
            lu = self.parameters.get('dKpcUnit', 1000.)
            # In cosmological runs, lengths are stored as length*scale_factor
            self.length_unit = self.quan(lu, 'kpc')*self.scale_factor
            self.mass_unit = self.quan(mu, 'Msun')
            density_unit = self.mass_unit/ (self.length_unit/self.scale_factor)**3
            # Gasoline's hubble constant, dHubble0, is stored units of proper code time.
            self.hubble_constant *= np.sqrt(G.in_units('kpc**3*Msun**-1*s**-2')*density_unit).value/(3.2407793e-18)  
            cosmo = Cosmology(self.hubble_constant,
                              self.omega_matter, self.omega_lambda)
            self.current_time = cosmo.hubble_time(self.current_redshift)
        else:
            mu = self.parameters.get('dMsolUnit', 1.0)
            self.mass_unit = self.quan(mu, 'Msun')
            lu = self.parameters.get('dKpcUnit', 1.0)
            self.length_unit = self.quan(lu, 'kpc')
            density_unit = self.mass_unit / self.length_unit**3
        self.time_unit = 1.0 / np.sqrt(G * density_unit)

    @staticmethod
    def _validate_header(filename):
        '''
        This method automatically detects whether the tipsy file is big/little endian
        and is not corrupt/invalid.  It returns a tuple of (Valid, endianswap) where
        Valid is a boolean that is true if the file is a tipsy file, and endianswap is 
        the endianness character '>' or '<'.
        '''
        try:
            f = open(filename,'rb')
        except:
            return False, 1
        try:
            f.seek(0, os.SEEK_END)
            fs = f.tell()
            f.seek(0, os.SEEK_SET)
            #Read in the header
            t, n, ndim, ng, nd, ns = struct.unpack("<diiiii", f.read(28))
        except IOError:
            return False, 1
        endianswap = "<"
        #Check Endianness
        if (ndim < 1 or ndim > 3):
            endianswap = ">"
            f.seek(0)
            t, n, ndim, ng, nd, ns = struct.unpack(">diiiii", f.read(28))
        # File is borked if this is true.  The header is 28 bytes, and may
        # Be followed by a 4 byte pad.  Next comes gas particles, which use
        # 48 bytes, followed by 36 bytes per dark matter particle, and 44 bytes
        # per star particle.  If positions are stored as doubles, each of these
        # sizes is increased by 12 bytes.
        if (fs != 28+48*ng+36*nd+44*ns and fs != 28+60*ng+48*nd+56*ns and
                fs != 32+48*ng+36*nd+44*ns and fs != 32+60*ng+48*nd+56*ns):
            f.close()
            return False, 0
        f.close()
        return True, endianswap

    @classmethod
    def _is_valid(self, *args, **kwargs):
        return TipsyDataset._validate_header(args[0])[0]

class HTTPParticleFile(ParticleFile):
    pass

class HTTPStreamDataset(ParticleDataset):
    _index_class = ParticleIndex
    _file_class = HTTPParticleFile
    _field_info_class = SPHFieldInfo
    _particle_mass_name = "Mass"
    _particle_coordinates_name = "Coordinates"
    _particle_velocity_name = "Velocities"
    filename_template = ""
    
    def __init__(self, base_url,
                 dataset_type = "http_particle_stream",
                 n_ref = 64, over_refine_factor=1):
        if requests is None:
            raise RuntimeError
        self.base_url = base_url
        self.n_ref = n_ref
        self.over_refine_factor = over_refine_factor
        super(HTTPStreamDataset, self).__init__("", dataset_type)

    def __repr__(self):
        return self.base_url

    def _parse_parameter_file(self):
        self.dimensionality = 3
        self.refine_by = 2
        self.parameters["HydroMethod"] = "sph"

        # Here's where we're going to grab the JSON index file
        hreq = requests.get(self.base_url + "/yt_index.json")
        if hreq.status_code != 200:
            raise RuntimeError
        header = json.loads(hreq.content)
        header['particle_count'] = dict((int(k), header['particle_count'][k])
            for k in header['particle_count'])
        self.parameters = header

        # Now we get what we need
        self.domain_left_edge = np.array(header['domain_left_edge'], "float64")
        self.domain_right_edge = np.array(header['domain_right_edge'], "float64")
        nz = 1 << self.over_refine_factor
        self.domain_dimensions = np.ones(3, "int32") * nz
        self.periodicity = (True, True, True)

        self.current_time = header['current_time']
        self.unique_identifier = header.get("unique_identifier", time.time())
        self.cosmological_simulation = int(header['cosmological_simulation'])
        for attr in ('current_redshift', 'omega_lambda', 'omega_matter',
                     'hubble_constant'):
            setattr(self, attr, float(header[attr]))

        self.file_count = header['num_files']

    def _set_units(self):
        length_unit = float(self.parameters['units']['length'])
        time_unit = float(self.parameters['units']['time'])
        mass_unit = float(self.parameters['units']['mass'])
        density_unit = mass_unit / length_unit ** 3
        velocity_unit = length_unit / time_unit
        self._unit_base = {}
        self._unit_base['cm'] = 1.0/length_unit
        self._unit_base['s'] = 1.0/time_unit
        super(HTTPStreamDataset, self)._set_units()
        self.conversion_factors["velocity"] = velocity_unit
        self.conversion_factors["mass"] = mass_unit
        self.conversion_factors["density"] = density_unit

    @classmethod
    def _is_valid(self, *args, **kwargs):
        if not args[0].startswith("http://"):
            return False
        hreq = requests.get(args[0] + "/yt_index.json")
        if hreq.status_code == 200:
            return True
        return False
