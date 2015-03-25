"""
Data structures for Tipsy frontend




"""
from __future__ import print_function

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
import stat
import struct
import glob
import os

from yt.frontends.sph.data_structures import \
    ParticleDataset
from yt.geometry.particle_geometry_handler import \
    ParticleIndex
from yt.data_objects.static_output import \
    ParticleFile
from yt.utilities.cosmology import \
    Cosmology
from yt.utilities.physical_constants import \
    G, \
    cm_per_kpc

from .fields import \
    TipsyFieldInfo

class TipsyFile(ParticleFile):
    def __init__(self, ds, io, filename, file_id):
        # To go above 1 domain, we need to include an indexing step in the
        # IOHandler, rather than simply reading from a single file.
        assert file_id == 0
        super(TipsyFile, self).__init__(ds, io, filename, file_id)
        io._create_dtypes(self)
        io._update_domain(self)#Check automatically what the domain size is

    def _calculate_offsets(self, field_list):
        self.field_offsets = self.io._calculate_particle_offsets(self)

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
                 n_ref=64, over_refine_factor=1,
                 bounding_box = None,
                 units_override=None):
        self.n_ref = n_ref
        self.over_refine_factor = over_refine_factor
        if field_dtypes is None:
            field_dtypes = {}
        success, self.endian = self._validate_header(filename)
        if not success:
            print("SOMETHING HAS GONE WRONG.  NBODIES != SUM PARTICLES.")
            print("%s != (%s == %s + %s + %s)" % (
                self.parameters['nbodies'],
                tot,
                self.parameters['nsph'],
                self.parameters['ndark'],
                self.parameters['nstar']))
            print("Often this can be fixed by changing the 'endian' parameter.")
            print("This defaults to '>' but may in fact be '<'.")
            raise RuntimeError
        self.storage_filename = None

        # My understanding is that dtypes are set on a field by field basis,
        # not on a (particle type, field) basis
        self._field_dtypes = field_dtypes

        self._unit_base = unit_base or {}

        if bounding_box is not None:
            bbox = np.array(bounding_box, dtype="float64")
            if bbox.shape == (2, 3):
                bbox = bbox.transpose()
            self.domain_left_edge = bbox[:,0]
            self.domain_right_edge = bbox[:,1]
        else:
            self.domain_left_edge = self.domain_right_edge = None
 

        self._cosmology_parameters = cosmology_parameters
        if parameter_file is not None:
            parameter_file = os.path.abspath(parameter_file)
        self._param_file = parameter_file
        filename = os.path.abspath(filename)
        if units_override is not None:
            raise RuntimeError("units_override is not supported for TipsyDataset. "+
                               "Use unit_base instead.")
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
