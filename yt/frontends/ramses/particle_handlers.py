import os
import yt.utilities.fortran_utils as fpu
import glob
from yt.funcs import mylog

from .io import _read_part_file_descriptor

PARTICLE_HANDLERS = []

def get_particle_handlers():
    return PARTICLE_HANDLERS

def register_particle_handler(ph):
    PARTICLE_HANDLERS.append(ph)


class ParticleFileHandler(object):
    '''Abstract class to handle particles in RAMSES.


    See `SinkParticleFileHandler` for an example implementation.'''

    ptype = None  # The name to give to the particle type
    fname = None  # The name of the file(s).

    attrs = None  # The attributes of the header
    known_fields = None  # A list of tuple containing the field name and its type

    @classmethod
    def any_exist(cls, ds):
        '''Return True if any file of this type is found

        Arguments
        ---------
        * ds: a Ramses Dataset
        '''
        raise NotImplementedError


    def __init__(self, ds, domain_id):
        self.ds = ds
        basename = os.path.abspath(
              os.path.dirname(ds.parameter_filename))
        iout = int(
            os.path.basename(ds.parameter_filename)
            .split(".")[0].
            split("_")[1])
        icpu = domain_id
        self.fname = os.path.join(
            basename,
            self.fname.format(iout=iout, icpu=icpu))

    @property
    def exists(self):
        '''Return True if the fname exists'''
        return os.path.exists(self.fname)

    def read_header(self):
        '''
        This function should read the header, compute the offsets of
        the file into self.offsets and store the fields found in
        self.fields.'''
        raise NotImplementedError

class DefaultParticleFileHandler(ParticleFileHandler):
    ptype = 'io'
    fname = 'part_{iout:05d}.out{icpu:05d}'
    file_descriptor = 'part_file_descriptor.txt'

    attrs = ( ('ncpu', 1, 'I'),
              ('ndim', 1, 'I'),
              ('npart', 1, 'I') )

    known_fields = [
        ("particle_position_x", "d"),
        ("particle_position_y", "d"),
        ("particle_position_z", "d"),
        ("particle_velocity_x", "d"),
        ("particle_velocity_y", "d"),
        ("particle_velocity_z", "d"),
        ("particle_mass", "d"),
        ("particle_identifier", "i"),
        ("particle_refinement_level", "I")]


    @classmethod
    def any_exist(cls, ds):
        files = os.path.join(
            os.path.split(ds.parameter_filename)[0],
            'part_?????.out?????')
        ret = len(glob.glob(files)) > 0
        return ret

    def read_header(self):
        f = open(self.fname, "rb")
        f.seek(0, os.SEEK_END)
        flen = f.tell()
        f.seek(0)
        hvals = {}
        attrs = ( ('ncpu', 1, 'I'),
                  ('ndim', 1, 'I'),
                  ('npart', 1, 'I') )
        hvals.update(fpu.read_attrs(f, attrs))
        fpu.read_vector(f, 'I')

        attrs = ( ('nstar_tot', 1, 'I'),
                  ('mstar_tot', 1, 'd'),
                  ('mstar_lost', 1, 'd'),
                  ('nsink', 1, 'I') )
        hvals.update(fpu.read_attrs(f, attrs))
        self.header = hvals
        self.local_particle_count = hvals['npart']

        if self.has_part_descriptor:
            particle_fields = (
                _read_part_file_descriptor(self.file_descriptor)
            )
        else:
            particle_fields = self.known_fields.copy()

            if self.ds._extra_particle_fields is not None:
                particle_fields += self.ds._extra_particle_fields

        field_offsets = {}
        _pfields = {}

        ptype = self.ptype

        # Read offsets
        for field, vtype in particle_fields:
            if f.tell() >= flen: break
            field_offsets[ptype, field] = f.tell()
            _pfields[ptype, field] = vtype
            fpu.skip(f, 1)

        iextra = 0
        while f.tell() < flen:
            iextra += 1
            field, vtype = ('particle_extra_field_%i' % iextra, 'd')
            particle_fields.append((field, vtype))

            field_offsets[ptype, field] = f.tell()
            _pfields[ptype, field] = vtype
            fpu.skip(f, 1)

        if iextra > 0 and not self.ds._warn_extra_fields:
            self.ds._warn_extra_fields = True
            w = ("Detected %s extra particle fields assuming kind "
                 "`double`. Consider using the `extra_particle_fields` "
                 "keyword argument if you have unexpected behavior.")
            mylog.warning(w % iextra)

        self.field_offsets = field_offsets
        self.field_types = _pfields

    @property
    def has_part_descriptor(self):
        '''
        Does the output include particle file descriptor?
        '''
        return os.path.exists(self.file_descriptor)



class SinkParticleFileHandler(ParticleFileHandler):
    '''Handle sink files'''
    ptype = 'sink'
    fname = 'sink_{iout:05d}.out{icpu:05d}'

    attrs = (('nsink', 1, 'I'),
             ('nindsink', 1, 'I'))

    known_fields = [
        ("particle_identifier", "i"),
        ("particle_mass", "d"),
        ("particle_position_x", "d"),
        ("particle_position_y", "d"),
        ("particle_position_z", "d"),
        ("particle_velocity_x", "d"),
        ("particle_velocity_y", "d"),
        ("particle_velocity_z", "d"),
        ("particle_age", "d"),
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
        ("BH_efficiency", "d")]

    @classmethod
    def any_exist(cls, ds):
        files = os.path.join(
            os.path.split(ds.parameter_filename)[0],
            'sink_?????.out?????')
        ret = len(glob.glob(files)) > 0
        return ret

    def read_header(self):
        f = open(self.fname, "rb")
        f.seek(0, os.SEEK_END)
        flen = f.tell()
        f.seek(0)
        hvals = {}
        # Read the header of the file
        attrs = self.attrs
        hvals.update(fpu.read_attrs(f, attrs))
        self._header = hvals
        self._sink_count = hvals['nsink']

        # Read the fields + add the sink properties
        fields = self.known_fields.copy()
        for i in range(self.ds.dimensionality*2+1):
            for j in range(self.ds.max_level, self.ds.min_level):
                fields.append((
                    "particle_prop_%s_%s" % (i, j), "d"
                ))

        field_offsets = {}
        _pfields = {}

        # Fill the fields, offsets and types
        self.fields = []
        for field, vtype in fields:
            self.fields.append(field)
            if f.tell() >= flen: break
            field_offsets[self.ptype, field] = f.tell()
            _pfields[self.ptype, field] = vtype
            fpu.skip(f, 1)
        self.field_offsets = field_offsets
        self.field_types = _pfields

register_particle_handler(DefaultParticleFileHandler)
register_particle_handler(SinkParticleFileHandler)
