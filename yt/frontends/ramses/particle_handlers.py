import os
import yt.utilities.fortran_utils as fpu
import glob
from yt.extern.six import add_metaclass
from yt.funcs import mylog
from yt.config import ytcfg

from .io import _read_part_file_descriptor

PARTICLE_HANDLERS = set()

def get_particle_handlers():
    return PARTICLE_HANDLERS

def register_particle_handler(ph):
    PARTICLE_HANDLERS.add(ph)


class RAMSESParticleFileHandlerRegistry(type):
    """
    This is a base class that on instantiation registers the file
    handler into the list. Used as a metaclass.
    """
    def __new__(meta, name, bases, class_dict):
        cls = type.__new__(meta, name, bases, class_dict)
        if cls.ptype is not None:
            register_particle_handler(cls)
        return cls


@add_metaclass(RAMSESParticleFileHandlerRegistry)
class ParticleFileHandler(object):
    '''
    Abstract class to handle particles in RAMSES. Each instance
    represents a single file (one domain).

    To add support to a new particle file, inherit from this class and
    implement all functions containing a `NotImplementedError`.

    See `SinkParticleFileHandler` for an example implementation.'''

    # These properties are static properties
    ptype = None  # The name to give to the particle type
    fname = None  # The name of the file(s).
    file_descriptor = None # The name of the file descriptor (if any)

    attrs = None  # The attributes of the header
    known_fields = None  # A list of tuple containing the field name and its type
    config_field = None  # Name of the config section (if any)

    # These properties are computed dynamically
    field_offsets = None     # Mapping from field to offset in file
    field_types = None       # Mapping from field to the type of the data (float, integer, ...)
    local_particle_count = None  # The number of particle in the domain


    def __init__(self, ds, domain_id):
        '''
        Initalize an instance of the class. This automatically sets
        the full path to the file. This is not intended to be
        overriden in most cases.

        If you need more flexibility, rewrite this function to your
        need in the inherited class.
        '''
        self.ds = ds
        self.domain_id = domain_id
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

        if self.file_descriptor is not None:
            self.file_descriptor = os.path.join(
                basename,
                self.file_descriptor)

        # Attempt to read the list of fields from the config file
        if self.config_field and ytcfg.has_section(self.config_field):
            cfg = ytcfg.get(self.config_field, 'fields')
            known_fields = []
            for c in (_.strip() for _ in cfg.split('\n') if _.strip() != ''):
                field, field_type = (_.strip() for _ in c.split(','))
                known_fields.append((field, field_type))
            self.known_fields = known_fields

    @property
    def exists(self):
        '''
        This function should return True if the *file* the instance
        exists. It is called for each file of the type found on the
        disk.

        By default, it just returns whether the file exists. Override
        it for more complex cases.
        '''
        return os.path.exists(self.fname)

    @property
    def has_part_descriptor(self):
        '''
        This function should return True if a *file descriptor*
        exists.

        By default, it just returns whether the file exists. Override
        it for more complex cases.
        '''
        return os.path.exists(self.file_descriptor)

    @classmethod
    def any_exist(cls, ds):
        '''
        This function should return True if the kind of particle
        represented by the class exists in the dataset. It takes as
        argument the class itself -not an instance- and a dataset.

        Arguments
        ---------
        * ds: a Ramses Dataset

        Note
        ----
        This function is usually called once at the initialization of
        the RAMSES Dataset structure to determine if the particle type
        (e.g. regular particles) exists.
        '''
        # this function must be implemented by subclasses
        raise NotImplementedError


    def read_header(self):
        '''
        This function is called once per file. It should:
        * read the header of the file and store any relevant information
        * detect the fields in the file
        * compute the offsets (location in the file) of each field

        It is in charge of setting `self.field_offsets` and `self.field_types`.
        * `field_offsets`: dictionary: tuple -> integer
           A dictionary that maps `(type, field_name)` to their
           location in the file (integer)
        * `field_types`: dictionary: tuple -> character
           A dictionary that maps `(type, field_name)` to their type
           (character), following Python's struct convention.
        '''
        # this function must be implemented by subclasses
        raise NotImplementedError


class DefaultParticleFileHandler(ParticleFileHandler):
    ptype = 'io'
    fname = 'part_{iout:05d}.out{icpu:05d}'
    file_descriptor = 'part_file_descriptor.txt'
    config_field = 'ramses-particles'

    attrs = ( ('ncpu', 1, 'I'),
              ('ndim', 1, 'I'),
              ('npart', 1, 'I'),
              ('localseed', 4, 'I'),
              ('nstar_tot', 1, 'I'),
              ('mstar_tot', 1, 'd'),
              ('mstar_lost', 1, 'd'),
              ('nsink', 1, 'I') )

    known_fields = [
        ("particle_position_x", "d"),
        ("particle_position_y", "d"),
        ("particle_position_z", "d"),
        ("particle_velocity_x", "d"),
        ("particle_velocity_y", "d"),
        ("particle_velocity_z", "d"),
        ("particle_mass", "d"),
        ("particle_identity", "i"),
        ("particle_refinement_level", "I")]

    @classmethod
    def any_exist(cls, ds):
        files = os.path.join(
            os.path.split(ds.parameter_filename)[0],
            'part_?????.out?????')
        ret = len(glob.glob(files)) > 0
        return ret

    def read_header(self):
        if not self.exists:
            self.field_offsets = {}
            self.field_types = {}
            self.local_particle_count = 0
            return
        f = open(self.fname, "rb")
        f.seek(0, os.SEEK_END)
        flen = f.tell()
        f.seek(0)
        hvals = {}
        attrs = self.attrs
        hvals.update(fpu.read_attrs(f, attrs))
        self.header = hvals
        self.local_particle_count = hvals['npart']
        extra_particle_fields = self.ds._extra_particle_fields

        if self.has_part_descriptor:
            particle_fields = (
                _read_part_file_descriptor(self.file_descriptor)
            )
        else:
            particle_fields = list(self.known_fields)

            if extra_particle_fields is not None:
                particle_fields += extra_particle_fields

        if hvals["nstar_tot"] > 0 and extra_particle_fields is not None:
            particle_fields += [("particle_birth_time", "d"),
                                ("particle_metallicity", "d")]

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

        if iextra > 0 and not self.ds._warned_extra_fields['io']:
            w = ("Detected %s extra particle fields assuming kind "
                 "`double`. Consider using the `extra_particle_fields` "
                 "keyword argument if you have unexpected behavior.")
            mylog.warning(w % iextra)
            self.ds._warned_extra_fields['io'] = True

        self.field_offsets = field_offsets
        self.field_types = _pfields


class SinkParticleFileHandler(ParticleFileHandler):
    '''Handle sink files'''
    ptype = 'sink'
    fname = 'sink_{iout:05d}.out{icpu:05d}'
    file_descriptor = 'sink_file_descriptor.txt'
    config_field = 'ramses-sink-particles'

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
        ("BH_efficiency", "d")]

    @classmethod
    def any_exist(cls, ds):
        files = os.path.join(
            os.path.split(ds.parameter_filename)[0],
            'sink_?????.out?????')
        ret = len(glob.glob(files)) > 0
        return ret

    def read_header(self):
        if not self.exists:
            self.field_offsets = {}
            self.field_types = {}
            self.local_particle_count = 0
            return
        f = open(self.fname, "rb")
        f.seek(0, os.SEEK_END)
        flen = f.tell()
        f.seek(0)
        hvals = {}
        # Read the header of the file
        attrs = self.attrs

        hvals.update(fpu.read_attrs(f, attrs))
        self._header = hvals

        # This is somehow a trick here: we only want one domain to
        # be read, as ramses writes all the sinks in all the
        # domains. Here, we set the local_particle_count to 0 except
        # for the first domain to be red.
        if getattr(self.ds, '_sink_file_flag', False):
            self.local_particle_count = 0
        else:
            self.ds._sink_file_flag = True
            self.local_particle_count = hvals['nsink']

        # Read the fields + add the sink properties
        if self.has_part_descriptor:
            fields = (
                _read_part_file_descriptor(self.file_descriptor)
            )
        else:
            fields = list(self.known_fields)

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
