import os
import yt.utilities.fortran_utils as fpu
import glob
from yt.extern.six import add_metaclass
from yt.funcs import mylog
import numpy as np

FIELD_HANDLERS = set()

def get_field_handlers():
    return FIELD_HANDLERS

def register_field_handler(ph):
    FIELD_HANDLERS.add(ph)


class RAMSESFieldFileHandlerRegister(type):
    """
    This is a base class that on instantiation registers the file
    handler into the list. Used as a metaclass.
    """
    def __new__(meta, name, bases, class_dict):
        cls = type.__new__(meta, name, bases, class_dict)
        if cls.ftype is not None:
            register_field_handler(cls)
        return cls


@add_metaclass(RAMSESFieldFileHandlerRegister)
class FieldFileHandler(object):
    '''
    Abstract class to handle particles in RAMSES. Each instance
    represents a single file (one domain).

    To add support to a new particle file, inherit from this class and
    implement all functions containing a `NotImplementedError`.

    See `SinkParticleFileHandler` for an example implementation.'''

    # These properties are static properties
    ftype = None  # The name to give to the field type
    fname = None  # The name of the file(s).
    attrs = None  # The attributes of the header
    known_fields = None  # A list of tuple containing the field name and its type

    # These properties are computed dynamically
    field_offsets = None     # Mapping from field to offset in file
    field_types = None       # Mapping from field to the type of the data (float, integer, …)

    def __init__(self, domain):
        '''
        Initalize an instance of the class. This automatically sets
        the full path to the file. This is not intended to be
        overriden in most cases.

        If you need more flexibility, rewrite this function to your
        need in the inherited class.
        '''
        self.domain = domain
        self.domain_id = domain.domain_id
        ds = domain.ds
        basename = os.path.abspath(
              os.path.dirname(ds.parameter_filename))
        iout = int(
            os.path.basename(ds.parameter_filename)
            .split(".")[0].
            split("_")[1])

        self.fname = os.path.join(
            basename,
            self.fname.format(iout=iout, icpu=domain.domain_id))

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

    @classmethod
    def any_exist(cls, ds):
        '''
        This function should return True if the kind of particle
        represented by the class exists in the dataset. It takes as
        argument the class itself —not an instance— and a dataset.

        Arguments
        ---------
        * ds: a Ramses Dataset

        Note
        ----
        This function is usually called once at the initialization of
        the RAMSES Dataset structure to determine if the particle type
        (e.g. regular particles) exists.
        '''
        raise NotImplementedError

    @classmethod
    def detect_fields(cls, ds):
        raise NotImplementedError


    @property
    def offset(self):
        raise NotImplementedError

    @property
    def level_count(self):
        raise NotImplementedError


class HydroFieldFileHandler(FieldFileHandler):
    ftype = 'ramses'
    fname = 'hydro_{iout:05d}.out{icpu:05d}'
    attrs = ( ('ncpu', 1, 'i'),
              ('nvar', 1, 'i'),
              ('ndim', 1, 'i'),
              ('nlevelmax', 1, 'i'),
              ('nboundary', 1, 'i'),
              ('gamma', 1, 'd'))

    @classmethod
    def any_exist(cls, ds):
        files = os.path.join(
            os.path.split(ds.parameter_filename)[0],
            'hydro_?????.out?????')
        ret = len(glob.glob(files)) > 0
        return ret

    @classmethod
    def detect_fields(cls, ds):
        if getattr(cls, 'field_list', None) is not None:
            return cls.field_list

        num = os.path.basename(ds.parameter_filename).split("."
                )[0].split("_")[1]
        testdomain = 1 # Just pick the first domain file to read
        basename = "%s/%%s_%s.out%05i" % (
            os.path.abspath(
              os.path.dirname(ds.parameter_filename)),
            num, testdomain)
        fname = basename % "hydro"

        f = open(fname, 'rb')
        attrs = cls.attrs
        hvals = fpu.read_attrs(f, attrs)

        # Store some metadata
        ds.gamma = hvals['gamma']
        nvar = cls.nvar = hvals['nvar']

        if ds._fields_in_file is not None:
            fields = [('ramses', f) for f in ds._fields_in_file]
        else:
            foldername  = os.path.abspath(os.path.dirname(ds.parameter_filename))
            rt_flag = any(glob.glob(os.sep.join([foldername, 'info_rt_*.txt'])))
            if rt_flag: # rt run
                if nvar < 10:
                    mylog.info('Detected RAMSES-RT file WITHOUT IR trapping.')
                    fields = ["Density", "x-velocity", "y-velocity", "z-velocity", "Pressure",
                              "Metallicity", "HII", "HeII", "HeIII"]
                else:
                    mylog.info('Detected RAMSES-RT file WITH IR trapping.')
                    fields = ["Density", "x-velocity", "y-velocity", "z-velocity", "Pres_IR",
                              "Pressure", "Metallicity", "HII", "HeII", "HeIII"]
            else:
                if nvar < 5:
                    mylog.debug("nvar=%s is too small! YT doesn't currently support 1D/2D runs in RAMSES %s")
                    raise ValueError
                # Basic hydro runs
                if nvar == 5:
                    fields = ["Density",
                              "x-velocity", "y-velocity", "z-velocity",
                              "Pressure"]
                if nvar > 5 and nvar < 11:
                    fields = ["Density",
                              "x-velocity", "y-velocity", "z-velocity",
                              "Pressure", "Metallicity"]
                # MHD runs - NOTE: THE MHD MODULE WILL SILENTLY ADD 3 TO THE NVAR IN THE MAKEFILE
                if nvar == 11:
                    fields = ["Density",
                              "x-velocity", "y-velocity", "z-velocity",
                              "x-Bfield-left", "y-Bfield-left", "z-Bfield-left",
                              "x-Bfield-right", "y-Bfield-right", "z-Bfield-right",
                              "Pressure"]
                if nvar > 11:
                    fields = ["Density",
                              "x-velocity", "y-velocity", "z-velocity",
                              "x-Bfield-left", "y-Bfield-left", "z-Bfield-left",
                              "x-Bfield-right", "y-Bfield-right", "z-Bfield-right",
                              "Pressure","Metallicity"]
            mylog.debug("No fields specified by user; automatically setting fields array to %s"
                        % str(fields))

        # Allow some wiggle room for users to add too many variables
        count_extra = 0
        while len(fields) < nvar:
            fields.append("var"+str(len(fields)))
            count_extra += 1
        if count_extra > 0:
            mylog.debug('Detected %s extra fluid fields.' % count_extra)
        cls.field_list = [(cls.ftype, f) for f in fields]

        return fields

    @property
    def offset(self):
        if getattr(self, '_offset', None) is not None:
            return self._offset

        with open(self.fname, 'rb') as f:
            # Skip header
            fpu.skip(f, 6)

            # It goes: level, CPU, 8-variable (1 cube)
            min_level = self.domain.ds.min_level
            n_levels = self.domain.amr_header['nlevelmax'] - min_level
            offset = np.zeros(n_levels, dtype='int64')
            offset -= 1
            level_count = np.zeros(n_levels, dtype='int64')
            skipped = []
            amr_header = self.domain.amr_header
            for level in range(amr_header['nlevelmax']):
                for cpu in range(amr_header['nboundary'] +
                                 amr_header['ncpu']):
                    header = ( ('file_ilevel', 1, 'I'),
                               ('file_ncache', 1, 'I') )
                    try:
                        hvals = fpu.read_attrs(f, header, "=")
                    except AssertionError:
                        mylog.error(
                            "You are running with the wrong number of fields. "
                            "If you specified these in the load command, check the array length. "
                            "In this file there are %s hydro fields." % skipped)
                        raise
                    if hvals['file_ncache'] == 0: continue
                    assert(hvals['file_ilevel'] == level+1)
                    if cpu + 1 == self.domain_id and level >= min_level:
                        offset[level - min_level] = f.tell()
                        level_count[level - min_level] = hvals['file_ncache']
                    skipped = fpu.skip(f, 8 * self.nvar)
        self._offset = offset
        self._level_count = level_count
        return self._offset

    @property
    def level_count(self):
        if getattr(self, '_level_count', None) is not None:
            return self._level_count
        self.offset

        return self._level_count
