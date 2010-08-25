"""
Generalized Enzo output objects, both static and time-series.

Presumably at some point EnzoRun will be absorbed into here.
Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2009 Matthew Turk, J. S. Oishi.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from yt.lagos import *
from yt.fido import ParameterFileStore, NoParameterShelf
from yt.funcs import *
import string, re, gc, time, os, os.path

# We want to support the movie format in the future.
# When such a thing comes to pass, I'll move all the stuff that is contant up
# to here, and then have it instantiate EnzoStaticOutputs as appropriate.

_cached_pfs = weakref.WeakValueDictionary()
_pf_store = ParameterFileStore()

class StaticOutput(object):
    class __metaclass__(type):
        def __init__(cls, name, b, d):
            type.__init__(cls, name, b, d)
            output_type_registry[name]=cls
            mylog.debug("Registering: %s as %s", name, cls)

    def __new__(cls, filename=None, *args, **kwargs):
        if not isinstance(filename, types.StringTypes): 
            obj = object.__new__(cls)
            obj.__init__(filename, *args, **kwargs)
            return obj
        apath = os.path.abspath(filename)
        if not os.path.exists(apath): raise IOError(filename)
        if apath not in _cached_pfs:
            obj = object.__new__(cls)
            obj.__init__(filename, *args, **kwargs)
            _cached_pfs[apath] = obj
            if ytcfg.getboolean('lagos','serialize'):
                try:
                    _pf_store.check_pf(obj)
                except NoParameterShelf:
                    pass
        return _cached_pfs[apath]

    def __init__(self, filename, data_style=None):
        """
        Base class for generating new output types.  Principally consists of
        a *filename* and a *data_style* which will be passed on to children.
        """
        self.data_style = data_style
        self.parameter_filename = str(filename)
        self.basename = os.path.basename(filename)
        self.directory = os.path.expanduser(os.path.dirname(filename))
        self.fullpath = os.path.abspath(self.directory)
        self._instantiated = time.time()
        if len(self.directory) == 0:
            self.directory = "."
        self.conversion_factors = {}
        self.parameters = {}
        self._parse_parameter_file()
        self._set_units()
        # These can be taken out if you so desire

    def __reduce__(self):
        args = (self._hash(),)
        return (_reconstruct_pf, args)

    def __repr__(self):
        return self.basename

    def _hash(self):
        s = "%s;%s;%s" % (self.basename,
            self["InitialTime"], self["CurrentTimeIdentifier"])
        try:
            import hashlib
            return hashlib.md5(s).hexdigest()
        except ImportError:
            return s.replace(";", "*")

    @classmethod
    def _is_valid(cls, *args, **kwargs):
        return False

    def __getitem__(self, key):
        """
        Returns _units, parameters, or _conversion_factors in that order
        """
        for d in [self.units, self.time_units, self.parameters, \
                  self.conversion_factors]:
            if key in d: return d[key]
        raise KeyError(key)

    def keys(self):
        """
        Returns a list of possible keys, from _units, parameters and
        _conversion_factors
        """
        return self.units.keys() \
             + self.time_units.keys() \
             + self.parameters.keys() \
             + self.conversion_factors.keys()

    def __iter__(self):
        for ll in [self.units, self.time_units,
                   self.parameters, self.conversion_factors]:
            for i in ll.keys(): yield i

    def get_smallest_appropriate_unit(self, v):
        max_nu = 1e30
        good_u = None
        for unit in ['mpc','kpc','pc','au','rsun','cm']:
            vv = v*self[unit]
            if vv < max_nu and vv > 1.0:
                good_u = unit
                max_nu = v*self[unit]
        return good_u

    def has_key(self, key):
        """
        Returns true or false
        """
        return key in self.units or \
               key in self.time_units or \
               key in self.parameters or \
               key in self.conversion_factors

    def _get_hierarchy(self):
        if self.__hierarchy == None:
            if self._hierarchy_class == None:
                raise RuntimeError("You should not instantiate StaticOutput.")
            self.__hierarchy = self._hierarchy_class(self, data_style=self.data_style)
        return self.__hierarchy

    def _set_hierarchy(self, newh):
        if self.__hierarchy != None:
            mylog.warning("Overriding hierarchy attribute!  This is probably unwise!")
        self.__hierarchy = newh

    __hierarchy = None
    hierarchy = property(_get_hierarchy, _set_hierarchy)
    h = property(_get_hierarchy, _set_hierarchy)

def _reconstruct_pf(*args, **kwargs):
    pfs = ParameterFileStore()
    pf = pfs.get_pf_hash(*args)
    return pf

class GadgetStaticOutput(StaticOutput):
    _hierarchy_class = GadgetHierarchy
    _fieldinfo_class = GadgetFieldContainer
    def __init__(self, h5filename,storage_filename=None) :
        StaticOutput.__init__(self, h5filename, 'gadget_hdf5')
        self.storage_filename = storage_filename #Don't know what this is
        self.field_info = self._fieldinfo_class()
        x = self._get_param('maxlevel')**2
        self.max_grid_size = (x,x,x)
        self.parameters["InitialTime"] = 0.0
        # These should be explicitly obtained from the file, but for now that
        # will wait until a reorganization of the source tree and better
        # generalization.
        self.parameters["TopGridRank"] = 3
        self.parameters["RefineBy"] = 2
        self.parameters["DomainLeftEdge"] = self.leftedge
        self.parameters["DomainRightEdge"] = self.rightedge
        
        
    def _parse_parameter_file(self):
        # read the units in from the hdf5 file 
        #fill in self.units dict
        #fill in self.time_units dict (keys: 'days','years', '1')
        
        #import all of the parameter file params 
        #this is NOT originally from the gadget snapshot but instead
        #from the paramfile starting the sim
        skips = ('TITLE','CLASS','VERSION') #these are just hdf5 crap
        fh = h5py.File(self.parameter_filename)
        for kw in fh['root'].attrs.keys():
            if any([skip in kw for skip in skips]):
                continue
            val = fh['root'].attrs[kw]
            if type(val)==type(''):
                try:    val = cPickle.loads(val)
                except: pass
            #also, includes unit info
            setattr(self,kw,val)
            
    def _get_param(self,kw,location='/root'):
        fh = h5py.File(self.parameter_filename)
        val = fh[location].attrs[kw]
        try:    val = cPickle.loads(val)
        except: pass
        return val
            
    def _set_units(self):
        #check out the unit params from _parse_parameter_file and use them
        #code below is all filler
        self.units = {}
        self.time_units = {}
        self.conversion_factors = defaultdict(lambda: 1.0)
        self.time_units['1'] = 1
        self.units['1'] = 1.0
        self.units['unitary'] = 1.0
        self.units['cm'] = 1.0
        seconds = 1 #self["Time"]
        self.time_units['years'] = seconds / (365*3600*24.0)
        self.time_units['days']  = seconds / (3600*24.0)
        for key in yt2orionFieldsDict:
            self.conversion_factors[key] = 1.0
        
        
    @classmethod
    def _is_valid(cls, *args, **kwargs):
        # check for a /root to exist in the h5 file
        try:
            h5f=h5py.File(self.h5filename)
            valid = 'root' in h5f.items()[0]
            h5f.close()
            return valid
        except:
            pass
        return False

class TigerStaticOutput(StaticOutput):
    _hierarchy_class = TigerHierarchy
    _fieldinfo_class = TigerFieldContainer

    def __init__(self, rhobname, root_size, max_grid_size=128,
                 data_style='tiger', storage_filename = None):
        StaticOutput.__init__(self, rhobname, data_style)
        self.storage_filename = storage_filename
        self.basename = rhobname[:-4]
        if not os.path.exists(self.basename + "rhob"):
            print "%s doesn't exist, don't know how to handle this!" % (
                        self.basename + "rhob")
            raise IOError
        if not iterable(root_size): root_size = (root_size,) * 3
        self.root_size = root_size
        if not iterable(max_grid_size): max_grid_size = (max_grid_size,) * 3
        self.max_grid_size = max_grid_size

        self.field_info = self._fieldinfo_class()

        # We assume that we have basename + "rhob" and basename + "temp"
        # to get at our various parameters.

        # First we get our our header:
        
        header = [
            ('i', 'dummy0'),
            ('f', 'ZR'),
            ('f', 'OMEGA0'),
            ('f', 'FLAM0'),
            ('f', 'OMEGAB'),
            ('f', 'H0'),
            ('f', 'BOXL0'),
            ('i', 'dummy1'),
            ]

        h_fmt, h_key = zip(*header)
        header_string = "".join(h_fmt)

        fs = open(self.basename + "rhob")
        header_raw = read_struct(fs, header_string)
        self.parameters.update(dict(zip(h_key, header_raw)))

        if "InitialTime" not in self.parameters:
            self.parameters["InitialTime"] = 0.0
        self.parameters["CurrentTimeIdentifier"] = \
            int(os.stat(self.parameter_filename)[ST_CTIME])
        self.parameters['TopGridDimensions'] = root_size
        self.parameters['TopGridRank'] = 3
        self.units["Density"] = 1.0
        self.parameters['RefineBy'] = 2

    def _set_units(self):
        self.parameters["DomainLeftEdge"] = na.zeros(3, dtype='float64')
        self.parameters["DomainRightEdge"] = na.ones(3, dtype='float64')
        self.units = {}
        self.time_units = {}
        self.time_units['1'] = 1
        self.units['1'] = 1.0
        self.units['cm'] = 1.0 # This is just plain false
        self.units['unitary'] = 1.0 / (self["DomainRightEdge"] - self["DomainLeftEdge"]).max()

    def _parse_parameter_file(self):
        pass

    @classmethod
    def _is_valid(self, *args, **kwargs):
        return os.path.exists(args[0] + "rhob")

