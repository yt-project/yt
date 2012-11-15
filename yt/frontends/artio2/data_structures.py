#from _artio_reader import artio_is_valid
from _artio_reader import artio_is_valid, artio_fileset

import os
import stat
import numpy as np

from yt.data_objects.static_output import \
    StaticOutput

from .fields import Artio2FieldInfo, add_artio2_field, KnownArtio2Fields

class Artio2StaticOutput(StaticOutput):
    _handle = None
    _fieldinfo_fallback = Artio2FieldInfo
    _fieldinfo_known = KnownArtio2Fields

    def __init__(self, filename, data_style='artio2'):
        if self._handle is not None : return
        self._filename = filename
        self._fileset_prefix = filename[:-4]
        self._handle = artio_fileset(self._fileset_prefix)
        print self._handle
        StaticOutput.__init__(self, filename, data_style)

        self.dimensionality = 3
        self.refine_by = 2
        self.parameters["HydroMethod"] = 'art'
        
    def _parse_parameter_file(self) :
        self.unique_identifier = \
            int(os.stat(self._filename)[stat.ST_CTIME])       

        self.domain_dimensions = np.ones(3,dtype='int64') * \
            self._handle.parameters["num_root_cells"][0]**(1./3.)

        self.domain_left_edge = np.zeros(3, dtype="float64")
        self.domain_right_edge = self.domain_dimensions

        self.min_level = 0
        self.max_level = self._handle.parameters["max_refinement_level"][0]

        self.current_time = self._handle.parameters["tl"][0]
  
        # detect cosmology
        if self._handle.parameters["abox"] :
            self.cosmological_simulation = True
            self.omega_lambda = self._handle.parameters["OmegaL"][0]
            self.omega_matter = self._handle.parameters["OmegaM"][0]
            self.hubble_constant = self._handle.parameters["hubble"][0]
            self.current_redshift = 1.0/self._handle.parameters["abox"][0] - 1.

            self.parameters["initial_redshift"] = \
                self._handle.parameters["auni_init"][0]

    def _set_units(self) :
        self.units = {}
        self.time_units = {}
            

    @classmethod
    def _is_valid(self, *args, **kwargs) :
        # a valid artio header file starts with a prefix and ends with .art
        if not args[0].endswith(".art"): return False
        return artio_is_valid(args[0][:-4])
