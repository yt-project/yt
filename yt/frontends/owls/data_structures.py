"""
Data structures for OWLS frontend




"""
from __future__ import print_function

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os

from yt.utilities.on_demand_imports import _h5py as h5py

import yt.units
from yt.frontends.gadget.data_structures import \
    GadgetHDF5Dataset
from yt.utilities.definitions import \
    sec_conversion

from .fields import \
    OWLSFieldInfo

class OWLSDataset(GadgetHDF5Dataset):
    _particle_mass_name = "Mass"
    _field_info_class = OWLSFieldInfo
    _time_readin = "Time_GYR"


    def _parse_parameter_file(self):

        # read values from header
        hvals = self._get_hvals()
        self.parameters = hvals

        # set features common to OWLS and Eagle
        self._set_owls_eagle()

        # Set time from value in header
        self.current_time = hvals[self._time_readin] * \
                            sec_conversion["Gyr"] * yt.units.s


    def _set_code_unit_attributes(self):
        self._set_owls_eagle_units()


    @classmethod
    def _is_valid(self, *args, **kwargs):
        need_groups = ['Constants', 'Header', 'Parameters', 'Units']
        veto_groups = ['SUBFIND', 'FOF',
                       'PartType0/ChemistryAbundances', 
                       'PartType0/ChemicalAbundances',
                       'RuntimePars', 'HashTable']
        valid = True
        valid_fname = args[0]
        # If passed arg is a directory, look for the .0 file in that dir
        if os.path.isdir(args[0]):
            valid_files = []
            for f in os.listdir(args[0]):
                fname = os.path.join(args[0], f)
                if ('.0' in f) and ('.ewah' not in f) and os.path.isfile(fname):
                    valid_files.append(fname)
            if len(valid_files) == 0:
                valid = False
            elif len(valid_files) > 1:
                valid = False
            else:
                valid_fname = valid_files[0]
        try:
            fileh = h5py.File(valid_fname, mode='r')
            for ng in need_groups:
                if ng not in fileh["/"]:
                    valid = False
            for vg in veto_groups:
                if vg in fileh["/"]:
                    valid = False                    
            fileh.close()
        except:
            valid = False
            pass
        return valid
