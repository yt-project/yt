"""
The field detector.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from collections import defaultdict
from yt.utilities.units import Unit
from yt.data_objects.yt_array import YTArray
from .field_exceptions import \
    ValidationException, \
    NeedsGridType, \
    NeedsOriginalGrid, \
    NeedsDataField, \
    NeedsProperty, \
    NeedsParameter, \
    FieldUnitsError

class FieldDetector(defaultdict):
    Level = 1
    NumberOfParticles = 1
    _read_exception = None
    _id_offset = 0

    def __init__(self, nd = 16, pf = None, flat = False):
        self.nd = nd
        self.flat = flat
        self._spatial = not flat
        self.ActiveDimensions = [nd,nd,nd]
        self.shape = tuple(self.ActiveDimensions)
        self.size = np.prod(self.ActiveDimensions)
        self.LeftEdge = [0.0, 0.0, 0.0]
        self.RightEdge = [1.0, 1.0, 1.0]
        self.dds = np.ones(3, "float64")
        class fake_parameter_file(defaultdict):
            pass

        if pf is None:
            # required attrs
            pf = fake_parameter_file(lambda: 1)
            pf["Massarr"] = np.ones(6)
            pf.current_redshift = pf.omega_lambda = pf.omega_matter = \
                pf.cosmological_simulation = 0.0
            pf.gamma = 5./3.0
            pf.hubble_constant = 0.7
            pf.domain_left_edge = np.zeros(3, 'float64')
            pf.domain_right_edge = np.ones(3, 'float64')
            pf.dimensionality = 3
            pf.periodicity = (True, True, True)
        self.pf = pf

        class fake_hierarchy(object):
            class fake_io(object):
                def _read_data_set(io_self, data, field):
                    return self._read_data(field)
                _read_exception = RuntimeError
            io = fake_io()
            def get_smallest_dx(self):
                return 1.0

        self.hierarchy = fake_hierarchy()
        self.requested = []
        self.requested_parameters = []
        if not self.flat:
            defaultdict.__init__(self,
                lambda: np.ones((nd, nd, nd), dtype='float64')
                + 1e-4*np.random.random((nd, nd, nd)))
        else:
            defaultdict.__init__(self,
                lambda: np.ones((nd * nd * nd), dtype='float64')
                + 1e-4*np.random.random((nd * nd * nd)))

    def _reshape_vals(self, arr):
        if not self._spatial: return arr
        if len(arr.shape) == 3: return arr
        return arr.reshape(self.ActiveDimensions, order="C")

    def __missing__(self, item):
        if hasattr(self.pf, "field_info"):
            if not isinstance(item, tuple):
                field = ("unknown", item)
                finfo = self.pf._get_field_info(*field)
                #mylog.debug("Guessing field %s is %s", item, finfo.name)
            else:
                field = item
            finfo = self.pf._get_field_info(*field)
            # For those cases where we are guessing the field type, we will
            # need to re-update -- otherwise, our item will always not have the
            # field type.  This can lead to, for instance, "unknown" particle
            # types not getting correctly identified.
            # Note that the *only* way this works is if we also fix our field
            # dependencies during checking.  Bug #627 talks about this.
            item = self.pf._last_freq
        else:
            FI = getattr(self.pf, "field_info", FieldInfo)
            if item in FI:
                finfo = FI[item]
            else:
                finfo = None
        if finfo is not None and finfo._function.func_name != 'NullFunc':
            try:
                vv = finfo(self)
            except NeedsGridType as exc:
                ngz = exc.ghost_zones
                nfd = FieldDetector(self.nd + ngz * 2, pf = self.pf)
                nfd._num_ghost_zones = ngz
                vv = finfo(nfd)
                if ngz > 0: vv = vv[ngz:-ngz, ngz:-ngz, ngz:-ngz]
                for i in nfd.requested:
                    if i not in self.requested: self.requested.append(i)
                for i in nfd.requested_parameters:
                    if i not in self.requested_parameters:
                        self.requested_parameters.append(i)
            if vv is not None:
                if not self.flat: self[item] = vv
                else: self[item] = vv.ravel()
                return self[item]
        elif finfo is not None and finfo.particle_type:
            if item == "Coordinates" or item[1] == "Coordinates" or \
               item == "Velocities" or item[1] == "Velocities" or \
               item == "Velocity" or item[1] == "Velocity":
                # A vector
                self[item] = \
                  YTArray(np.ones((self.NumberOfParticles, 3)),
                          finfo.units, registry=self.pf.unit_registry)
            else:
                # Not a vector
                self[item] = \
                  YTArray(np.ones(self.NumberOfParticles),
                          finfo.units, registry=self.pf.unit_registry)
            self.requested.append(item)
            return self[item]
        self.requested.append(item)
        if item not in self:
            self[item] = self._read_data(item)
        return self[item]

    def deposit(self, *args, **kwargs):
        return np.random.random((self.nd, self.nd, self.nd))

    def _read_data(self, field_name):
        self.requested.append(field_name)
        if hasattr(self.pf, "field_info"):
            finfo = self.pf._get_field_info(*field_name)
        else:
            finfo = FieldInfo[field_name]
        if finfo.particle_type:
            self.requested.append(field_name)
            return np.ones(self.NumberOfParticles)
        return YTArray(defaultdict.__missing__(self, field_name),
                       input_units=finfo.units,
                       registry=self.pf.unit_registry)

    fp_units = {
        'bulk_velocity' : 'cm/s',
        'center' : 'cm',
        'normal' : '',
        'cp_x_vec': '',
        'cp_y_vec': '',
        'cp_z_vec': '',
        }

    def get_field_parameter(self, param, default = None):
        self.requested_parameters.append(param)
        if param in ['bulk_velocity', 'center', 'normal']:
            return YTArray(np.random.random(3) * 1e-2, self.fp_units[param])
        elif param in ['axis']:
            return 0
        elif param.startswith("cp_"):
            ax = param[3]
            rv = YTArray((0.0, 0.0, 0.0), self.fp_units[param])
            rv['xyz'.index(ax)] = 1.0
            return rv
        else:
            return 0.0

    _num_ghost_zones = 0
    id = 1

    def has_field_parameter(self, param):
        return True

    @property
    def fcoords(self):
        fc = np.array(np.mgrid[0:1:self.nd*1j,
                               0:1:self.nd*1j,
                               0:1:self.nd*1j])
        if self.flat:
            fc.shape = (self.nd*self.nd*self.nd, 3)
        else:
            fc = fc.transpose()
        return fc

    @property
    def icoords(self):
        ic = np.mgrid[0:self.nd-1:self.nd*1j,
                      0:self.nd-1:self.nd*1j,
                      0:self.nd-1:self.nd*1j]
        if self.flat:
            ic.shape = (self.nd*self.nd*self.nd, 3)
        else:
            ic = ic.transpose()
        return ic

    @property
    def ires(self):
        ir = np.ones(self.nd**3, dtype="int64")
        if not self.flat:
            ir.shape = (self.nd, self.nd, self.nd)
        return ir

    @property
    def fwidth(self):
        fw = np.ones(self.nd**3, dtype="float64") / self.nd
        if not self.flat:
            fw.shape = (self.nd, self.nd, self.nd)
        return fw

