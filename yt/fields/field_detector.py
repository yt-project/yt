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
from yt.units.yt_array import YTArray
from .field_exceptions import \
    NeedsGridType

class FieldDetector(defaultdict):
    Level = 1
    NumberOfParticles = 1
    _read_exception = None
    _id_offset = 0
    domain_id = 0

    def __init__(self, nd = 16, ds = None, flat = False, field_parameters=None):
        self.nd = nd
        self.flat = flat
        self._spatial = not flat
        self.ActiveDimensions = [nd,nd,nd]
        self.shape = tuple(self.ActiveDimensions)
        self.size = np.prod(self.ActiveDimensions)
        self.LeftEdge = [0.0, 0.0, 0.0]
        self.RightEdge = [1.0, 1.0, 1.0]
        self.dds = np.ones(3, "float64")
        if field_parameters is None:
            self.field_parameters = {}
        else:
            self.field_parameters = field_parameters
        class fake_dataset(defaultdict):
            pass

        if ds is None:
            # required attrs
            ds = fake_dataset(lambda: 1)
            ds["Massarr"] = np.ones(6)
            ds.current_redshift = ds.omega_lambda = ds.omega_matter = \
                ds.cosmological_simulation = 0.0
            ds.gamma = 5./3.0
            ds.hubble_constant = 0.7
            ds.domain_left_edge = np.zeros(3, 'float64')
            ds.domain_right_edge = np.ones(3, 'float64')
            ds.dimensionality = 3
            ds.periodicity = (True, True, True)
        self.ds = ds

        class fake_index(object):
            class fake_io(object):
                def _read_data_set(io_self, data, field):
                    return self._read_data(field)
                _read_exception = RuntimeError
            io = fake_io()
            def get_smallest_dx(self):
                return 1.0

        self.index = fake_index()
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
        if not isinstance(item, tuple):
            field = ("unknown", item)
        else:
            field = item
        finfo = self.ds._get_field_info(*field)
        params, permute_params = finfo._get_needed_parameters(self)
        self.field_parameters.update(params)
        # For those cases where we are guessing the field type, we will
        # need to re-update -- otherwise, our item will always not have the
        # field type.  This can lead to, for instance, "unknown" particle
        # types not getting correctly identified.
        # Note that the *only* way this works is if we also fix our field
        # dependencies during checking.  Bug #627 talks about this.
        item = self.ds._last_freq
        if finfo is not None and finfo._function.__name__ != 'NullFunc':
            try:
                for param, param_v in permute_params.items():
                    for v in param_v:
                        self.field_parameters[param] = v
                        vv = finfo(self)
                if not permute_params:
                    vv = finfo(self)
            except NeedsGridType as exc:
                ngz = exc.ghost_zones
                nfd = FieldDetector(
                    self.nd + ngz * 2, ds = self.ds,
                    field_parameters=self.field_parameters.copy())
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
        elif finfo is not None and finfo.sampling_type == "particle":
            if "particle_position" in (item, item[1]) or \
               "particle_velocity" in (item, item[1]) or \
               "particle_magnetic_field" in (item, item[1]) or \
               "Velocity" in (item, item[1]) or \
               "Velocities" in (item, item[1]) or \
               "Coordinates" in (item, item[1]) or \
               "MagneticField" in (item, item[1]):
                # A vector
                self[item] = \
                  YTArray(np.ones((self.NumberOfParticles, 3)),
                          finfo.units, registry=self.ds.unit_registry)
            elif "FourMetalFractions" in (item, item[1]):
                self[item] = \
                  YTArray(np.ones((self.NumberOfParticles, 4)),
                          finfo.units, registry=self.ds.unit_registry)
            else:
                # Not a vector
                self[item] = \
                  YTArray(np.ones(self.NumberOfParticles),
                          finfo.units, registry=self.ds.unit_registry)
            if item == ('STAR', 'BIRTH_TIME'):
                # hack for the artio frontend so we pass valid times to
                # the artio functions for calculating physical times
                # from internal times
                self[item] *= -0.1
            self.requested.append(item)
            return self[item]
        self.requested.append(item)
        if item not in self:
            self[item] = self._read_data(item)
        return self[item]

    def _debug(self):
        # We allow this to pass through.
        return

    def deposit(self, *args, **kwargs):
        from yt.frontends.stream.data_structures import StreamParticlesDataset
        from yt.data_objects.static_output import ParticleDataset
        if kwargs['method'] == 'mesh_id':
            if isinstance(self.ds, (StreamParticlesDataset, ParticleDataset)):
                raise ValueError
        return np.random.random((self.nd, self.nd, self.nd))

    def smooth(self, *args, **kwargs):
        tr = np.random.random((self.nd, self.nd, self.nd))
        if kwargs['method'] == "volume_weighted":
            return [tr]
        return tr

    def particle_operation(self, *args, **kwargs):
        return None

    def _read_data(self, field_name):
        self.requested.append(field_name)
        finfo = self.ds._get_field_info(*field_name)
        if finfo.sampling_type == "particle":
            self.requested.append(field_name)
            return np.ones(self.NumberOfParticles)
        return YTArray(defaultdict.__missing__(self, field_name),
                       input_units=finfo.units,
                       registry=self.ds.unit_registry)

    fp_units = {
        'bulk_velocity' : 'cm/s',
        'bulk_magnetic_field': 'G',
        'center' : 'cm',
        'normal' : '',
        'cp_x_vec': '',
        'cp_y_vec': '',
        'cp_z_vec': '',
        'x_hat': '',
        'y_hat': '',
        'z_hat': '',
        'omega_baryon': '',
        'virial_radius': 'cm',
        'observer_redshift': '',
        'source_redshift': '',
        }

    def get_field_parameter(self, param, default = 0.0):
        if self.field_parameters and param in self.field_parameters:
            return self.field_parameters[param]
        self.requested_parameters.append(param)
        if param in ['center', 'normal'] or param.startswith('bulk'):
            return self.ds.arr(
                np.random.random(3) * 1e-2, self.fp_units[param])
        elif param in ['surface_height']:
            return self.ds.quan(0.0, 'code_length')
        elif param in ['axis']:
            return 0
        elif param.startswith("cp_"):
            ax = param[3]
            rv = self.ds.arr((0.0, 0.0, 0.0), self.fp_units[param])
            rv['xyz'.index(ax)] = 1.0
            return rv
        elif param.endswith("_hat"):
            ax = param[0]
            rv = YTArray((0.0, 0.0, 0.0), self.fp_units[param])
            rv['xyz'.index(ax)] = 1.0
            return rv
        elif param == "fof_groups":
            return None
        elif param == "mu":
            return 1.0
        else:
            return default

    _num_ghost_zones = 0
    id = 1

    def apply_units(self, arr, units):
        return self.ds.arr(arr, input_units = units)

    def has_field_parameter(self, param):
        return param in self.field_parameters

    @property
    def fcoords(self):
        fc = np.array(np.mgrid[0:1:self.nd*1j,
                               0:1:self.nd*1j,
                               0:1:self.nd*1j])
        if self.flat:
            fc.shape = (self.nd*self.nd*self.nd, 3)
        else:
            fc = fc.transpose()
        return self.ds.arr(fc, input_units = "code_length")

    @property
    def fcoords_vertex(self):
        fc = np.random.random((self.nd, self.nd, self.nd, 8, 3))
        if self.flat:
            fc.shape = (self.nd*self.nd*self.nd, 8, 3)
        return self.ds.arr(fc, input_units = "code_length")

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
        fw = np.ones((self.nd**3, 3), dtype="float64") / self.nd
        if not self.flat:
            fw.shape = (self.nd, self.nd, self.nd, 3)
        return self.ds.arr(fw, input_units = "code_length")

