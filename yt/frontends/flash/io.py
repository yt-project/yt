"""
FLASH-specific IO functions

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD
Homepage: http://yt-project.org/
License:
  Copyright (C) 2010-2011 Matthew Turk.  All Rights Reserved.

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

import numpy as na
import h5py

from yt.utilities.io_handler import \
    BaseIOHandler

def particles_validator_region(x, y, z, args) :

    left_edge = args[0]
    right_edge = args[1]
    periodic = args[2]
    DLE = args[3]
    DRE = args[4]

    xx = x
    yy = y
    zz = z

    if periodic == 1 : 

        DW = DRE - DLE
        xx[x < left_edge[0]] = x + DW[0]
        xx[x > right_edge[0]] = x - DW[0]
        yy[y < left_edge[1]] = y + DW[1]
        yy[y > right_edge[1]] = y - DW[1]
        zz[z < left_edge[2]] = z + DW[2]
        zz[z > right_edge[2]] = z - DW[2]

    idxx = na.logical_and(xx >= left_edge[0], xx <= right_edge[0])
    idxy = na.logical_and(yy >= left_edge[1], yy <= right_edge[1])
    idxz = na.logical_and(zz >= left_edge[2], zz <= right_edge[2])

    idxs = na.logical_and(idxx, idxy)
    idxs = na.logical_and(idxz, idxs)

    return idxs

def particles_validator_sphere(x, y, z, args) :
    
    center = args[0]
    radius = args[1]
    periodic = args[2]
    DLE = args[3]
    DRE = args[4]

    xx = na.abs(x-center[0])
    yy = na.abs(y-center[1])
    zz = na.abs(z-center[2])

    if periodic == 1 : 

        DW = DRE - DLE

        xx = na.minimum(xx,DW[0]-xx)
        yy = na.minimum(yy,DW[1]-yy)
        zz = na.minimum(zz,DW[2]-zz)

    r = na.sqrt(xx*xx+yy*yy+zz*zz)

    return r <= radius

def particles_validator_disk(x, y, z, args) :
    
    center = args[0]
    normal = args[1]
    radius = args[2]
    height = args[3]

    d = -na.dot(normal*center)

    ph = na.abs(x*normal[0] + y*normal[1] + z*normal[2] + d)
    pd2 = (x-center[0])**2+(y-center[1])**2+(z-center[2])**2

    pr = na.sqrt(pd2-ph*ph)

    return na.logical_and(pr <= radius, ph <= height)

class IOHandlerFLASH(BaseIOHandler):
    _particle_reader = True
    _data_style = "flash_hdf5"

    def __init__(self, pf, *args, **kwargs):
        self._num_per_stride = kwargs.pop("num_per_stride", 1000000)
        BaseIOHandler.__init__(self, *args, **kwargs)
        # Now we cache the particle fields
        self.pf = pf
        self._handle = pf._handle
        try :
            particle_fields = [s[0].strip() for s in
                               self._handle["/particle names"][:]]
            self._particle_fields = dict([("particle_" + s, i) for i, s in
                                          enumerate(particle_fields)])
        except KeyError:
            self._particle_fields = {}

    def _read_particles(self, fields_to_read, type, args, grid_list,
            count_list, conv_factors):
        f = self._handle
        particles = []
        _particles = f["/tracer particles"][:,:]
        fx = self._particle_fields["particle_posx"]
        fy = self._particle_fields["particle_posy"]
        fz = self._particle_fields["particle_posz"]
        posx = _particles[:,fx]
        posy = _particles[:,fy]
        posz = _particles[:,fz]
        if type == 0 :
            idxs = particles_validator_region(posx,posy,posz,args)
        elif type == 1 :
            idxs = particles_validator_sphere(posx,posy,posz,args)
        elif type == 2 :
            idxs = particles_validator_disk(posx,posy,posz,args)
        for field in fields_to_read :
            fi = self._particle_fields[field]
            particles.append(_particles[idxs,fi])
        del _particles
        return particles

    def _read_data_set(self, grid, field):
        f = self._handle
        if field in self._particle_fields:
            if grid.NumberOfParticles == 0: return na.array([], dtype='float64')
            start = self.pf.h._particle_indices[grid.id - grid._id_offset]
            end = self.pf.h._particle_indices[grid.id - grid._id_offset + 1]
            fi = self._particle_fields[field]
            tr = f["/tracer particles"][start:end, fi]
        else:
            tr = f["/%s" % field][grid.id - grid._id_offset,:,:,:].transpose()
        return tr.astype("float64")

    def _read_data_slice(self, grid, field, axis, coord):
        sl = [slice(None), slice(None), slice(None)]
        sl[axis] = slice(coord, coord + 1)
        f = self._handle
        tr = f["/%s" % field][grid.id - grid._id_offset].transpose()[sl]
        return tr.astype("float64")

