"""
Import the components of the volume rendering extension

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2009 Matthew Turk.  All Rights Reserved.

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
from grid_partitioner import HomogenizedVolume
from TransferFunction import ProjectionTransferFunction
from yt.funcs import *
import yt.amr_utils as au

class Camera(object):
    def __init__(self, center, normal_vector, width,
                 resolution, transfer_function,
                 north_vector = None,
                 volume = None, fields = None,
                 log_fields = None,
                 sub_samples = 5, pf = None):
        if pf is not None: self.pf = pf
        if not iterable(resolution):
            resolution = (resolution, resolution)
        self.resolution = resolution
        self.sub_samples = sub_samples
        if not iterable(width):
            width = (width, width, width) # front/back, left/right, top/bottom
        self.width = width
        self.center = center
        if fields is None: fields = ["Density"]
        self.fields = fields
        if transfer_function is None:
            transfer_function = ProjectionTransferFunction()
        self.transfer_function = transfer_function
        self._setup_normalized_vectors(normal_vector, north_vector)
        self.log_fields = log_fields
        if volume is None:
            volume = HomogenizedVolume(fields, pf = self.pf,
                                       log_fields = log_fields)
        self.volume = volume

    def _setup_normalized_vectors(self, normal_vector, north_vector):
        # Now we set up our various vectors
        normal_vector /= na.sqrt( na.dot(normal_vector, normal_vector))
        if north_vector is None:
            vecs = na.identity(3)
            t = na.cross(normal_vector, vecs).sum(axis=1)
            ax = t.argmax()
            north_vector = na.cross(vecs[ax,:], normal_vector).ravel()
        north_vector /= na.sqrt(na.dot(north_vector, north_vector))
        east_vector = -na.cross(north_vector, normal_vector).ravel()
        east_vector /= na.sqrt(na.dot(east_vector, east_vector))
        self.unit_vectors = [north_vector, east_vector, normal_vector]
        self.box_vectors = na.array([self.unit_vectors[0]*self.width[0],
                                     self.unit_vectors[1]*self.width[1],
                                     self.unit_vectors[2]*self.width[2]])

        self.origin = self.center - 0.5*self.width[0]*self.unit_vectors[0] \
                                  - 0.5*self.width[1]*self.unit_vectors[1] \
                                  - 0.5*self.width[2]*self.unit_vectors[2]
        self.back_center = self.center - 0.5*self.width[0]*self.unit_vectors[2]
        self.front_center = self.center + 0.5*self.width[0]*self.unit_vectors[2]
        self.inv_mat = na.linalg.pinv(self.unit_vectors)

    def look_at(self, new_center, north_vector = None):
        normal_vector = self.front_center - new_center
        self._setup_normalized_vectors(normal_vector, north_vector)

    def get_vector_plane(self, image):
        # We should move away from pre-generation of vectors like this and into
        # the usage of on-the-fly generation in the VolumeIntegrator module
        # We might have a different width and back_center
        px = na.linspace(-self.width[0]/2.0, self.width[0]/2.0,
                         self.resolution[0])[:,None]
        py = na.linspace(-self.width[1]/2.0, self.width[1]/2.0,
                         self.resolution[1])[None,:]
        inv_mat = self.inv_mat
        bc = self.back_center
        vectors = na.zeros((self.resolution[0], self.resolution[1], 3),
                          dtype='float64', order='C')
        vectors[:,:,0] = inv_mat[0,0]*px+inv_mat[0,1]*py+self.back_center[0]
        vectors[:,:,1] = inv_mat[1,0]*px+inv_mat[1,1]*py+self.back_center[1]
        vectors[:,:,2] = inv_mat[2,0]*px+inv_mat[2,1]*py+self.back_center[2]
        bounds = (px.min(), px.max(), py.min(), py.max())
        vector_plane = au.VectorPlane(vectors, self.box_vectors[2],
                                      self.back_center, bounds, image,
                                      self.unit_vectors[0],
                                      self.unit_vectors[1])
        return vector_plane

    def snapshot(self):
        image = na.zeros((self.resolution[0], self.resolution[1], 3),
                         dtype='float64', order='C')
        vector_plane = self.get_vector_plane(image)
        tfp = au.TransferFunctionProxy(self.transfer_function) # Reset it every time
        tfp.ns = self.sub_samples
        self.volume.initialize_source()
        pbar = get_pbar("Ray casting",
                        (self.volume.brick_dimensions + 1).prod(axis=-1).sum())
        total_cells = 0
        for brick in self.volume.traverse(self.back_center, self.front_center):
            brick.cast_plane(tfp, vector_plane)
            total_cells += na.prod(brick.my_data[0].shape)
            pbar.update(total_cells)
        pbar.finish()
        return image

    def zoom(self, factor):
        self.width = [w / factor for w in self.width]
        self._setup_normalized_vectors(
                self.unit_vectors[2], self.unit_vectors[0])

    def zoomin(self, final, n_steps):
        f = final**(1.0/n_steps)
        for i in xrange(n_steps):
            self.zoom(f)
            yield self.snapshot()

class StereoPairCamera(Camera):
    def __init__(self, original_camera, relative_separation):
        self.original_camera = original_camera
        self.relative_separation = relative_separation

    def split(self):
        oc = self.original_camera
        uv = oc.unit_vectors
        c = oc.center
        wx = oc.width[0]
        wy = oc.width[1]
        left_center = c + uv[1] * 0.5*self.relative_separation * wx 
        right_center = c - uv[1] * 0.5*self.relative_separation * wx
        left_camera = Camera(left_center, uv[2], oc.width,
                             oc.resolution, oc.transfer_function, uv[0],
                             oc.volume, oc.fields, oc.log_fields,
                             oc.sub_samples, oc.pf)
        right_camera = Camera(right_center, uv[2], oc.width,
                             oc.resolution, oc.transfer_function, uv[0],
                             oc.volume, oc.fields, oc.log_fields,
                             oc.sub_samples, oc.pf)
        return (left_camera, right_camera)
