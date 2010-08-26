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
from transfer_functions import ProjectionTransferFunction
from yt.funcs import *
import yt.amr_utils as au

class Camera(object):
    def __init__(self, center, normal_vector, width,
                 resolution, transfer_function,
                 north_vector = None,
                 volume = None, fields = None,
                 log_fields = None,
                 sub_samples = 5, pf = None):
        r"""A viewpoint into a volume, for volume rendering.

        The camera represents the eye of an observer, which will be used to
        generate ray-cast volume renderings of the domain.

        Parameters
        ----------
        center : array_like
            The current "center" of the view port -- the focal point for the
            camera.
        normal_vector : array_like
            The vector between the camera position and the center.
        width : float or list of floats
            The current width of the image.  If a single float, the volume is
            cubical, but if not, it is front/back, left/right, top/bottom.
        resolution : int or list of ints
            The number of pixels in each direction.
        north_vector : array_like, optional
            The "up" direction for the plane of rays.  If not specific, calculated
            automatically.
        volume : `yt.extensions.volume_rendering.HomogenizedVolume`, optional
            The volume to ray cast through.  Can be specified for finer-grained
            control, but otherwise will be automatically generated.
        fields : list of fields, optional
            This is the list of fields we want to volume render; defaults to
            Density.
        log_fields : list of bool, optional
            Whether we should take the log of the fields before supplying them to
            the volume rendering mechanism.
        sub_samples : int, optional
            The number of samples to take inside every cell per ray.
        pf : `~yt.lagos.StaticOutput`
            For now, this is a require parameter!  But in the future it will become
            optional.  This is the parameter file to volume render.

        Examples
        --------

        >>> cam = vr.Camera(c, L, W, (N,N), transfer_function = tf, pf = pf)
        >>> image = cam.snapshot()
        """
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
        r"""Change the view direction based on a new focal point.

        This will recalculate all the necessary vectors and vector planes related
        to a camera to point at a new location.

        Parameters
        ----------
        new_center : array_like
            The new "center" of the view port -- the focal point for the
            camera.
        north_vector : array_like, optional
            The "up" direction for the plane of rays.  If not specific,
            calculated automatically.
        """
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
        positions = na.zeros((self.resolution[0], self.resolution[1], 3),
                          dtype='float64', order='C')
        positions[:,:,0] = inv_mat[0,0]*px+inv_mat[0,1]*py+self.back_center[0]
        positions[:,:,1] = inv_mat[1,0]*px+inv_mat[1,1]*py+self.back_center[1]
        positions[:,:,2] = inv_mat[2,0]*px+inv_mat[2,1]*py+self.back_center[2]
        bounds = (px.min(), px.max(), py.min(), py.max())
        vector_plane = au.VectorPlane(positions, self.box_vectors[2],
                                      self.back_center, bounds, image,
                                      self.unit_vectors[0],
                                      self.unit_vectors[1])
        return vector_plane

    def snapshot(self):
        r"""Ray-cast the camera.

        This method instructs the camera to take a snapshot -- i.e., call the ray
        caster -- based on its current settings.

        Returns
        -------
        image : array
            An (N,M,3) array of the final returned values, in float64 form.
        """
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
        r"""Change the distance to the focal point.

        This will zoom the camera in by some `factor` toward the focal point,
        along the current view direction, modifying the left/right and up/down
        extents as well.

        Parameters
        ----------
        factor : float
            The factor by which to reduce the distance to the focal point.


        Notes
        -----

        You will need to call snapshot() again to get a new image.

        """
        self.width = [w / factor for w in self.width]
        self._setup_normalized_vectors(
                self.unit_vectors[2], self.unit_vectors[0])

    def zoomin(self, final, n_steps):
        r"""Loop over a zoomin and return snapshots along the way.

        This will yield `n_steps` snapshots until the current view has been
        zooming in to a final factor of `final`.

        Parameters
        ----------
        final : float
            The zoom factor, with respect to current, desired at the end of the
            sequence.
        n_steps : int
            The number of zoom snapshots to make.


        Examples
        --------

        >>> for i, snapshot in enumerate(cam.zoomin(100.0, 10)):
        ...     iw.write_bitmap(snapshot, "zoom_%04i.png" % i)
        """
        f = final**(1.0/n_steps)
        for i in xrange(n_steps):
            self.zoom(f)
            yield self.snapshot()


class PerspectiveCamera(Camera):
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
        positions = na.zeros((self.resolution[0], self.resolution[1], 3),
                          dtype='float64', order='C')
        positions[:,:,0] = inv_mat[0,0]*px+inv_mat[0,1]*py+self.back_center[0]
        positions[:,:,1] = inv_mat[1,0]*px+inv_mat[1,1]*py+self.back_center[1]
        positions[:,:,2] = inv_mat[2,0]*px+inv_mat[2,1]*py+self.back_center[2]
        bounds = (px.min(), px.max(), py.min(), py.max())
        
        # We are likely adding on an odd cutting condition here
        vectors = self.front_center - positions
        n = vectors * vectors
        n = na.sum(n, axis=2)**0.5
        vectors /= n[:,:,None]
        print vectors.shape, vectors.dtype, vectors.flags

        vector_plane = au.VectorPlane(positions, vectors,
                                      self.back_center, bounds, image,
                                      self.unit_vectors[0],
                                      self.unit_vectors[1])
        return vector_plane

class StereoPairCamera(Camera):
    def __init__(self, original_camera, relative_separation = 0.005):
        self.original_camera = original_camera
        self.relative_separation = relative_separation

    def split(self):
        oc = self.original_camera
        uv = oc.unit_vectors
        c = oc.center
        fc = oc.front_center
        wx, wy, wz = oc.width
        left_normal = fc + uv[1] * 0.5*self.relative_separation * wx - c
        right_normal = fc - uv[1] * 0.5*self.relative_separation * wx - c
        left_camera = Camera(c, left_normal, oc.width,
                             oc.resolution, oc.transfer_function, uv[0],
                             oc.volume, oc.fields, oc.log_fields,
                             oc.sub_samples, oc.pf)
        right_camera = Camera(c, right_normal, oc.width,
                             oc.resolution, oc.transfer_function, uv[0],
                             oc.volume, oc.fields, oc.log_fields,
                             oc.sub_samples, oc.pf)
        return (left_camera, right_camera)
