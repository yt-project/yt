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

from yt.funcs import *

from .grid_partitioner import HomogenizedVolume
from .transfer_functions import ProjectionTransferFunction

from yt.utilities.amr_utils import TransferFunctionProxy, VectorPlane
from yt.visualization.image_writer import write_bitmap
from yt.data_objects.data_containers import data_object_registry
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface
from yt.utilities.amr_kdtree.api import AMRKDTree

class Camera(ParallelAnalysisInterface):
    def __init__(self, center, normal_vector, width,
                 resolution, transfer_function,
                 north_vector = None, steady_north=False,
                 volume = None, fields = None,
                 log_fields = None,
                 sub_samples = 5, pf = None,
                 use_kd=True, l_max=None, no_ghost=False,
                 tree_type='domain',expand_factor=1.0,
                 le=None, re=None):
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
            The 'up' direction for the plane of rays.  If not specific, calculated
            automatically.
        steady_north : bool, optional
            Boolean to control whether to normalize the north_vector
            by subtracting off the dot product of it and the normal
            vector.  Makes it easier to do rotations along a single
            axis.  If north_vector is specifies, is switched to
            True. Default: False
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
        use_kd: bool, optional
            Specifies whether or not to use a kd-Tree framework for
            the Homogenized Volume and ray-casting.  Default to True.
        l_max: int, optional
            Specifies the maximum level to be rendered.  Also
            specifies the maximum level used in the kd-Tree
            construction.  Defaults to None (all levels), and only
            applies if use_kd=True.
        no_ghost: bool, optional
            Optimization option.  If True, homogenized bricks will
            extrapolate out from grid instead of interpolating from
            ghost zones that have to first be calculated.  This can
            lead to large speed improvements, but at a loss of
            accuracy/smoothness in resulting image.  The effects are
            less notable when the transfer function is smooth and
            broad. Default: False
        tree_type: string, optional
            Specifies the type of kd-Tree to be constructed/cast.
            There are three options, the default being 'domain'. Only
            affects parallel rendering.  'domain' is suggested.

            'domain' - Tree construction/casting is load balanced by
            splitting up the domain into the first N subtrees among N
            processors (N must be a power of 2).  Casting then
            proceeds with each processor rendering their subvolume,
            and final image is composited on the root processor.  The
            kd-Tree is never combined, reducing communication and
            memory overhead. The viewpoint can be changed without
            communication or re-partitioning of the data, making it
            ideal for rotations/spins.

            'breadth' - kd-Tree is first constructed as in 'domain',
            but then combined among all the subtrees.  Rendering is
            then split among N processors (again a power of 2), based
            on the N most expensive branches of the tree.  As in
            'domain', viewpoint can be changed without re-partitioning
            or communication.

            'depth' - kd-Tree is first constructed as in 'domain', but
            then combined among all subtrees.  Rendering is then load
            balanced in a back-to-front manner, splitting up the cost
            as evenly as possible.  If the viewpoint changes,
            additional data might have to be partitioned.  Is also
            prone to longer data IO times.  If all the data can fit in
            memory on each cpu, this can be the fastest option for
            multiple ray casts on the same dataset.
        expand_factor: float, optional
            A parameter to be used with the PerspectiveCamera.
            Controls how much larger a volume to render, which is
            currently difficult to gauge for the PerspectiveCamera.
            For full box renders, values in the 2.0-3.0 range seem to
            produce desirable results. Default: 1.0
        le: array_like, optional
            Specifies the left edge of the volume to be rendered.
            Currently only works with use_kd=True.
        re: array_like, optional
            Specifies the right edge of the volume to be rendered.
            Currently only works with use_kd=True.

        Examples
        --------

        >>> cam = vr.Camera(c, L, W, (N,N), transfer_function = tf, pf = pf)
        >>> image = cam.snapshot()

        >>> from yt.mods import *
        >>> import yt.visualization.volume_rendering.api as vr
        
        >>> pf = EnzoStaticOutput('DD1701') # Load pf
        >>> c = [0.5]*3 # Center
        >>> L = [1.0,1.0,1.0] # Viewpoint
        >>> W = na.sqrt(3) # Width
        >>> N = 1024 # Pixels (1024^2)

        # Get density min, max
        >>> mi, ma = pf.h.all_data().quantities['Extrema']('Density')[0]
        >>> mi, ma = na.log10(mi), na.log10(ma)

        # Construct transfer function
        >>> tf = vr.ColorTransferFunction((mi-2, ma+2))
        # Sample transfer function with 5 gaussians.  Use new col_bounds keyword.
        >>> tf.add_layers(5,w=0.05, col_bounds = (mi+1,ma), colormap='spectral')
        
        # Create the camera object
        >>> cam = vr.Camera(c, L, W, (N,N), transfer_function=tf, pf=pf) 
        
        # Ray cast, and save the image.
        >>> image = cam.snapshot(fn='my_rendering.png')

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
        self.steady_north = steady_north
        self.expand_factor = expand_factor
        # This seems to be necessary for now.  Not sure what goes wrong when not true.
        if north_vector is not None: self.steady_north=True
        self.north_vector = north_vector
        self.rotation_vector = north_vector
        if fields is None: fields = ["Density"]
        self.fields = fields
        if transfer_function is None:
            transfer_function = ProjectionTransferFunction()
        self.transfer_function = transfer_function
        self._setup_normalized_vectors(normal_vector, north_vector)
        self.log_fields = log_fields
        self.use_kd = use_kd
        self.l_max = l_max
        self.no_ghost = no_ghost
        self.tree_type = tree_type
        if volume is None:
            if self.use_kd:
                volume = AMRKDTree(self.pf, l_max=l_max, fields=self.fields, no_ghost=no_ghost, tree_type=tree_type,
                                   log_fields = log_fields, le=le, re=re)
            else:
                volume = HomogenizedVolume(fields, pf = self.pf,
                                           log_fields = log_fields)
        else:
            self.use_kd = isinstance(volume, AMRKDTree)
        self.volume = volume

    def _setup_normalized_vectors(self, normal_vector, north_vector):
        # Now we set up our various vectors
        normal_vector /= na.sqrt( na.dot(normal_vector, normal_vector))
        if north_vector is None:
            vecs = na.identity(3)
            t = na.cross(normal_vector, vecs).sum(axis=1)
            ax = t.argmax()
            north_vector = na.cross(vecs[ax,:], normal_vector).ravel()
            if self.rotation_vector is None:
                self.rotation_vector=north_vector
        else:
            if self.steady_north:
                north_vector = north_vector - na.dot(north_vector,normal_vector)*normal_vector
        north_vector /= na.sqrt(na.dot(north_vector, north_vector))
        east_vector = -na.cross(north_vector, normal_vector).ravel()
        east_vector /= na.sqrt(na.dot(east_vector, east_vector))
        self.normal_vector = normal_vector
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

    def switch_view(self, normal_vector=None, width=None, center=None, north_vector=None):
        r"""Change the view direction based on any of the view parameters.

        This will recalculate all the necessary vectors and vector planes related
        to a camera with new normal vectors, widths, centers, or north vectors.

        Parameters (All Optional)
        ----------
        normal_vector: array_like, optional
            The new looking vector.
        width: float or array of floats, optional
            The new width.  Can be a single value W -> [W,W,W] or an
            array [W1, W2, W3]
        center: array_like, optional
            Specifies the new center.
        north_vector : array_like, optional
            The 'up' direction for the plane of rays.  If not specific,
            calculated automatically.
        """
        if width is None:
            width = self.width
        if not iterable(width):
            width = (width, width, width) # front/back, left/right, top/bottom
        self.width = width
        if center is not None:
            self.center = center
        if north_vector is None:
            north_vector = self.north_vector
        if normal_vector is None:
            normal_vector = self.front_center-self.center
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
        vector_plane = VectorPlane(positions, self.box_vectors[2],
                                      self.back_center, bounds, image,
                                      self.unit_vectors[0],
                                      self.unit_vectors[1])
        return vector_plane

    def snapshot(self, fn = None):
        r"""Ray-cast the camera.

        This method instructs the camera to take a snapshot -- i.e., call the ray
        caster -- based on its current settings.

        Parameters
        ----------
        fn : string, optional
            If supplied, the image will be saved out to this before being
            returned.  Scaling will be to the maximum value.

        Returns
        -------
        image : array
            An (N,M,3) array of the final returned values, in float64 form.
        """
        image = na.zeros((self.resolution[0], self.resolution[1], 3),
                         dtype='float64', order='C')
        vector_plane = self.get_vector_plane(image)
        tfp = TransferFunctionProxy(self.transfer_function) # Reset it every time
        tfp.ns = self.sub_samples
        self.volume.initialize_source()
        if self.use_kd:
            self.volume.reset_cast()
            image = self.volume.kd_ray_cast(image, tfp, vector_plane,
                                            self.back_center, self.front_center)
        else:
            pbar = get_pbar("Ray casting",
                            (self.volume.brick_dimensions + 1).prod(axis=-1).sum())
            total_cells = 0
            for brick in self.volume.traverse(self.back_center, self.front_center):
                brick.cast_plane(tfp, vector_plane)
                total_cells += na.prod(brick.my_data[0].shape)
                pbar.update(total_cells)
            pbar.finish()

        if self._mpi_get_rank() is 0 and fn is not None:
            write_bitmap(image, fn)
            
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

    def move_to(self, final, n_steps, final_width=None):
        r"""Loop over a look_at

        This will yield `n_steps` snapshots until the current view has been
        moved to a final center of `final` with a final width of final_width.

        Parameters
        ----------
        final : array_like
            The final center to move to after `n_steps`
        n_steps : int
            The number of look_at snapshots to make.
        final_width: float or array_like, optional
            Specifies the final width after `n_steps`.  Useful for
            moving and zooming at the same time.
            
        Examples
        --------

        >>> for i, snapshot in enumerate(cam.move_to([0.2,0.3,0.6], 10)):
        ...     iw.write_bitmap(snapshot, "move_%04i.png" % i)
        """
        self.center = na.array(self.center)
        dW = None
        if final_width is not None:
            if not iterable(final_width):
                width = na.array([final_width, final_width, final_width]) # front/back, left/right, top/bottom
            dW = (1.0*final_width-na.array(self.width))/n_steps
        dx = (na.array(final)-self.center)*1.0/n_steps
        for i in xrange(n_steps):
            self.switch_view(center=self.center+dx, width=self.width+dW)
            yield self.snapshot()

    def rotate(self, theta, rot_vector=None):
        r"""Rotate by a given angle

        Rotate the view.  If `rot_vector` is None, rotation will occur
        around the `north_vector`.

        Parameters
        ----------
        theta : float, in radians
             Angle (in radians) by which to rotate the view.
        rot_vector  : array_like, optional
            Specify the rotation vector around which rotation will
            occur.  Defaults to None, which sets rotation around
            `north_vector`

        Examples
        --------

        >>> cam.rotate(na.pi/4)
        """
        if rot_vector is None:
            rot_vector = self.rotation_vector
            
        ux = rot_vector[0]
        uy = rot_vector[1]
        uz = rot_vector[2]
        cost = na.cos(theta)
        sint = na.sin(theta)
        
        R = na.array([[cost+ux**2*(1-cost), ux*uy*(1-cost)-uz*sint, ux*uz*(1-cost)+uy*sint],
                      [uy*ux*(1-cost)+uz*sint, cost+uy**2*(1-cost), uy*uz*(1-cost)-ux*sint],
                      [uz*ux*(1-cost)-uy*sint, uz*uy*(1-cost)+ux*sint, cost+uz**2*(1-cost)]])

        normal_vector = self.front_center-self.center

        self.switch_view(normal_vector=na.dot(R,normal_vector))


    def rotation(self, theta, n_steps, rot_vector=None):
        r"""Loop over rotate, creating a rotation

        This will yield `n_steps` snapshots until the current view has been
        rotated by an angle `theta`

        Parameters
        ----------
        theta : float, in radians
            Angle (in radians) by which to rotate the view.
        n_steps : int
            The number of look_at snapshots to make.
        rot_vector  : array_like, optional
            Specify the rotation vector around which rotation will
            occur.  Defaults to None, which sets rotation around the
            original `north_vector`

        Examples
        --------

        >>> for i, snapshot in enumerate(cam.rotation(na.pi, 10)):
        ...     iw.write_bitmap(snapshot, "rotation_%04i.png" % i)
        """

        dtheta = (1.0*theta)/n_steps
        for i in xrange(n_steps):
            self.rotate(dtheta, rot_vector=rot_vector)
            yield self.snapshot()

data_object_registry["camera"] = Camera

class PerspectiveCamera(Camera):
    def get_vector_plane(self, image):
        # We should move away from pre-generation of vectors like this and into
        # the usage of on-the-fly generation in the VolumeIntegrator module
        # We might have a different width and back_center
        dl = (self.back_center - self.front_center)
        self.front_center += dl
        self.back_center -= dl
        px = self.expand_factor*na.linspace(-self.width[0]/2.0, self.width[0]/2.0,
                         self.resolution[0])[:,None]
        py = self.expand_factor*na.linspace(-self.width[1]/2.0, self.width[1]/2.0,
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
        positions = self.front_center - 2.0*(((self.back_center-self.front_center)**2).sum())**0.5*vectors
        vectors = (self.front_center - positions)
        print vectors.shape, vectors.dtype, vectors.flags

        vector_plane = VectorPlane(positions, vectors,
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
