"""
Import the components of the volume rendering extension

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
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

import __builtin__
import numpy as na

from yt.funcs import *
from yt.utilities.math_utils import *

from .grid_partitioner import HomogenizedVolume
from .transfer_functions import ProjectionTransferFunction

from yt.utilities.lib import \
    arr_vec2pix_nest, arr_pix2vec_nest, \
    arr_ang2pix_nest, arr_fisheye_vectors
from yt.utilities.math_utils import get_rotation_matrix
from yt.utilities.orientation import Orientation
from yt.visualization.image_writer import write_bitmap, write_image
from yt.data_objects.data_containers import data_object_registry
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface, ProcessorPool
from yt.utilities.amr_kdtree.api import AMRKDTree

from yt.utilities.lib import \
    PartitionedGrid, ProjectionSampler, VolumeRenderSampler, \
    LightSourceRenderSampler, InterpolatedProjectionSampler, \
    arr_vec2pix_nest, arr_pix2vec_nest, arr_ang2pix_nest, \
    pixelize_healpix, arr_fisheye_vectors

class Camera(ParallelAnalysisInterface):

    _sampler_object = VolumeRenderSampler

    def __init__(self, center, normal_vector, width,
                 resolution, transfer_function,
                 north_vector = None, steady_north=False,
                 volume = None, fields = None,
                 log_fields = None,
                 sub_samples = 5, pf = None,
                 use_kd=True, l_max=None, no_ghost=True,
                 tree_type='domain',
                 le=None, re=None, use_light=False):
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
            cubical, but if not, it is left/right, top/bottom, front/back.
        resolution : int or list of ints
            The number of pixels in each direction.
        north_vector : array_like, optional
            The 'up' direction for the plane of rays.  If not specific, calculated
            automatically.
        steady_north : bool, optional
            Boolean to control whether to normalize the north_vector
            by subtracting off the dot product of it and the normal
            vector.  Makes it easier to do rotations along a single
            axis.  If north_vector is specified, is switched to
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
        pf : `~yt.data_objects.api.StaticOutput`
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
        ParallelAnalysisInterface.__init__(self)
        if pf is not None: self.pf = pf
        if not iterable(resolution):
            resolution = (resolution, resolution)
        self.resolution = resolution
        self.sub_samples = sub_samples
        if not iterable(width):
            width = (width, width, width) # left/right, top/bottom, front/back 
        self.orienter = Orientation(normal_vector, north_vector=north_vector, steady_north=steady_north)
        self.rotation_vector = self.orienter.north_vector
        self._setup_box_properties(width, center, self.orienter.unit_vectors)
        if fields is None: fields = ["Density"]
        self.fields = fields
        if transfer_function is None:
            transfer_function = ProjectionTransferFunction()
        self.transfer_function = transfer_function
        self.log_fields = log_fields
        self.use_kd = use_kd
        self.l_max = l_max
        self.no_ghost = no_ghost
        self.use_light = use_light
        self.light_dir = None
        self.light_rgba = None
        if self.no_ghost:
            mylog.info('Warning: no_ghost is currently True (default). This may lead to artifacts at grid boundaries.')
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

    def _setup_box_properties(self, width, center, unit_vectors):
        self.width = width
        self.center = center
        self.box_vectors = na.array([unit_vectors[0]*width[0],
                                     unit_vectors[1]*width[1],
                                     unit_vectors[2]*width[2]])
        self.origin = center - 0.5*na.dot(width,unit_vectors)
        self.back_center =  center - 0.5*width[2]*unit_vectors[2]
        self.front_center = center + 0.5*width[2]*unit_vectors[2]         

    def update_view_from_matrix(self, mat):
        pass

    def look_at(self, new_center, north_vector = None):
        r"""Change the view direction based on a new focal point.

        This will recalculate all the necessary vectors and vector planes to orient
        the image plane so that it points at a new location.

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
        self.orienter.switch_orientation(normal_vector=normal_vector,
                                         north_vector = north_vector)

    def switch_view(self, normal_vector=None, width=None, center=None, north_vector=None):
        r"""Change the view based on any of the view parameters.

        This will recalculate the orientation and width based on any of
        normal_vector, width, center, and north_vector.

        Parameters
        ----------
        normal_vector: array_like, optional
            The new looking vector.
        width: float or array of floats, optional
            The new width.  Can be a single value W -> [W,W,W] or an
            array [W1, W2, W3] (left/right, top/bottom, front/back)
        center: array_like, optional
            Specifies the new center.
        north_vector : array_like, optional
            The 'up' direction for the plane of rays.  If not specific,
            calculated automatically.
        """
        if width is None:
            width = self.width
        if not iterable(width):
            width = (width, width, width) # left/right, tom/bottom, front/back 
        self.width = width
        if center is not None:
            self.center = center
        if north_vector is None:
            north_vector = self.orienter.north_vector
        if normal_vector is None:
            normal_vector = self.orienter.normal_vector
        self.orienter.switch_orientation(normal_vector = normal_vector,
                                         north_vector = north_vector)
        self._setup_box_properties(width, self.center, self.orienter.unit_vectors)
    def new_image(self):
        image = na.zeros((self.resolution[0], self.resolution[1], 3), dtype='float64', order='C')
        return image

    def get_sampler_args(self, image):
        rotp = na.concatenate([self.orienter.inv_mat.ravel('F'), self.back_center.ravel()])
        args = (rotp, self.box_vectors[2], self.back_center,
                (-self.width[0]/2.0, self.width[0]/2.0,
                 -self.width[1]/2.0, self.width[1]/2.0),
                image, self.orienter.unit_vectors[0], self.orienter.unit_vectors[1],
                na.array(self.width), self.transfer_function, self.sub_samples)
        return args

    def get_sampler(self, args):
        if self.use_light:
            if self.light_dir is None:
                self.set_default_light_dir()
            temp_dir = na.empty(3,dtype='float64')
            temp_dir = self.light_dir[0] * self.orienter.unit_vectors[1] + \
                    self.light_dir[1] * self.orienter.unit_vectors[2] + \
                    self.light_dir[2] * self.orienter.unit_vectors[0]
            if self.light_rgba is None:
                self.set_default_light_rgba()
            sampler = LightSourceRenderSampler(*args, light_dir=temp_dir,
                    light_rgba=self.light_rgba)
        else:
            sampler = self._sampler_object(*args)
        return sampler

    def finalize_image(self, image):
        pass

    def _render(self, double_check, num_threads, image, sampler):
        pbar = get_pbar("Ray casting", (self.volume.brick_dimensions + 1).prod(axis=-1).sum())
        total_cells = 0
        if double_check:
            for brick in self.volume.bricks:
                for data in brick.my_data:
                    if na.any(na.isnan(data)):
                        raise RuntimeError

        view_pos = self.front_center + self.orienter.unit_vectors[2] * 1.0e6 * self.width[2]
        for brick in self.volume.traverse(view_pos, self.front_center, image):
            sampler(brick, num_threads=num_threads)
            total_cells += na.prod(brick.my_data[0].shape)
            pbar.update(total_cells)

        pbar.finish()
        image = sampler.aimage
        self.finalize_image(image)
        return image

    def save_image(self, fn, clip_ratio, image):
        if self.comm.rank is 0 and fn is not None:
            if clip_ratio is not None:
                write_bitmap(image, fn, clip_ratio * image.std())
            else:
                write_bitmap(image, fn)


    def initialize_source(self):
        return self.volume.initialize_source()

    def snapshot(self, fn = None, clip_ratio = None, double_check = False,
                 num_threads = 0):
        r"""Ray-cast the camera.

        This method instructs the camera to take a snapshot -- i.e., call the ray
        caster -- based on its current settings.

        Parameters
        ----------
        fn : string, optional
            If supplied, the image will be saved out to this before being
            returned.  Scaling will be to the maximum value.
        clip_ratio : float, optional
            If supplied, the 'max_val' argument to write_bitmap will be handed
            clip_ratio * image.std()
        double_check : bool, optional
            Optionally makes sure that the data contains only valid entries.
            Used for debugging.
        num_threads : int, optional
            If supplied, will use 'num_threads' number of OpenMP threads during
            the rendering.  Defaults to 0, which uses the environment variable
            OMP_NUM_THREADS.

        Returns
        -------
        image : array
            An (N,M,3) array of the final returned values, in float64 form.
        """
        if num_threads is None:
            num_threads=get_num_threads()
        image = self.new_image()
        args = self.get_sampler_args(image)
        sampler = self.get_sampler(args)
        self.initialize_source()
        image = self._render(double_check, num_threads, image, sampler)
        self.save_image(fn, clip_ratio, image)
        return image

    def show(self, clip_ratio = None):
        r"""This will take a snapshot and display the resultant image in the
        IPython notebook.

        If yt is being run from within an IPython session, and it is able to
        determine this, this function will snapshot and send the resultant
        image to the IPython notebook for display.

        If yt can't determine if it's inside an IPython session, it will raise
        YTNotInsideNotebook.

        Parameters
        ----------
        clip_ratio : float, optional
            If supplied, the 'max_val' argument to write_bitmap will be handed
            clip_ratio * image.std()

        Examples
        --------

        >>> cam.show()

        """
        if "__IPYTHON__" in dir(__builtin__):
            from IPython.core.displaypub import publish_display_data
            image = self.snapshot()
            if clip_ratio is not None: clip_ratio *= image.std()
            data = write_bitmap(image, None, clip_ratio)
            publish_display_data(
                'yt.visualization.volume_rendering.camera.Camera',
                {'image/png' : data}
            )
        else:
            raise YTNotInsideNotebook


    def set_default_light_dir(self):
        self.light_dir = [1.,1.,1.]

    def set_default_light_rgba(self):
        self.light_rgba = [1.,1.,1.,1.]

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
        self._setup_box_properties(self.width, self.center, self.orienter.unit_vectors)

    def zoomin(self, final, n_steps, clip_ratio = None):
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
        clip_ratio : float, optional
            If supplied, the 'max_val' argument to write_bitmap will be handed
            clip_ratio * image.std()


        Examples
        --------

        >>> for i, snapshot in enumerate(cam.zoomin(100.0, 10)):
        ...     iw.write_bitmap(snapshot, "zoom_%04i.png" % i)
        """
        f = final**(1.0/n_steps)
        for i in xrange(n_steps):
            self.zoom(f)
            yield self.snapshot(clip_ratio = clip_ratio)

    def move_to(self, final, n_steps, final_width=None, exponential=False, clip_ratio = None):
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
        exponential : boolean
            Specifies whether the move/zoom transition follows an
            exponential path toward the destination or linear
        clip_ratio : float, optional
            If supplied, the 'max_val' argument to write_bitmap will be handed
            clip_ratio * image.std()
            
        Examples
        --------

        >>> for i, snapshot in enumerate(cam.move_to([0.2,0.3,0.6], 10)):
        ...     iw.write_bitmap(snapshot, "move_%04i.png" % i)
        """
        self.center = na.array(self.center)
        dW = None
        if exponential:
            if final_width is not None:
                if not iterable(final_width):
                    width = na.array([final_width, final_width, final_width]) 
                    # left/right, top/bottom, front/back 
                if (self.center == 0.0).all():
                    self.center += (na.array(final) - self.center) / (10. * n_steps)
                final_zoom = final_width/na.array(self.width)
                dW = final_zoom**(1.0/n_steps)
            else:
                dW = na.array([1.0,1.0,1.0])
            position_diff = (na.array(final)/self.center)*1.0
            dx = position_diff**(1.0/n_steps)
        else:
            if final_width is not None:
                if not iterable(final_width):
                    width = na.array([final_width, final_width, final_width]) 
                    # left/right, top/bottom, front/back
                dW = (1.0*final_width-na.array(self.width))/n_steps
            else:
                dW = na.array([0.0,0.0,0.0])
            dx = (na.array(final)-self.center)*1.0/n_steps
        for i in xrange(n_steps):
            if exponential:
                self.switch_view(center=self.center*dx, width=self.width*dW)
            else:
                self.switch_view(center=self.center+dx, width=self.width+dW)
            yield self.snapshot(clip_ratio = clip_ratio)

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
          
        R = get_rotation_matrix(theta, rot_vector)

        normal_vector = self.front_center-self.center

        self.switch_view(normal_vector=na.dot(R,normal_vector))

    def roll(self, theta):
        r"""Roll by a given angle

        Roll the view.

        Parameters
        ----------
        theta : float, in radians
             Angle (in radians) by which to roll the view.

        Examples
        --------

        >>> cam.roll(na.pi/4)
        """
        rot_vector = self.orienter.normal_vector
        R = get_rotation_matrix(theta, rot_vector)
        north_vector = self.orienter.north_vector
        self.switch_view(north_vector=na.dot(R, north_vector))

    def rotation(self, theta, n_steps, rot_vector=None, clip_ratio = None):
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
        clip_ratio : float, optional
            If supplied, the 'max_val' argument to write_bitmap will be handed
            clip_ratio * image.std()

        Examples
        --------

        >>> for i, snapshot in enumerate(cam.rotation(na.pi, 10)):
        ...     iw.write_bitmap(snapshot, 'rotation_%04i.png' % i)
        """

        dtheta = (1.0*theta)/n_steps
        for i in xrange(n_steps):
            self.rotate(dtheta, rot_vector=rot_vector)
            yield self.snapshot(clip_ratio = clip_ratio)

data_object_registry["camera"] = Camera

class InteractiveCamera(Camera):
    frames = []

    def snapshot(self, fn=None, clip_ratio=None):
        import matplotlib.pylab as pylab
        pylab.figure(2)
        self.transfer_function.show()
        pylab.draw()
        im = Camera.snapshot(self, fn, clip_ratio)
        pylab.figure(1)
        pylab.imshow(im / im.max())
        pylab.draw()
        self.frames.append(im)

    def rotation(self, theta, n_steps, rot_vector=None):
        for frame in Camera.rotation(self, theta, n_steps, rot_vector):
            if frame is not None:
                self.frames.append(frame)
                
    def zoomin(self, final, n_steps):
        for frame in Camera.zoomin(self, final, n_steps):
            if frame is not None:
                self.frames.append(frame)
                
    def clear_frames(self):
        del self.frames
        self.frames = []
        
    def save_frames(self, basename, clip_ratio=None):
        for i, frame in enumerate(self.frames):
            fn = basename + '_%04i.png'%i
            if clip_ratio is not None:
                write_bitmap(frame, fn, clip_ratio*image.std())
            else:
                write_bitmap(frame, fn)

data_object_registry["interactive_camera"] = InteractiveCamera

class PerspectiveCamera(Camera):
    expand_factor = 1.0
    def __init__(self, *args, **kwargs):
        expand_factor = kwargs.pop('expand_factor', 1.0)
        Camera.__init__(self, *args, **kwargs)

    def get_sampler_args(self, image):
        # We should move away from pre-generation of vectors like this and into
        # the usage of on-the-fly generation in the VolumeIntegrator module
        # We might have a different width and back_center
        dl = (self.back_center - self.front_center)
        self.front_center += self.expand_factor*dl
        self.back_center -= dl

        px = na.linspace(-self.width[0]/2.0, self.width[0]/2.0,
                         self.resolution[0])[:,None]
        py = na.linspace(-self.width[1]/2.0, self.width[1]/2.0,
                         self.resolution[1])[None,:]
        inv_mat = self.orienter.inv_mat
        positions = na.zeros((self.resolution[0], self.resolution[1], 3),
                          dtype='float64', order='C')
        positions[:,:,0] = inv_mat[0,0]*px+inv_mat[0,1]*py+self.back_center[0]
        positions[:,:,1] = inv_mat[1,0]*px+inv_mat[1,1]*py+self.back_center[1]
        positions[:,:,2] = inv_mat[2,0]*px+inv_mat[2,1]*py+self.back_center[2]
        bounds = (px.min(), px.max(), py.min(), py.max())

        # We are likely adding on an odd cutting condition here
        vectors = self.front_center - positions
        positions = self.front_center - 1.0*(((self.back_center-self.front_center)**2).sum())**0.5*vectors
        vectors = (self.front_center - positions)

        uv = na.ones(3, dtype='float64')
        image.shape = (self.resolution[0]**2,1,3)
        vectors.shape = (self.resolution[0]**2,1,3)
        positions.shape = (self.resolution[0]**2,1,3)
        args = (positions, vectors, self.back_center, 
                (0.0,1.0,0.0,1.0),
                image, uv, uv,
                na.zeros(3, dtype='float64'), 
                self.transfer_function, self.sub_samples)
        return args

    def finalize_image(self, image):
        image.shape = self.resolution[0], self.resolution[0], 3

def corners(left_edge, right_edge):
    return na.array([
      [left_edge[:,0], left_edge[:,1], left_edge[:,2]],
      [right_edge[:,0], left_edge[:,1], left_edge[:,2]],
      [right_edge[:,0], right_edge[:,1], left_edge[:,2]],
      [right_edge[:,0], right_edge[:,1], right_edge[:,2]],
      [left_edge[:,0], right_edge[:,1], right_edge[:,2]],
      [left_edge[:,0], left_edge[:,1], right_edge[:,2]],
      [right_edge[:,0], left_edge[:,1], right_edge[:,2]],
      [left_edge[:,0], right_edge[:,1], left_edge[:,2]],
    ], dtype='float64')

class HEALpixCamera(Camera):
    def __init__(self, center, radius, nside,
                 transfer_function = None, fields = None,
                 sub_samples = 5, log_fields = None, volume = None,
                 pf = None, use_kd=True, no_ghost=False, use_light=False):
        ParallelAnalysisInterface.__init__(self)
        if pf is not None: self.pf = pf
        self.center = na.array(center, dtype='float64')
        self.radius = radius
        self.nside = nside
        self.use_kd = use_kd
        if transfer_function is None:
            transfer_function = ProjectionTransferFunction()
        self.transfer_function = transfer_function
        if fields is None: fields = ["Density"]
        self.fields = fields
        self.sub_samples = sub_samples
        self.log_fields = log_fields
        self.use_light = use_light
        self.light_dir = None
        self.light_rgba = None
        if volume is None:
            volume = AMRKDTree(self.pf, fields=self.fields, no_ghost=no_ghost,
                               log_fields=log_fields)
        self.use_kd = isinstance(volume, AMRKDTree)
        self.volume = volume

    def new_image(self):
        image = na.zeros((12 * self.nside ** 2, 1, 3), dtype='float64', order='C')
        return image

    def get_sampler_args(self, image):
        nv = 12 * self.nside ** 2
        vs = arr_pix2vec_nest(self.nside, na.arange(nv))
        vs *= self.radius
        vs.shape = nv, 1, 3
        uv = na.ones(3, dtype='float64')
        positions = na.ones((nv, 1, 3), dtype='float64') * self.center
        args = (positions, vs, self.center,
                (0.0, 1.0, 0.0, 1.0),
                image, uv, uv,
                na.zeros(3, dtype='float64'),
                self.transfer_function, self.sub_samples)
        return args
 

    def _render(self, double_check, num_threads, image, sampler):
        pbar = get_pbar("Ray casting", (self.volume.brick_dimensions + 1).prod(axis=-1).sum())
        total_cells = 0
        if double_check:
            for brick in self.volume.bricks:
                for data in brick.my_data:
                    if na.any(na.isnan(data)):
                        raise RuntimeError
        
        view_pos = self.center
        for brick in self.volume.traverse(view_pos, None, image):
            sampler(brick, num_threads=num_threads)
            total_cells += na.prod(brick.my_data[0].shape)
            pbar.update(total_cells)
        
        pbar.finish()
        image = sampler.aimage

        self.finalize_image(image)

        return image

    def snapshot(self, fn = None, clip_ratio = None, double_check = False,
                 num_threads = 0, clim = None, label = None):
        r"""Ray-cast the camera.

        This method instructs the camera to take a snapshot -- i.e., call the ray
        caster -- based on its current settings.

        Parameters
        ----------
        fn : string, optional
            If supplied, the image will be saved out to this before being
            returned.  Scaling will be to the maximum value.
        clip_ratio : float, optional
            If supplied, the 'max_val' argument to write_bitmap will be handed
            clip_ratio * image.std()

        Returns
        -------
        image : array
            An (N,M,3) array of the final returned values, in float64 form.
        """
        if num_threads is None:
            num_threads=get_num_threads()
        image = self.new_image()
        args = self.get_sampler_args(image)
        sampler = self.get_sampler(args)
        self.volume.initialize_source()
        image = self._render(double_check, num_threads, image, sampler)
        self.save_image(fn, clim, image, label = label)
        return image

    def save_image(self, fn, clim, image, label = None):
        if self.comm.rank is 0 and fn is not None:
            # This assumes Density; this is a relatively safe assumption.
            import matplotlib.figure
            import matplotlib.backends.backend_agg
            phi, theta = na.mgrid[0.0:2*na.pi:800j, 0:na.pi:800j]
            pixi = arr_ang2pix_nest(self.nside, theta.ravel(), phi.ravel())
            image *= self.radius * self.pf['cm']
            img = na.log10(image[:,0,0][pixi]).reshape((800,800))

            fig = matplotlib.figure.Figure((10, 5))
            ax = fig.add_subplot(1,1,1,projection='hammer')
            implot = ax.imshow(img, extent=(-na.pi,na.pi,-na.pi/2,na.pi/2), clip_on=False, aspect=0.5)
            cb = fig.colorbar(implot, orientation='horizontal')

            if label == None:
                cb.set_label("Projected %s" % self.fields[0])
            else:
                cb.set_label(label)
            if clim is not None: cb.set_clim(*clim)
            ax.xaxis.set_ticks(())
            ax.yaxis.set_ticks(())
            canvas = matplotlib.backends.backend_agg.FigureCanvasAgg(fig)
            canvas.print_figure(fn)


class AdaptiveHEALpixCamera(Camera):
    def __init__(self, center, radius, nside,
                 transfer_function = None, fields = None,
                 sub_samples = 5, log_fields = None, volume = None,
                 pf = None, use_kd=True, no_ghost=False,
                 rays_per_cell = 0.1, max_nside = 8192):
        ParallelAnalysisInterface.__init__(self)
        if pf is not None: self.pf = pf
        self.center = na.array(center, dtype='float64')
        self.radius = radius
        self.use_kd = use_kd
        if transfer_function is None:
            transfer_function = ProjectionTransferFunction()
        self.transfer_function = transfer_function
        if fields is None: fields = ["Density"]
        self.fields = fields
        self.sub_samples = sub_samples
        self.log_fields = log_fields
        if volume is None:
            volume = AMRKDTree(self.pf, fields=self.fields, no_ghost=no_ghost,
                               log_fields=log_fields)
        self.use_kd = isinstance(volume, AMRKDTree)
        self.volume = volume
        self.initial_nside = nside
        self.rays_per_cell = rays_per_cell
        self.max_nside = max_nside

    def snapshot(self, fn = None):
        tfp = TransferFunctionProxy(self.transfer_function)
        tfp.ns = self.sub_samples
        self.volume.initialize_source()
        mylog.info("Adaptively rendering.")
        pbar = get_pbar("Ray casting",
                        (self.volume.brick_dimensions + 1).prod(axis=-1).sum())
        total_cells = 0
        bricks = [b for b in self.volume.traverse(None, self.center, None)][::-1]
        left_edges = na.array([b.LeftEdge for b in bricks])
        right_edges = na.array([b.RightEdge for b in bricks])
        min_dx = min(((b.RightEdge[0] - b.LeftEdge[0])/b.my_data[0].shape[0]
                     for b in bricks))
        # We jitter a bit if we're on a boundary of our initial grid
        for i in range(3):
            if bricks[0].LeftEdge[i] == self.center[i]:
                self.center += 1e-2 * min_dx
            elif bricks[0].RightEdge[i] == self.center[i]:
                self.center -= 1e-2 * min_dx
        ray_source = AdaptiveRaySource(self.center, self.rays_per_cell,
                                       self.initial_nside, self.radius,
                                       bricks, left_edges, right_edges, self.max_nside)
        for i,brick in enumerate(bricks):
            ray_source.integrate_brick(brick, tfp, i, left_edges, right_edges,
                                       bricks)
            total_cells += na.prod(brick.my_data[0].shape)
            pbar.update(total_cells)
        pbar.finish()
        info, values = ray_source.get_rays()
        return info, values

class StereoPairCamera(Camera):
    def __init__(self,original_camera,
                 auto_focus=False,
                 focal_length=None,
                 frac_near_plane = 0.90, 
                 frac_far_plane  = 1.10,
                 frac_eye_separation=0.05,
                 aperture = 60.0,
                 relative_separation=0.005):
        """
        Auto-focus is adapted from a guide & code at :
        http://paulbourke.net/miscellaneous/stereographics/stereorender/
        """
        ParallelAnalysisInterface.__init__(self)
        self.original_camera = original_camera
        oc = self.original_camera
        if self.auto_focus:
            dist = lambda x,y: na.sqrt(na.sum((x-y)**2.0))
            if self.focal_length is None:
                self.focal_length = dist(oc.normal_vector,0.0)
            self.focal_far  = oc.center + frac_far_plane*oc.normal_vector
            self.focal_near = oc.center + frac_near_plane*oc.normal_vector
            self.wh_ratio = oc.resolution[0]/oc.resolution[1]
            self.eye_sep  = self.focal_length*frac_eye_separation
            self.aperture = aperture
            self.frac_eye_separation = frac_eye_separation
            self.center_eye_pos = oc.center + oc.normal_vector
        else:
            #default to old separation
            self.relative_separation = relative_separation
    
    def finalize_image(self,image):
        if self.auto_focus:
            #we have extra frustum pixels on the left and right
            #cameras
            left_trim,right_trim = self.trim[0],self.trim[1]
            left = abs(left_trim)
            right = image.shapae[0]-abs(right_trim)
            image = image[left:right,:]
            return image

	def auto_split(self):
		"""We must calculate the new camera centers, as well
        as the extended frustum pixels."""
        oc = self.original_camera
        nv = oc.orienter.normal_vector
        up = oc.north_vector
        c = oc.center
        px = resolution[0] #pixel width
        norm = lambda x: na.sqrt(na.dot(x,x.conj()))
        between_eyes = na.cross(nv,up)
        between_eyes /= norm(between_eyes)
        between_eyes *= eye_sep/2.0
        le_norm = nv-between_eyes 
        le_c= c-between_eyes 
        re_norm = nv+between_eyes 
        re_c = c+between_eyes 
        angular_aperture = na.tan(self.aperture/360.0*2.0*na.pi/2.0)
        delta = na.rint(px*self.frac_eye_separation/(2.0*(angular_aperture)))
        delta = delta.astype('int')
        eresolution = resolution[0]+delta
        left_camera = Camera(le_c, le_norm, oc.width,
                     eresolution, oc.transfer_function, north_vector=up,
                     volume=oc.volume, fields=oc.fields, 
                     log_fields=oc.log_fields,
                     sub_samples=oc.sub_samples, pf=oc.pf)
        left_camera.trim = [-delta,0]
        right_camera = Camera(re_c, re_norm, oc.width,
                     eresolution, oc.transfer_function, north_vector=up,
                     volume=oc.volume, fields=oc.fields, 
                     log_fields=oc.log_fields,
                     sub_samples=oc.sub_samples, pf=oc.pf)
        right_camera.trim = [0,-delta]
        return (left_camera, right_camera)

    def split(self):
        if self.auto_focus:
            return self.auto_split()
        else:
            return self.default_split()
    
    def default_split(self):
        oc = self.original_camera
        uv = oc.orienter.unit_vectors
        c = oc.center
        fc = oc.front_center
        wx, wy, wz = oc.width
        left_normal = fc + uv[1] * 0.5*self.relative_separation * wx - c
        right_normal = fc - uv[1] * 0.5*self.relative_separation * wx - c
        left_camera = Camera(c, left_normal, oc.width,
                             oc.resolution, oc.transfer_function, north_vector=uv[0],
                             volume=oc.volume, fields=oc.fields, log_fields=oc.log_fields,
                             sub_samples=oc.sub_samples, pf=oc.pf)
        right_camera = Camera(c, right_normal, oc.width,
                             oc.resolution, oc.transfer_function, north_vector=uv[0],
                             volume=oc.volume, fields=oc.fields, log_fields=oc.log_fields,
                             sub_samples=oc.sub_samples, pf=oc.pf)
        return (left_camera, right_camera)



        

class FisheyeCamera(Camera):
    def __init__(self, center, radius, fov, resolution,
                 transfer_function = None, fields = None,
                 sub_samples = 5, log_fields = None, volume = None,
                 pf = None, no_ghost=False, rotation = None, use_light=False):
        ParallelAnalysisInterface.__init__(self)
        self.use_light = use_light
        self.light_dir = None
        self.light_rgba = None
        if rotation is None: rotation = na.eye(3)
        self.rotation_matrix = rotation
        if pf is not None: self.pf = pf
        self.center = na.array(center, dtype='float64')
        self.radius = radius
        self.fov = fov
        if iterable(resolution):
            raise RuntimeError("Resolution must be a single int")
        self.resolution = resolution
        if transfer_function is None:
            transfer_function = ProjectionTransferFunction()
        self.transfer_function = transfer_function
        if fields is None: fields = ["Density"]
        self.fields = fields
        self.sub_samples = sub_samples
        self.log_fields = log_fields
        if volume is None:
            volume = AMRKDTree(self.pf, fields=self.fields, no_ghost=no_ghost,
                               log_fields=log_fields)
        self.volume = volume

    def new_image(self):
        image = na.zeros((self.resolution**2,1,3), dtype='float64', order='C')
        return image
        
    def get_sampler_args(self, image):
        vp = arr_fisheye_vectors(self.resolution, self.fov)
        vp.shape = (self.resolution**2,1,3)
        vp2 = vp.copy()
        for i in range(3):
            vp[:,:,i] = (vp2 * self.rotation_matrix[:,i]).sum(axis=2)
        del vp2
        vp *= self.radius
        uv = na.ones(3, dtype='float64')
        positions = na.ones((self.resolution**2, 1, 3), dtype='float64') * self.center

        args = (positions, vp, self.center,
                (0.0, 1.0, 0.0, 1.0),
                image, uv, uv,
                na.zeros(3, dtype='float64'),
                self.transfer_function, self.sub_samples)
        return args


    def finalize_image(self, image):
        image.shape = self.resolution, self.resolution, 3

    def _render(self, double_check, num_threads, image, sampler):
        pbar = get_pbar("Ray casting", (self.volume.brick_dimensions + 1).prod(axis=-1).sum())
        total_cells = 0
        if double_check:
            for brick in self.volume.bricks:
                for data in brick.my_data:
                    if na.any(na.isnan(data)):
                        raise RuntimeError
        
        view_pos = self.center
        for brick in self.volume.traverse(view_pos, None, image):
            sampler(brick, num_threads=num_threads)
            total_cells += na.prod(brick.my_data[0].shape)
            pbar.update(total_cells)
        
        pbar.finish()
        image = sampler.aimage

        self.finalize_image(image)

        return image

class MosaicFisheyeCamera(Camera):
    def __init__(self, center, radius, fov, resolution, focal_center=None,
                 transfer_function=None, fields=None,
                 sub_samples=5, log_fields=None, volume=None,
                 pf=None, l_max=None, no_ghost=False,nimx=1, nimy=1, procs_per_wg=None,
                 rotation=None):
        r"""A fisheye lens camera, taking adantage of image plane decomposition
        for parallelism..

        The camera represents the eye of an observer, which will be used to
        generate ray-cast volume renderings of the domain. In this case, the
        rays are defined by a fisheye lens

        Parameters
        ----------
        center : array_like
            The current "center" of the observer, from which the rays will be
            cast
        radius : float
            The radial distance to cast to
        resolution : int
            The number of pixels in each direction.  Must be a single int.
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
        pf : `~yt.data_objects.api.StaticOutput`
            For now, this is a require parameter!  But in the future it will become
            optional.  This is the parameter file to volume render.
        l_max: int, optional
            Specifies the maximum level to be rendered.  Also
            specifies the maximum level used in the AMRKDTree
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
        nimx: int, optional
            The number by which to decompose the image plane into in the x
            direction.  Must evenly divide the resolution.
        nimy: int, optional
            The number by which to decompose the image plane into in the y 
            direction.  Must evenly divide the resolution.
        procs_per_wg: int, optional
            The number of processors to use on each sub-image. Within each
            subplane, the volume will be decomposed using the AMRKDTree with
            procs_per_wg processors.  

        Notes
        -----
            The product of nimx*nimy*procs_per_wg must be equal to or less than
            the total number of mpi processes.  

            Unlike the non-Mosaic camera, this will only return each sub-image
            to the root processor of each sub-image workgroup in order to save
            memory.  To save the final image, one must then call
            MosaicFisheyeCamera.save_image('filename')

        Examples
        --------

        >>> from yt.mods import *
        
        >>> pf = load('DD1717')
        
        >>> N = 512 # Pixels (1024^2)
        >>> c = (pf.domain_right_edge + pf.domain_left_edge)/2. # Center
        >>> radius = (pf.domain_right_edge - pf.domain_left_edge)/2.
        >>> fov = 180.0
        
        >>> field='Density'
        >>> mi,ma = pf.h.all_data().quantities['Extrema']('Density')[0]
        >>> mi,ma = na.log10(mi), na.log10(ma)
        
        # You may want to comment out the above lines and manually set the min and max
        # of the log of the Density field. For example:
        # mi,ma = -30.5,-26.5
        
        # Another good place to center the camera is close to the maximum density.
        # v,c = pf.h.find_max('Density')
        # c -= 0.1*radius
        
       
        # Construct transfer function
        >>> tf = ColorTransferFunction((mi-1, ma+1),nbins=1024)
        
        # Sample transfer function with Nc gaussians.  Use col_bounds keyword to limit
        # the color range to the min and max values, rather than the transfer function
        # bounds.
        >>> Nc = 5
        >>> tf.add_layers(Nc,w=0.005, col_bounds = (mi,ma), alpha=na.logspace(-2,0,Nc),
        >>>         colormap='RdBu_r')
        >>> 
        # Create the camera object. Use the keyword: no_ghost=True if a lot of time is
        # spent creating vertex-centered data. In this case I'm running with 8
        # processors, and am splitting the image plane into 4 pieces and using 2
        # processors on each piece.
        >>> cam = MosaicFisheyeCamera(c, radius, fov, N,
        >>>         transfer_function = tf, 
        >>>         sub_samples = 5, 
        >>>         pf=pf, 
        >>>         nimx=2,nimy=2,procs_per_wg=2)
        
        # Take a snapshot
        >>> im = cam.snapshot()
        
        # Save the image
        >>> cam.save_image('fisheye_mosaic.png')

        """

        ParallelAnalysisInterface.__init__(self)
        self.image_decomp = self.comm.size>1
        if self.image_decomp:
            PP = ProcessorPool()
            npatches = nimy*nimx
            if procs_per_wg is None:
                if (PP.size % npatches):
                    raise RuntimeError("Cannot evenly divide %i procs to %i patches" % (PP.size,npatches))
                else:
                    procs_per_wg = PP.size / npatches
            if (PP.size != npatches*procs_per_wg):
               raise RuntimeError("You need %i processors to utilize %i procs per one patch in [%i,%i] grid" 
                     % (npatches*procs_per_wg,procs_per_wg,nimx,nimy))
 
            for j in range(nimy):
                for i in range(nimx):
                    PP.add_workgroup(size=procs_per_wg, name='%04i_%04i'%(i,j))
                    
            for wg in PP.workgroups:
                if self.comm.rank in wg.ranks:
                    my_wg = wg
            
            self.global_comm = self.comm
            self.comm = my_wg.comm
            self.wg = my_wg
            self.imi = int(self.wg.name[0:4])
            self.imj = int(self.wg.name[5:9])
            mylog.info('My new communicator has the name %s' % self.wg.name)
            self.nimx = nimx
            self.nimy = nimy
        else:
            self.imi = 0
            self.imj = 0
            self.nimx = 1
            self.nimy = 1
        if pf is not None: self.pf = pf
        
        if rotation is None: rotation = na.eye(3)
        self.rotation_matrix = rotation
        
        self.normal_vector = na.array([0.,0.,1])
        self.north_vector = na.array([1.,0.,0.])
        self.east_vector = na.array([0.,1.,0.])
        self.rotation_vector = self.north_vector

        if iterable(resolution):
            raise RuntimeError("Resolution must be a single int")
        self.resolution = resolution
        self.center = na.array(center, dtype='float64')
        self.focal_center = focal_center
        self.radius = radius
        self.fov = fov
        if transfer_function is None:
            transfer_function = ProjectionTransferFunction()
        self.transfer_function = transfer_function
        if fields is None: fields = ["Density"]
        self.fields = fields
        self.sub_samples = sub_samples
        self.log_fields = log_fields
        if volume is None:
            volume = AMRKDTree(self.pf, fields=self.fields, no_ghost=no_ghost,
                               log_fields=log_fields,l_max=l_max)
        self.volume = volume
        self.vp = None
        self.image = None 

    def get_vector_plane(self):
        if self.focal_center is not None:
            rvec =  na.array(self.focal_center) - na.array(self.center)
            rvec /= (rvec**2).sum()**0.5
            angle = na.arccos( (self.normal_vector*rvec).sum()/( (self.normal_vector**2).sum()**0.5 *
                (rvec**2).sum()**0.5))
            rot_vector = na.cross(rvec, self.normal_vector)
            rot_vector /= (rot_vector**2).sum()**0.5
            
            self.rotation_matrix = get_rotation_matrix(angle,rot_vector)
            self.normal_vector = na.dot(self.rotation_matrix,self.normal_vector)
            self.north_vector = na.dot(self.rotation_matrix,self.north_vector)
            self.east_vector = na.dot(self.rotation_matrix,self.east_vector)
        else:
            self.focal_center = self.center + self.radius*self.normal_vector  
        dist = ((self.focal_center - self.center)**2).sum()**0.5
        # We now follow figures 4-7 of:
        # http://paulbourke.net/miscellaneous/domefisheye/fisheye/
        # ...but all in Cython.
        
        self.vp = arr_fisheye_vectors(self.resolution, self.fov, self.nimx, 
                self.nimy, self.imi, self.imj)
        
        self.vp = rotate_vectors(self.vp, self.rotation_matrix)

        self.center = self.focal_center - dist*self.normal_vector
        self.vp *= self.radius
        nx, ny = self.vp.shape[0], self.vp.shape[1]
        self.vp.shape = (nx*ny,1,3)

    def snapshot(self):
        if self.vp is None:
            self.get_vector_plane()

        nx,ny = self.resolution/self.nimx, self.resolution/self.nimy
        image = na.zeros((nx*ny,1,3), dtype='float64', order='C')
        uv = na.ones(3, dtype='float64')
        positions = na.ones((nx*ny, 1, 3), dtype='float64') * self.center
        vector_plane = VectorPlane(positions, self.vp, self.center,
                        (0.0, 1.0, 0.0, 1.0), image, uv, uv)
        tfp = TransferFunctionProxy(self.transfer_function)
        tfp.ns = self.sub_samples
        self.volume.initialize_source()
        mylog.info("Rendering fisheye of %s^2", self.resolution)
        pbar = get_pbar("Ray casting",
                        (self.volume.brick_dimensions + 1).prod(axis=-1).sum())

        total_cells = 0
        for brick in self.volume.traverse(None, self.center, image):
            brick.cast_plane(tfp, vector_plane)
            total_cells += na.prod(brick.my_data[0].shape)
            pbar.update(total_cells)
        pbar.finish()
        image.shape = (nx, ny, 3)

        if self.image is not None:
            del self.image
        self.image = image
       
        return image

    def save_image(self, fn, clip_ratio=None):
        if '.png' not in fn:
            fn = fn + '.png'
        
        try:
            image = self.image
        except:
            mylog.error('You must first take a snapshot')
            raise(UserWarning)
        
        image = self.image
        nx,ny = self.resolution/self.nimx, self.resolution/self.nimy
        if self.image_decomp:
            if self.comm.rank == 0:
                if self.global_comm.rank == 0:
                    final_image = na.empty((nx*self.nimx, 
                        ny*self.nimy, 3),
                        dtype='float64',order='C')
                    final_image[:nx, :ny, :] = image
                    for j in range(self.nimy):
                        for i in range(self.nimx):
                            if i==0 and j==0: continue
                            arr = self.global_comm.recv_array((self.wg.size)*(j*self.nimx + i), tag = (self.wg.size)*(j*self.nimx + i))

                            final_image[i*nx:(i+1)*nx, j*ny:(j+1)*ny,:] = arr
                            del arr
                    if clip_ratio is not None:
                        write_bitmap(final_image, fn, clip_ratio*final_image.std())
                    else:
                        write_bitmap(final_image, fn)
                else:
                    self.global_comm.send_array(image, 0, tag = self.global_comm.rank)
        else:
            if self.comm.rank == 0:
                if clip_ratio is not None:
                    write_bitmap(image, fn, clip_ratio*image.std())
                else:
                    write_bitmap(image, fn)
        return

    def rotate(self, theta, rot_vector=None, keep_focus=True):
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
            rot_vector = self.north_vector

        dist = ((self.focal_center - self.center)**2).sum()**0.5

        R = get_rotation_matrix(theta, rot_vector)

        self.vp = rotate_vectors(self.vp, R)
        self.normal_vector = na.dot(R,self.normal_vector)
        self.north_vector = na.dot(R,self.north_vector)
        self.east_vector = na.dot(R,self.east_vector)

        if keep_focus:
            self.center = self.focal_center - dist*self.normal_vector

    def rotation(self, theta, n_steps, rot_vector=None, keep_focus=True):
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
        ...     iw.write_bitmap(snapshot, 'rotation_%04i.png' % i)
        """

        dtheta = (1.0*theta)/n_steps
        for i in xrange(n_steps):
            self.rotate(dtheta, rot_vector=rot_vector, keep_focus=keep_focus)
            yield self.snapshot()

    def move_to(self,final,n_steps,exponential=False):
        r"""Loop over a look_at

        This will yield `n_steps` snapshots until the current view has been
        moved to a final center of `final`.

        Parameters
        ----------
        final : array_like
            The final center to move to after `n_steps`
        n_steps : int
            The number of look_at snapshots to make.
        exponential : boolean
            Specifies whether the move/zoom transition follows an
            exponential path toward the destination or linear

        Examples
        --------

        >>> for i, snapshot in enumerate(cam.move_to([0.2,0.3,0.6], 10)):
        ...     cam.save_image('move_%04i.png' % i)
        """
        if exponential:
            position_diff = (na.array(final)/self.center)*1.0
            dx = position_diff**(1.0/n_steps)
        else:
            dx = (na.array(final) - self.center)*1.0/n_steps
        for i in xrange(n_steps):
            if exponential:
                self.center *= dx
            else:
                self.center += dx
            yield self.snapshot()

def allsky_projection(pf, center, radius, nside, field, weight = None,
                      inner_radius = 10, rotation = None):
    r"""Project through a parameter file, through an allsky-method
    decomposition from HEALpix, and return the image plane.

    This function will accept the necessary items to integrate through a volume
    over 4pi and return the integrated field of view to the user.  Note that if
    a weight is supplied, it will multiply the pre-interpolated values
    together.

    Parameters
    ----------
    pf : `~yt.data_objects.api.StaticOutput`
        This is the parameter file to volume render.
    center : array_like
        The current "center" of the view port -- the focal point for the
        camera.
    radius : float or list of floats
        The radius to integrate out to of the image.
    nside : int
        The HEALpix degree.  The number of rays integrated is 12*(Nside**2)
        Must be a power of two!
    field : string
        The field to project through the volume
    weight : optional, default None
        If supplied, the field will be pre-multiplied by this, then divided by
        the integrated value of this field.  This returns an average rather
        than a sum.
    inner_radius : optional, float, defaults to 0.05
        The radius of the inner clipping plane, in units of the dx at the point
        at which the volume rendering is centered.  This avoids unphysical
        effects of nearby cells.
    rotation : optional, 3x3 array
        If supplied, the vectors will be rotated by this.  You can construct
        this by, for instance, calling na.array([v1,v2,v3]) where those are the
        three reference planes of an orthogonal frame (see ortho_find).

    Returns
    -------
    image : array
        An ((Nside**2)*12,1,3) array of the final integrated values, in float64 form.

    Examples
    --------

    >>> image = allsky_projection(pf, [0.5, 0.5, 0.5], 1.0/pf['mpc'],
                      32, "Temperature", "Density")
    >>> plot_allsky_healpix(image, 32, "healpix.png")

    """
    # We manually modify the ProjectionTransferFunction to get it to work the
    # way we want, with a second field that's also passed through.
    fields = [field]
    center = na.array(center, dtype='float64')
    if weight is not None:
        # This is a temporary field, which we will remove at the end.
        def _make_wf(f, w):
            def temp_weightfield(a, b):
                tr = b[f].astype("float64") * b[w]
                return tr
            return temp_weightfield
        pf.field_info.add_field("temp_weightfield",
            function=_make_wf(field, weight))
        fields = ["temp_weightfield", weight]
    nv = 12*nside**2
    image = na.zeros((nv,1,3), dtype='float64', order='C')
    vs = arr_pix2vec_nest(nside, na.arange(nv))
    vs.shape = (nv,1,3)
    if rotation is not None:
        vs2 = vs.copy()
        for i in range(3):
            vs[:,:,i] = (vs2 * rotation[:,i]).sum(axis=2)
    else:
        vs += 1e-8
    positions = na.ones((nv, 1, 3), dtype='float64', order='C') * center
    dx = min(g.dds.min() for g in pf.h.find_point(center)[0])
    positions += inner_radius * dx * vs
    vs *= radius
    uv = na.ones(3, dtype='float64')
    grids = pf.h.sphere(center, radius)._grids
    sampler = ProjectionSampler(positions, vs, center, (0.0, 0.0, 0.0, 0.0),
                                image, uv, uv, na.zeros(3, dtype='float64'))
    pb = get_pbar("Sampling ", len(grids))
    for i,grid in enumerate(grids):
        data = [grid[field] * grid.child_mask.astype('float64')
                for field in fields]
        pg = PartitionedGrid(
            grid.id, data,
            grid.LeftEdge, grid.RightEdge,
            grid.ActiveDimensions.astype("int64"))
        grid.clear_data()
        sampler(pg)
        pb.update(i)
    pb.finish()
    image = sampler.aimage
    if weight is None:
        dl = radius * pf.units[pf.field_info[field].projection_conversion]
        image *= dl
    else:
        image[:,:,0] /= image[:,:,1]
        pf.field_info.pop("temp_weightfield")
        for g in pf.h.grids:
            if "temp_weightfield" in g.keys():
                del g["temp_weightfield"]
    return image[:,0,0]

def plot_allsky_healpix(image, nside, fn, label = "", rotation = None,
                        take_log = True, resolution=512, cmin=None, cmax=None):
    import matplotlib.figure
    import matplotlib.backends.backend_agg
    if rotation is None: rotation = na.eye(3).astype("float64")

    img, count = pixelize_healpix(nside, image, resolution, resolution, rotation)

    fig = matplotlib.figure.Figure((10, 5))
    ax = fig.add_subplot(1,1,1,projection='aitoff')
    if take_log: func = na.log10
    else: func = lambda a: a
    implot = ax.imshow(func(img), extent=(-na.pi,na.pi,-na.pi/2,na.pi/2),
                       clip_on=False, aspect=0.5, vmin=cmin, vmax=cmax)
    cb = fig.colorbar(implot, orientation='horizontal')
    cb.set_label(label)
    ax.xaxis.set_ticks(())
    ax.yaxis.set_ticks(())
    canvas = matplotlib.backends.backend_agg.FigureCanvasAgg(fig)
    canvas.print_figure(fn)
    return img, count

class ProjectionCamera(Camera):
    def __init__(self, center, normal_vector, width, resolution,
            field, weight=None, volume=None, no_ghost = False, 
            le=None, re=None,
            north_vector=None, pf=None, interpolated=False):

        if not interpolated:
            volume = 1

        self.interpolated = interpolated
        self.field = field
        self.weight = weight
        self.resolution = resolution

        fields = [field]
        if self.weight is not None:
            # This is a temporary field, which we will remove at the end.
            def _make_wf(f, w):
                def temp_weightfield(a, b):
                    tr = b[f].astype("float64") * b[w]
                    return tr
                return temp_weightfield
            pf.field_info.add_field("temp_weightfield",
                function=_make_wf(self.field, self.weight))
            fields = ["temp_weightfield", self.weight]
        
        self.fields = fields
        self.log_fields = [False]*len(self.fields)
        Camera.__init__(self, center, normal_vector, width, resolution, None,
                fields = fields, pf=pf, volume=volume,
                log_fields=self.log_fields, 
                le=le, re=re, north_vector=north_vector,
                no_ghost=no_ghost)

    def get_sampler(self, args):
        if self.interpolated:
            sampler = InterpolatedProjectionSampler(*args)
        else:
            sampler = ProjectionSampler(*args)
        return sampler

    def initialize_source(self):
        if self.interpolated:
            Camera.initialize_source(self)
        else:
            pass

    def get_sampler_args(self, image):
        rotp = na.concatenate([self.orienter.inv_mat.ravel('F'), self.back_center.ravel()])
        args = (rotp, self.box_vectors[2], self.back_center,
            (-self.width[0]/2, self.width[0]/2,
             -self.width[1]/2, self.width[1]/2),
            image, self.orienter.unit_vectors[0], self.orienter.unit_vectors[1],
                na.array(self.width), self.sub_samples)
        return args

    def finalize_image(self,image):
        pf = self.pf
        if self.weight is None:
            dl = self.width[2] * pf.units[pf.field_info[self.field].projection_conversion]
            image *= dl
        else:
            image[:,:,0] /= image[:,:,1]
        return image[:,:,0]


    def _render(self, double_check, num_threads, image, sampler):
        # Calculate the eight corners of the box
        # Back corners ...
        if self.interpolated:
            return Camera._render(self, double_check, num_threads, image,
                    sampler)
        pf = self.pf
        width = self.width[2]
        north_vector = self.orienter.unit_vectors[0]
        east_vector = self.orienter.unit_vectors[1]
        normal_vector = self.orienter.unit_vectors[2]
        fields = self.fields

        mi = pf.domain_right_edge.copy()
        ma = pf.domain_left_edge.copy()
        for off1 in [-1, 1]:
            for off2 in [-1, 1]:
                for off3 in [-1, 1]:
                    this_point = (self.center + width/2. * off1 * north_vector
                                         + width/2. * off2 * east_vector
                                         + width/2. * off3 * normal_vector)
                    na.minimum(mi, this_point, mi)
                    na.maximum(ma, this_point, ma)
        # Now we have a bounding box.
        grids = pf.h.region(self.center, mi, ma)._grids

        pb = get_pbar("Sampling ", len(grids))
        for i,grid in enumerate(grids):
            data = [(grid[field] * grid.child_mask).astype("float64")
                    for field in fields]
            pg = PartitionedGrid(
                grid.id, data,
                grid.LeftEdge, grid.RightEdge, grid.ActiveDimensions.astype("int64"))
            grid.clear_data()
            sampler(pg, num_threads = num_threads)
            pb.update(i)
        pb.finish()

        image = sampler.aimage
        self.finalize_image(image)
        return image

    def save_image(self, fn, clip_ratio, image):
        if self.pf.field_info[self.field].take_log:
            im = na.log10(image)
        else:
            im = image
        if self.comm.rank is 0 and fn is not None:
            if clip_ratio is not None:
                write_image(im, fn)
            else:
                write_image(im, fn)

    def snapshot(self, fn = None, clip_ratio = None, double_check = False,
                 num_threads = 0):

        if num_threads is None:
            num_threads=get_num_threads()

        fields = [self.field]
        resolution = self.resolution

        image = self.new_image()

        args = self.get_sampler_args(image)

        sampler = self.get_sampler(args)

        self.initialize_source()

        image = self._render(double_check, num_threads, image, sampler)

        self.save_image(fn, clip_ratio, image)

        return image
    snapshot.__doc__ = Camera.snapshot.__doc__

data_object_registry["projection_camera"] = ProjectionCamera

def off_axis_projection(pf, center, normal_vector, width, resolution,
                        field, weight = None, 
                        volume = None, no_ghost = False, interpolated = False):
    r"""Project through a parameter file, off-axis, and return the image plane.

    This function will accept the necessary items to integrate through a volume
    at an arbitrary angle and return the integrated field of view to the user.
    Note that if a weight is supplied, it will multiply the pre-interpolated
    values together, then create cell-centered values, then interpolate within
    the cell to conduct the integration.

    Parameters
    ----------
    pf : `~yt.data_objects.api.StaticOutput`
        This is the parameter file to volume render.
    center : array_like
        The current 'center' of the view port -- the focal point for the
        camera.
    normal_vector : array_like
        The vector between the camera position and the center.
    width : float or list of floats
        The current width of the image.  If a single float, the volume is
        cubical, but if not, it is left/right, top/bottom, front/back
    resolution : int or list of ints
        The number of pixels in each direction.
    field : string
        The field to project through the volume
    weight : optional, default None
        If supplied, the field will be pre-multiplied by this, then divided by
        the integrated value of this field.  This returns an average rather
        than a sum.
    volume : `yt.extensions.volume_rendering.HomogenizedVolume`, optional
        The volume to ray cast through.  Can be specified for finer-grained
        control, but otherwise will be automatically generated.
    no_ghost: bool, optional
        Optimization option.  If True, homogenized bricks will
        extrapolate out from grid instead of interpolating from
        ghost zones that have to first be calculated.  This can
        lead to large speed improvements, but at a loss of
        accuracy/smoothness in resulting image.  The effects are
        less notable when the transfer function is smooth and
        broad. Default: True
    interpolated : optional, default False
        If True, the data is first interpolated to vertex-centered data, 
        then tri-linearly interpolated along the ray. Not suggested for 
        quantitative studies.

    Returns
    -------
    image : array
        An (N,N) array of the final integrated values, in float64 form.

    Examples
    --------

    >>> image = off_axis_projection(pf, [0.5, 0.5, 0.5], [0.2,0.3,0.4],
                      0.2, N, "Temperature", "Density")
    >>> write_image(na.log10(image), "offaxis.png")

    """
    projcam = ProjectionCamera(center, normal_vector, width, resolution,
            field, weight=weight, pf=pf, volume=volume,
            no_ghost=no_ghost, interpolated=interpolated)
    image = projcam.snapshot()
    if weight is not None:
        pf.field_info.pop("temp_weightfield")
    del projcam
    return image[:,:,0]

