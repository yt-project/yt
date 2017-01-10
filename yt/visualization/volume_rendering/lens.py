"""
Lens Classes



"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import division
from yt.funcs import mylog
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface
from yt.data_objects.image_array import ImageArray
from yt.units.yt_array import unorm, uvstack, uhstack
from yt.utilities.math_utils import get_rotation_matrix
import numpy as np

from yt.utilities.lib.grid_traversal import \
    arr_fisheye_vectors


class Lens(ParallelAnalysisInterface):
    """A Lens is used to define the set of rays for rendering."""
    def __init__(self, ):
        super(Lens, self).__init__()
        self.viewpoint = None
        self.sub_samples = 5
        self.num_threads = 0
        self.box_vectors = None
        self.origin = None
        self.back_center = None
        self.front_center = None
        self.sampler = None

    def set_camera(self, camera):
        """Set the properties of the lens based on the camera.

        This is a proxy for setup_box_properties
        """
        self.setup_box_properties(camera)

    def new_image(self, camera):
        """Initialize a new ImageArray to be used with this lens."""
        self.current_image = ImageArray(
            np.zeros((camera.resolution[0], camera.resolution[1],
                      4), dtype='float64', order='C'),
            info={'imtype': 'rendering'})
        return self.current_image

    def setup_box_properties(self, camera):
        """Set up the view and stage based on the properties of the camera."""
        unit_vectors = camera.unit_vectors
        width = camera.width
        center = camera.focus

        self.box_vectors = camera.scene.arr(
            [unit_vectors[0] * width[0],
             unit_vectors[1] * width[1],
             unit_vectors[2] * width[2]])
        self.origin = center - 0.5 * width.dot(unit_vectors)
        self.back_center = center - 0.5 * width[2] * unit_vectors[2]
        self.front_center = center + 0.5 * width[2] * unit_vectors[2]
        self.set_viewpoint(camera)

    def set_viewpoint(self, camera):
        """
        Set the viewpoint used for AMRKDTree traversal such that you yield
        bricks from back to front or front to back from with respect to this
        point.  Must be implemented for each Lens type.
        """
        raise NotImplementedError("Need to choose viewpoint for this class")


class PlaneParallelLens(Lens):
    r"""The lens for orthographic projections.

    All rays emerge parallel to each other, arranged along a plane.

    The initializer takes no parameters.
    """

    def __init__(self, ):
        super(PlaneParallelLens, self).__init__()

    def _get_sampler_params(self, camera, render_source):
        if render_source.zbuffer is not None:
            image = render_source.zbuffer.rgba
        else:
            image = self.new_image(camera)

        vp_pos = np.concatenate(
            [camera.inv_mat.ravel('F').d,
             self.back_center.ravel().in_units('code_length').d])

        sampler_params =\
            dict(vp_pos=vp_pos,
                 vp_dir=self.box_vectors[2],  # All the same
                 center=self.back_center,
                 bounds=(-camera.width[0] / 2.0, camera.width[0] / 2.0,
                         -camera.width[1] / 2.0, camera.width[1] / 2.0),
                 x_vec=camera.unit_vectors[0],
                 y_vec=camera.unit_vectors[1],
                 width=np.array(camera.width, dtype='float64'),
                 image=image, lens_type="plane-parallel")
        return sampler_params

    def set_viewpoint(self, camera):
        """Set the viewpoint based on the camera"""
        # This is a hack that should be replaced by an alternate plane-parallel
        # traversal. Put the camera really far away so that the effective
        # viewpoint is infinitely far away, making for parallel rays.
        self.viewpoint = self.front_center + \
            camera.unit_vectors[2] * 1.0e6 * camera.width[2]

    def project_to_plane(self, camera, pos, res=None):
        if res is None:
            res = camera.resolution

        origin = self.origin.in_units('code_length').d
        front_center = self.front_center.in_units('code_length').d
        width = camera.width.in_units('code_length').d

        dx = np.array(np.dot(pos - origin, camera.unit_vectors[1]))
        dy = np.array(np.dot(pos - origin, camera.unit_vectors[0]))
        dz = np.array(np.dot(pos - front_center, -camera.unit_vectors[2]))
        # Transpose into image coords.

        py = (res[0]*(dx/width[0])).astype('int64')
        px = (res[1]*(dy/width[1])).astype('int64')
        return px, py, dz

    def __repr__(self):
        disp = ("<Lens Object>:\n\tlens_type:plane-parallel\n\tviewpoint:%s" %
                (self.viewpoint))
        return disp


class PerspectiveLens(Lens):
    r"""A lens for viewing a scene with a set of rays within an opening angle.

    The scene will have an element of perspective to it since the rays are not
    parallel.
    """

    def __init__(self):
        super(PerspectiveLens, self).__init__()

    def new_image(self, camera):
        self.current_image = ImageArray(
            np.zeros((camera.resolution[0], camera.resolution[1], 
                      4), dtype='float64', order='C'),
            info={'imtype': 'rendering'})
        return self.current_image

    def _get_sampler_params(self, camera, render_source):
        # Enforce width[1] / width[0] = resolution[1] / resolution[0]
        camera.width[1] = camera.width[0] \
            * (camera.resolution[1] / camera.resolution[0])

        if render_source.zbuffer is not None:
            image = render_source.zbuffer.rgba
        else:
            image = self.new_image(camera)

        east_vec = camera.unit_vectors[0]
        north_vec = camera.unit_vectors[1]
        normal_vec = camera.unit_vectors[2]

        px = np.mat(np.linspace(-.5, .5, camera.resolution[0]))
        py = np.mat(np.linspace(-.5, .5, camera.resolution[1]))

        sample_x = camera.width[0] * np.array(east_vec.reshape(3, 1) * px)
        sample_x = sample_x.transpose()
        sample_y = camera.width[1] * np.array(north_vec.reshape(3, 1) * py)
        sample_y = sample_y.transpose()

        vectors = np.zeros((camera.resolution[0], camera.resolution[1], 3),
                           dtype='float64', order='C')

        sample_x = np.repeat(sample_x.reshape(camera.resolution[0], 1, 3),
                             camera.resolution[1], axis=1)
        sample_y = np.repeat(sample_y.reshape(1, camera.resolution[1], 3),
                             camera.resolution[0], axis=0)

        normal_vecs = np.tile(
            normal_vec, camera.resolution[0] * camera.resolution[1])
        normal_vecs = normal_vecs.reshape(
            camera.resolution[0], camera.resolution[1], 3)

        # The maximum possible length of ray
        max_length = (unorm(camera.position - camera._domain_center)
                      + 0.5 * unorm(camera._domain_width))

        # Rescale the ray to be long enough to cover the entire domain
        vectors = (sample_x + sample_y + normal_vecs * camera.width[2]) * \
            (max_length / camera.width[2])

        positions = np.tile(
            camera.position, camera.resolution[0] * camera.resolution[1])
        positions = positions.reshape(
            camera.resolution[0], camera.resolution[1], 3)

        uv = np.ones(3, dtype='float64')

        image = self.new_image(camera)

        sampler_params =\
            dict(vp_pos=positions,
                 vp_dir=vectors,
                 center=self.back_center,
                 bounds=(0.0, 1.0, 0.0, 1.0),
                 x_vec=uv,
                 y_vec=uv,
                 width=np.zeros(3, dtype='float64'),
                 image=image,
                 lens_type="perspective")

        mylog.debug(positions)
        mylog.debug(vectors)

        return sampler_params

    def set_viewpoint(self, camera):
        """
        For a PerspectiveLens, the viewpoint is the front center.
        """
        self.viewpoint = self.front_center

    def project_to_plane(self, camera, pos, res=None):
        if res is None:
            res = camera.resolution

        width = camera.width.in_units('code_length').d
        position = camera.position.in_units('code_length').d

        width[1] = width[0] * res[1] / res[0]

        sight_vector = pos - position

        pos1 = sight_vector

        for i in range(0, sight_vector.shape[0]):
            sight_vector_norm = np.sqrt(np.dot(sight_vector[i], sight_vector[i]))
            if sight_vector_norm != 0:
                sight_vector[i] = sight_vector[i] / sight_vector_norm

        sight_center = camera.position + camera.width[2] * camera.unit_vectors[2]

        sight_center = sight_center.in_units('code_length').d

        for i in range(0, sight_vector.shape[0]):
            sight_angle_cos = np.dot(sight_vector[i], camera.unit_vectors[2])
            # clip sight_angle_cos since floating point noise might
            # go outside the domain of arccos
            sight_angle_cos = np.clip(sight_angle_cos, -1.0, 1.0)
            if np.arccos(sight_angle_cos) < 0.5 * np.pi:
                sight_length = width[2] / sight_angle_cos
            else:
                # If the corner is on the backwards, then we put it outside of
                # the image It can not be simply removed because it may connect
                # to other corner within the image, which produces visible
                # domain boundary line
                sight_length = np.sqrt(width[0]**2 + width[1]**2)
                sight_length = sight_length / np.sqrt(1 - sight_angle_cos**2)
            pos1[i] = position + sight_length * sight_vector[i]

        dx = np.dot(pos1 - sight_center, camera.unit_vectors[0])
        dy = np.dot(pos1 - sight_center, camera.unit_vectors[1])
        dz = np.dot(pos - position, camera.unit_vectors[2])

        # Transpose into image coords.
        px = (res[0] * 0.5 + res[0] / width[0] * dx).astype('int64')
        py = (res[1] * 0.5 + res[1] / width[1] * dy).astype('int64')

        return px, py, dz

    def __repr__(self):
        disp = ("<Lens Object>:\n\tlens_type:perspective\n\tviewpoint:%s" %
                (self.viewpoint))
        return disp


class StereoPerspectiveLens(Lens):
    """A lens that includes two sources for perspective rays, for 3D viewing

    The degree of differences between the left and right images is controlled by 
    the disparity (the maximum distance between cooresponding points in the left
    and right images). By default, the disparity is set to be 3 pixels.
    """

    def __init__(self):
        super(StereoPerspectiveLens, self).__init__()
        self.disparity = None

    def new_image(self, camera):
        """Initialize a new ImageArray to be used with this lens."""
        self.current_image = ImageArray(
            np.zeros((camera.resolution[0], camera.resolution[1], 4),
                     dtype='float64', order='C'),
            info={'imtype': 'rendering'})
        return self.current_image

    def _get_sampler_params(self, camera, render_source):
        # Enforce width[1] / width[0] = 2 * resolution[1] / resolution[0]
        # For stereo-type lens, images for left and right eye are pasted together,
        # so the resolution of single-eye image will be 50% of the whole one.
        camera.width[1] = camera.width[0] \
            * (2. * camera.resolution[1] / camera.resolution[0])

        if self.disparity is None:
            self.disparity = 3. * camera.width[0] / camera.resolution[0]

        if render_source.zbuffer is not None:
            image = render_source.zbuffer.rgba
        else:
            image = self.new_image(camera)

        vectors_left, positions_left = self._get_positions_vectors(
            camera, -self.disparity)
        vectors_right, positions_right = self._get_positions_vectors(
            camera, self.disparity)

        uv = np.ones(3, dtype='float64')

        image = self.new_image(camera)
        vectors_comb = uvstack([vectors_left, vectors_right])
        positions_comb = uvstack([positions_left, positions_right])

        image.shape = (camera.resolution[0], camera.resolution[1], 4)
        vectors_comb.shape = (camera.resolution[0], camera.resolution[1], 3)
        positions_comb.shape = (camera.resolution[0], camera.resolution[1], 3)

        sampler_params =\
            dict(vp_pos=positions_comb,
                 vp_dir=vectors_comb,
                 center=self.back_center,
                 bounds=(0.0, 1.0, 0.0, 1.0),
                 x_vec=uv,
                 y_vec=uv,
                 width=np.zeros(3, dtype='float64'),
                 image=image,
                 lens_type="stereo-perspective")

        return sampler_params

    def _get_positions_vectors(self, camera, disparity):

        single_resolution_x = int(np.floor(camera.resolution[0]) / 2)

        east_vec = camera.unit_vectors[0]
        north_vec = camera.unit_vectors[1]
        normal_vec = camera.unit_vectors[2]

        angle_disparity = - np.arctan2(disparity.in_units(camera.width.units),
                                       camera.width[2])
        R = get_rotation_matrix(angle_disparity, north_vec)

        east_vec_rot = np.dot(R, east_vec)
        normal_vec_rot = np.dot(R, normal_vec)

        px = np.mat(np.linspace(-.5, .5, single_resolution_x))
        py = np.mat(np.linspace(-.5, .5, camera.resolution[1]))

        sample_x = camera.width[0] * np.array(east_vec_rot.reshape(3, 1) * px)
        sample_x = sample_x.transpose()
        sample_y = camera.width[1] * np.array(north_vec.reshape(3, 1) * py)
        sample_y = sample_y.transpose()

        vectors = np.zeros((single_resolution_x, camera.resolution[1], 3),
                           dtype='float64', order='C')

        sample_x = np.repeat(sample_x.reshape(single_resolution_x, 1, 3),
                             camera.resolution[1], axis=1)
        sample_y = np.repeat(sample_y.reshape(1, camera.resolution[1], 3),
                             single_resolution_x, axis=0)

        normal_vecs = np.tile(
            normal_vec_rot, single_resolution_x * camera.resolution[1])
        normal_vecs = normal_vecs.reshape(
            single_resolution_x, camera.resolution[1], 3)
        east_vecs = np.tile(
            east_vec_rot, single_resolution_x * camera.resolution[1])
        east_vecs = east_vecs.reshape(
            single_resolution_x, camera.resolution[1], 3)

        # The maximum possible length of ray
        max_length = (unorm(camera.position - camera._domain_center)
                      + 0.5 * unorm(camera._domain_width)
                      + np.abs(self.disparity))
        # Rescale the ray to be long enough to cover the entire domain
        vectors = (sample_x + sample_y + normal_vecs * camera.width[2]) * \
            (max_length / camera.width[2])

        positions = np.tile(
            camera.position, single_resolution_x * camera.resolution[1])
        positions = positions.reshape(
            single_resolution_x, camera.resolution[1], 3)

        # Here the east_vecs is non-rotated one
        positions = positions + east_vecs * disparity

        mylog.debug(positions)
        mylog.debug(vectors)

        return vectors, positions

    def project_to_plane(self, camera, pos, res=None):
        if res is None:
            res = camera.resolution

        # Enforce width[1] / width[0] = 2 * resolution[1] / resolution[0]
        # For stereo-type lens, images for left and right eye are pasted together,
        # so the resolution of single-eye image will be 50% of the whole one.
        camera.width[1] = camera.width[0] * (2. * res[1] / res[0])

        if self.disparity is None:
            self.disparity = 3. * camera.width[0] / camera.resolution[0]

        px_left, py_left, dz_left = self._get_px_py_dz(
            camera, pos, res, -self.disparity)
        px_right, py_right, dz_right = self._get_px_py_dz(
            camera, pos, res, self.disparity)

        px = uvstack([px_left, px_right])
        py = uvstack([py_left, py_right])
        dz = uvstack([dz_left, dz_right])

        return px, py, dz

    def _get_px_py_dz(self, camera, pos, res, disparity):

        res0_h = np.floor(res[0]) / 2

        east_vec = camera.unit_vectors[0]
        north_vec = camera.unit_vectors[1]
        normal_vec = camera.unit_vectors[2]

        angle_disparity = - np.arctan2(disparity.d, camera.width[2].d)
        R = get_rotation_matrix(angle_disparity, north_vec)

        east_vec_rot = np.dot(R, east_vec)
        normal_vec_rot = np.dot(R, normal_vec)

        camera_position_shift = camera.position + east_vec * disparity
        camera_position_shift = camera_position_shift.in_units('code_length').d
        width = camera.width.in_units('code_length').d
        sight_vector = pos - camera_position_shift
        pos1 = sight_vector

        for i in range(0, sight_vector.shape[0]):
            sight_vector_norm = np.sqrt(np.dot(sight_vector[i], sight_vector[i]))
            sight_vector[i] = sight_vector[i] / sight_vector_norm
        sight_center = camera_position_shift + camera.width[2] * normal_vec_rot

        for i in range(0, sight_vector.shape[0]):
            sight_angle_cos = np.dot(sight_vector[i], normal_vec_rot)
            # clip sight_angle_cos since floating point noise might
            # cause it go outside the domain of arccos
            sight_angle_cos = np.clip(sight_angle_cos, -1.0, 1.0)
            if np.arccos(sight_angle_cos) < 0.5 * np.pi:
                sight_length = width[2] / sight_angle_cos
            else:
                # If the corner is on the backwards, then we put it outside of
                # the image It can not be simply removed because it may connect
                # to other corner within the image, which produces visible
                # domain boundary line
                sight_length = np.sqrt(width[0]**2 + width[1]**2)
                sight_length = sight_length / np.sqrt(1 - sight_angle_cos**2)
            pos1[i] = camera_position_shift + sight_length * sight_vector[i]

        dx = np.dot(pos1 - sight_center, east_vec_rot)
        dy = np.dot(pos1 - sight_center, north_vec)
        dz = np.dot(pos - camera_position_shift, normal_vec_rot)

        # Transpose into image coords.
        if disparity > 0:
            px = (res0_h * 0.5 + res0_h / camera.width[0].d * dx).astype('int')
        else:
            px = (res0_h * 1.5 + res0_h / camera.width[0].d * dx).astype('int')
        py = (res[1] * 0.5 + res[1] / camera.width[1].d * dy).astype('int')

        return px, py, dz

    def set_viewpoint(self, camera):
        """
        For a PerspectiveLens, the viewpoint is the front center.
        """
        self.viewpoint = self.front_center

    def __repr__(self):
        disp = ("<Lens Object>:\n\tlens_type:perspective\n\tviewpoint:%s" %
                (self.viewpoint))
        return disp


class FisheyeLens(Lens):
    r"""A lens for dome-based renderings

    This lens type accepts a field-of-view property, fov, that describes how
    wide an angle the fisheye can see. Fisheye images are typically used for
    dome-based presentations; the Hayden planetarium for instance has a field of
    view of 194.6.  The images returned by this camera will be flat pixel images
    that can and should be reshaped to the resolution.

    """

    def __init__(self):
        super(FisheyeLens, self).__init__()
        self.fov = 180.0
        self.radius = 1.0
        self.center = None
        self.rotation_matrix = np.eye(3)

    def setup_box_properties(self, camera):
        """Set up the view and stage based on the properties of the camera."""
        self.radius = camera.width.max()
        super(FisheyeLens, self).setup_box_properties(camera)
        self.set_viewpoint(camera)

    def new_image(self, camera):
        """Initialize a new ImageArray to be used with this lens."""
        self.current_image = ImageArray(
            np.zeros((camera.resolution[0], camera.resolution[0],
                      4), dtype='float64', order='C'),
            info={'imtype': 'rendering'})
        return self.current_image

    def _get_sampler_params(self, camera, render_source):
        vp = -arr_fisheye_vectors(camera.resolution[0], self.fov)
        vp.shape = (camera.resolution[0], camera.resolution[0], 3)
        vp = vp.dot(np.linalg.inv(self.rotation_matrix))
        vp *= self.radius
        uv = np.ones(3, dtype='float64')
        positions = np.ones((camera.resolution[0], camera.resolution[0], 3),
            dtype='float64') * camera.position

        if render_source.zbuffer is not None:
            image = render_source.zbuffer.rgba
        else:
            image = self.new_image(camera)

        sampler_params =\
            dict(vp_pos=positions,
                 vp_dir=vp,
                 center=self.center,
                 bounds=(0.0, 1.0, 0.0, 1.0),
                 x_vec=uv,
                 y_vec=uv,
                 width=np.zeros(3, dtype='float64'),
                 image=image,
                 lens_type="fisheye")

        return sampler_params

    def set_viewpoint(self, camera):
        """For a FisheyeLens, the viewpoint is the camera's position"""
        self.viewpoint = camera.position

    def __repr__(self):
        disp = ("<Lens Object>:\n\tlens_type:fisheye\n\tviewpoint:%s"
                "\nt\tfov:%s\n\tradius:%s" %
                (self.viewpoint, self.fov, self.radius))
        return disp

    def project_to_plane(self, camera, pos, res=None):
        if res is None:
            res = camera.resolution
        # the return values here need to be px, py, dz
        # these are the coordinates and dz for the resultant image.
        # Basically, what we need is an inverse projection from the fisheye
        # vectors back onto the plane.  arr_fisheye_vectors goes from px, py to
        # vector, and we need the reverse.
        # First, we transform lpos into *relative to the camera* coordinates.

        position = camera.position.in_units('code_length').d

        lpos = position - pos
        lpos = lpos.dot(self.rotation_matrix)
        mag = (lpos * lpos).sum(axis=1)**0.5

        # screen out NaN values that would result from dividing by mag
        mag[mag == 0] = 1
        lpos /= mag[:, None]

        dz = (mag / self.radius).in_units('1/code_length').d
        theta = np.arccos(lpos[:, 2])
        fov_rad = self.fov * np.pi / 180.0
        r = 2.0 * theta / fov_rad
        phi = np.arctan2(lpos[:, 1], lpos[:, 0])
        px = r * np.cos(phi)
        py = r * np.sin(phi)

        # dz is distance the ray would travel
        px = (px + 1.0) * res[0] / 2.0
        py = (py + 1.0) * res[1] / 2.0
        # px and py should be dimensionless
        px = np.rint(px).astype("int64")
        py = np.rint(py).astype("int64")
        return px, py, dz


class SphericalLens(Lens):
    r"""A lens for cylindrical-spherical projection.

    Movies rendered in this way can be displayed in head-tracking devices or
    in YouTube 360 view.

    """

    def __init__(self):
        super(SphericalLens, self).__init__()
        self.radius = 1.0
        self.center = None
        self.rotation_matrix = np.eye(3)

    def setup_box_properties(self, camera):
        """Set up the view and stage based on the properties of the camera."""
        self.radius = camera.width.max()
        super(SphericalLens, self).setup_box_properties(camera)
        self.set_viewpoint(camera)

    def _get_sampler_params(self, camera, render_source):
        px = np.linspace(-np.pi, np.pi, camera.resolution[0],
                         endpoint=True)[:, None]
        py = np.linspace(-np.pi/2., np.pi/2., camera.resolution[1],
                         endpoint=True)[None, :]

        vectors = np.zeros((camera.resolution[0], camera.resolution[1], 3),
                           dtype='float64', order='C')
        vectors[:, :, 0] = np.cos(px) * np.cos(py)
        vectors[:, :, 1] = np.sin(px) * np.cos(py)
        vectors[:, :, 2] = np.sin(py)

        # The maximum possible length of ray
        max_length = (unorm(camera.position - camera._domain_center)
                      + 0.5 * unorm(camera._domain_width))
        # Rescale the ray to be long enough to cover the entire domain
        vectors = vectors * max_length

        positions = np.tile(
            camera.position,
            camera.resolution[0] * camera.resolution[1]).reshape(
                camera.resolution[0], camera.resolution[1], 3)

        R1 = get_rotation_matrix(0.5*np.pi, [1, 0, 0])
        R2 = get_rotation_matrix(0.5*np.pi, [0, 0, 1])
        uv = np.dot(R1, camera.unit_vectors)
        uv = np.dot(R2, uv)
        vectors.reshape((camera.resolution[0]*camera.resolution[1], 3))
        vectors = np.dot(vectors, uv)
        vectors.reshape((camera.resolution[0], camera.resolution[1], 3))

        if render_source.zbuffer is not None:
            image = render_source.zbuffer.rgba
        else:
            image = self.new_image(camera)

        dummy = np.ones(3, dtype='float64')
        image.shape = (camera.resolution[0], camera.resolution[1], 4)
        vectors.shape = (camera.resolution[0], camera.resolution[1], 3)
        positions.shape = (camera.resolution[0], camera.resolution[1], 3)

        sampler_params = dict(
            vp_pos=positions,
            vp_dir=vectors,
            center=self.back_center,
            bounds=(0.0, 1.0, 0.0, 1.0),
            x_vec=dummy,
            y_vec=dummy,
            width=np.zeros(3, dtype="float64"),
            image=image,
            lens_type="spherical")
        return sampler_params

    def set_viewpoint(self, camera):
        """For a SphericalLens, the viewpoint is the camera's position"""
        self.viewpoint = camera.position

    def project_to_plane(self, camera, pos, res=None):
        if res is None:
            res = camera.resolution
        # Much of our setup here is the same as in the fisheye, except for the
        # actual conversion back to the px, py values.
        position = camera.position.in_units('code_length').d

        lpos = position - pos
        mag = (lpos * lpos).sum(axis=1)**0.5

        # screen out NaN values that would result from dividing by mag
        mag[mag == 0] = 1
        lpos /= mag[:, None]

        # originally:
        #  the x vector is cos(px) * cos(py)
        #  the y vector is sin(px) * cos(py)
        #  the z vector is sin(py)
        # y / x = tan(px), so arctan2(lpos[:,1], lpos[:,0]) => px
        # z = sin(py) so arcsin(z) = py
        # px runs from -pi to pi
        # py runs from -pi/2 to pi/2
        px = np.arctan2(lpos[:, 1], lpos[:, 0])
        py = np.arcsin(lpos[:, 2])
        dz = mag / self.radius
        # dz is distance the ray would travel
        px = ((-px + np.pi) / (2.0*np.pi)) * res[0]
        py = ((-py + np.pi/2.0) / np.pi) * res[1]
        # px and py should be dimensionless
        px = np.rint(px).astype("int64")
        py = np.rint(py).astype("int64")
        return px, py, dz


class StereoSphericalLens(Lens):
    r"""A lens for a stereo cylindrical-spherical projection.

    Movies rendered in this way can be displayed in VR devices or stereo youtube
    360 degree movies.
    """

    def __init__(self):
        super(StereoSphericalLens, self).__init__()
        self.radius = 1.0
        self.center = None
        self.disparity = None
        self.rotation_matrix = np.eye(3)

    def setup_box_properties(self, camera):
        self.radius = camera.width.max()
        super(StereoSphericalLens, self).setup_box_properties(camera)
        self.set_viewpoint(camera)

    def _get_sampler_params(self, camera, render_source):
        if self.disparity is None:
            self.disparity = camera.width[0] / 1000.

        single_resolution_y = int(np.floor(camera.resolution[1]) / 2)
        px = np.linspace(-np.pi, np.pi, camera.resolution[0],
                         endpoint=True)[:, None]
        py = np.linspace(-np.pi/2., np.pi/2., single_resolution_y,
                         endpoint=True)[None, :]

        vectors = np.zeros((camera.resolution[0], single_resolution_y, 3),
                           dtype='float64', order='C')
        vectors[:, :, 0] = np.cos(px) * np.cos(py)
        vectors[:, :, 1] = np.sin(px) * np.cos(py)
        vectors[:, :, 2] = np.sin(py)

        # The maximum possible length of ray
        max_length = (unorm(camera.position - camera._domain_center)
                      + 0.5 * unorm(camera._domain_width)
                      + np.abs(self.disparity))
        # Rescale the ray to be long enough to cover the entire domain
        vectors = vectors * max_length

        R1 = get_rotation_matrix(0.5*np.pi, [1, 0, 0])
        R2 = get_rotation_matrix(0.5*np.pi, [0, 0, 1])
        uv = np.dot(R1, camera.unit_vectors)
        uv = np.dot(R2, uv)

        vectors.reshape((camera.resolution[0]*single_resolution_y, 3))
        vectors = np.dot(vectors, uv)
        vectors.reshape((camera.resolution[0], single_resolution_y, 3))

        vectors2 = np.zeros((camera.resolution[0], single_resolution_y, 3),
                            dtype='float64', order='C')
        vectors2[:, :, 0] = -np.sin(px) * np.ones((1, single_resolution_y))
        vectors2[:, :, 1] = np.cos(px) * np.ones((1, single_resolution_y))
        vectors2[:, :, 2] = 0

        vectors2.reshape((camera.resolution[0]*single_resolution_y, 3))
        vectors2 = np.dot(vectors2, uv)
        vectors2.reshape((camera.resolution[0], single_resolution_y, 3))

        positions = np.tile(
            camera.position, camera.resolution[0] * single_resolution_y)
        positions = positions.reshape(
            camera.resolution[0], single_resolution_y, 3)

        # The left and right are switched here since VR is in LHS.
        positions_left = positions + vectors2 * self.disparity
        positions_right = positions + vectors2 * (-self.disparity)

        if render_source.zbuffer is not None:
            image = render_source.zbuffer.rgba
        else:
            image = self.new_image(camera)

        dummy = np.ones(3, dtype='float64')

        vectors_comb = uhstack([vectors, vectors])
        positions_comb = uhstack([positions_left, positions_right])

        image.shape = (camera.resolution[0], camera.resolution[1], 4)
        vectors_comb.shape = (camera.resolution[0], camera.resolution[1], 3)
        positions_comb.shape = (camera.resolution[0], camera.resolution[1], 3)

        sampler_params = dict(
            vp_pos=positions_comb,
            vp_dir=vectors_comb,
            center=self.back_center,
            bounds=(0.0, 1.0, 0.0, 1.0),
            x_vec=dummy,
            y_vec=dummy,
            width=np.zeros(3, dtype="float64"),
            image=image,
            lens_type = "stereo-spherical")
        return sampler_params

    def set_viewpoint(self, camera):
        """
        For a PerspectiveLens, the viewpoint is the front center.
        """
        self.viewpoint = camera.position

lenses = {'plane-parallel': PlaneParallelLens,
          'perspective': PerspectiveLens,
          'stereo-perspective': StereoPerspectiveLens,
          'fisheye': FisheyeLens,
          'spherical': SphericalLens,
          'stereo-spherical': StereoSphericalLens}
