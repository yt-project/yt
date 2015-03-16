"""
Lens Classes



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.funcs import mylog
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface
from yt.units.yt_array import YTArray
from yt.data_objects.image_array import ImageArray
from yt.utilities.math_utils import get_rotation_matrix
import numpy as np

from yt.utilities.lib.grid_traversal import \
    arr_fisheye_vectors


class Lens(ParallelAnalysisInterface):

    """docstring for Lens"""

    def __init__(self, ):
        #mylog.debug("Entering %s" % str(self))
        super(Lens, self).__init__()
        self.viewpoint = None
        self.sub_samples = 5
        self.num_threads = 0
        self.double_check = False
        self.box_vectors = None
        self.origin = None
        self.back_center = None
        self.front_center = None
        self.sampler = None

    def set_camera(self, camera):
        self.setup_box_properties(camera)

    def new_image(self, camera):
        self.current_image = ImageArray(
            np.zeros((camera.resolution[0], camera.resolution[1],
                      4), dtype='float64', order='C'),
            info={'imtype': 'rendering'})
        return self.current_image

    def setup_box_properties(self, camera):
        unit_vectors = camera.unit_vectors
        width = camera.width
        center = camera.focus
        self.box_vectors = YTArray([unit_vectors[0] * width[0],
                                    unit_vectors[1] * width[1],
                                    unit_vectors[2] * width[2]])
        self.origin = center - 0.5 * width.dot(YTArray(unit_vectors, ""))
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

    """docstring for PlaneParallelLens"""

    def __init__(self, ):
        #mylog.debug("Entering %s" % str(self))
        super(PlaneParallelLens, self).__init__()

    def _get_sampler_params(self, camera, render_source):
        if render_source.zbuffer is not None:
            image = render_source.zbuffer.rgba
        else:
            image = self.new_image(camera)

        sampler_params =\
            dict(vp_pos=np.concatenate([camera.inv_mat.ravel('F'),
                                        self.back_center.ravel()]),
                 vp_dir=self.box_vectors[2],  # All the same
                 center=self.back_center,
                 bounds=(-camera.width[0] / 2.0, camera.width[0] / 2.0,
                         -camera.width[1] / 2.0, camera.width[1] / 2.0),
                 x_vec=camera.unit_vectors[0],
                 y_vec=camera.unit_vectors[1],
                 width=np.array(camera.width, dtype='float64'),
                 image=image)
        return sampler_params

    def set_viewpoint(self, camera):
        # This is a hack that should be replaced by an alternate plane-parallel
        # traversal. Put the camera really far away so that the effective
        # viewpoint is infinitely far away, making for parallel rays.
        self.viewpoint = self.front_center + \
            camera.unit_vectors[2] * 1.0e6 * camera.width[2]

    def project_to_plane(self, camera, pos, res=None):
        if res is None:
            res = camera.resolution
        dx = np.dot(pos - self.origin.d, camera.unit_vectors[1])
        dy = np.dot(pos - self.origin.d, camera.unit_vectors[0])
        dz = np.dot(pos - self.front_center.d, -camera.unit_vectors[2])
        # Transpose into image coords.
        py = (res[0]*(dx/camera.width[0].d)).astype('int')
        px = (res[1]*(dy/camera.width[1].d)).astype('int')
        return px, py, dz

    def __repr__(self):
        disp = "<Lens Object>:\n\tlens_type:plane-parallel\n\tviewpoint:%s" %\
                (self.viewpoint)
        return disp

class PerspectiveLens(Lens):

    """docstring for PerspectiveLens"""

    def __init__(self):
        super(PerspectiveLens, self).__init__()
        self.expand_factor = 1.5

    def new_image(self, camera):
        self.current_image = ImageArray(
            np.zeros((camera.resolution[0]*camera.resolution[1], 1,
                      4), dtype='float64', order='C'),
            info={'imtype': 'rendering'})
        return self.current_image

    def _get_sampler_params(self, camera, render_source):
        # We should move away from pre-generation of vectors like this and into
        # the usage of on-the-fly generation in the VolumeIntegrator module
        # We might have a different width and back_center
        #dl = (self.back_center - self.front_center)
        #self.front_center += self.expand_factor*dl
        #self.back_center -= dl

        px = np.linspace(-camera.width[0]/2.0, camera.width[0]/2.0,
                         camera.resolution[0])[:, None].d
        py = np.linspace(-camera.width[1]/2.0, camera.width[1]/2.0,
                         camera.resolution[1])[None, :].d
        inv_mat = camera.inv_mat.d
        positions = np.zeros((camera.resolution[0], camera.resolution[1], 3),
                             dtype='float64', order='C')
        positions[:, :, 0] = inv_mat[0, 0]*px + \
            inv_mat[0, 1]*py + self.back_center.d[0]
        positions[:, :, 1] = inv_mat[1, 0]*px + \
            inv_mat[1, 1]*py + self.back_center.d[1]
        positions[:, :, 2] = inv_mat[2, 0]*px + \
            inv_mat[2, 1]*py + self.back_center.d[2]
        # Can we use bounds for anything here?
        # bounds = (px.min(), px.max(), py.min(), py.max())

        # We are likely adding on an odd cutting condition here
        vectors = self.front_center.d - positions
        vectors = vectors / (vectors**2).sum()**0.5
        positions = self.front_center.d - 1.0 * \
            (((self.back_center.d-self.front_center.d)**2).sum())**0.5*vectors
        vectors = (self.back_center.d - positions)

        uv = np.ones(3, dtype='float64')
        vectors.shape = (camera.resolution[0]**2, 1, 3)
        positions.shape = (camera.resolution[0]**2, 1, 3)
        image = self.new_image(camera)

        sampler_params =\
            dict(vp_pos=positions,
                 vp_dir=vectors,
                 center=self.back_center.d,
                 bounds=(0.0, 1.0, 0.0, 1.0),
                 x_vec=uv,
                 y_vec=uv,
                 width=np.zeros(3, dtype='float64'),
                 image=image
                 )

        mylog.debug(positions)
        mylog.debug(vectors)

        return sampler_params

    def set_viewpoint(self, camera):
        """
        For a PerspectiveLens, the viewpoint is the front center.
        """
        self.viewpoint = self.front_center

    def __repr__(self):
        disp = "<Lens Object>: lens_type:perspective viewpoint:%s" % (self.viewpoint)
        return disp


class FisheyeLens(Lens):

    """docstring for FisheyeLens"""

    def __init__(self):
        super(FisheyeLens, self).__init__()
        self.fov = 180.0
        self.radius = 1.0
        self.center = None
        self.rotation_matrix = np.eye(3)

    def setup_box_properties(self, camera):
        self.radius = camera.width.max()
        super(FisheyeLens, self).setup_box_properties(camera)
        self.set_viewpoint(camera)

    def new_image(self, camera):
        self.current_image = ImageArray(
            np.zeros((camera.resolution[0]**2, 1,
                      4), dtype='float64', order='C'),
            info={'imtype': 'rendering'})
        return self.current_image

    def _get_sampler_params(self, camera, render_source):
        vp = -arr_fisheye_vectors(camera.resolution[0], self.fov)
        vp.shape = (camera.resolution[0]**2, 1, 3)
        vp = vp.dot(np.linalg.inv(self.rotation_matrix))
        vp *= self.radius
        uv = np.ones(3, dtype='float64')
        positions = np.ones((camera.resolution[0]**2, 1, 3),
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
                 image=image
                 )

        return sampler_params

    def set_viewpoint(self, camera):
        """
        For a PerspectiveLens, the viewpoint is the front center.
        """
        self.viewpoint = camera.position

    def __repr__(self):
        disp = "<Lens Object>: lens_type:fisheye viewpoint:%s fov:%s radius:" %\
                (self.viewpoint, self.fov, self.radius)
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
        lpos = camera.position - pos
        inv_mat = np.linalg.inv(self.rotation_matrix)
        lpos = lpos.dot(self.rotation_matrix)
        #lpos = lpos.dot(self.rotation_matrix)
        mag = (lpos * lpos).sum(axis=1)**0.5
        lpos /= mag[:,None]
        dz = mag / self.radius
        theta = np.arccos(lpos[:,2])
        fov_rad = self.fov * np.pi / 180.0
        r = 2.0 * theta / fov_rad
        phi = np.arctan2(lpos[:,1], lpos[:,0])
        px = r * np.cos(phi)
        py = r * np.sin(phi)
        u = camera.focus.uq
        # dz is distance the ray would travel
        px = (px + 1.0) * res[0] / 2.0
        py = (py + 1.0) * res[1] / 2.0
        px = (u * np.rint(px)).astype("int64")
        py = (u * np.rint(py)).astype("int64")
        return px, py, dz

class SphericalLens(Lens):

    def __init__(self):
        super(SphericalLens, self).__init__()
        self.radius = 1.0
        self.center = None
        self.rotation_matrix = np.eye(3)

    def setup_box_properties(self, camera):
        self.radius = camera.width.max()
        super(SphericalLens, self).setup_box_properties(camera)
        self.set_viewpoint(camera)

    def new_image(self, camera):
        self.current_image = ImageArray(
            np.zeros((camera.resolution[0]*camera.resolution[1], 1,
                      4), dtype='float64', order='C'),
            info={'imtype': 'rendering'})
        return self.current_image

    def _get_sampler_params(self, camera, render_source):
        px = np.linspace(-np.pi, np.pi, camera.resolution[0], endpoint=True)[:,None]
        py = np.linspace(-np.pi/2., np.pi/2., camera.resolution[1], endpoint=True)[None,:]
        
        vectors = np.zeros((camera.resolution[0], camera.resolution[1], 3),
                           dtype='float64', order='C')
        vectors[:,:,0] = np.cos(px) * np.cos(py)
        vectors[:,:,1] = np.sin(px) * np.cos(py)
        vectors[:,:,2] = np.sin(py)

        vectors = vectors * camera.width[0]
        positions = camera.position*vectors.uq + (vectors * 0)
        R1 = get_rotation_matrix(0.5*np.pi, [1,0,0])
        R2 = get_rotation_matrix(-0.5*np.pi, [0,0,1])
        uv = np.dot(R1, camera.unit_vectors)
        uv = np.dot(R2, uv)
        vectors.reshape((camera.resolution[0]*camera.resolution[1], 3))
        vectors = np.dot(vectors, uv)
        #vectors = np.dot(vectors, self.rotation_matrix)
        vectors.reshape((camera.resolution[0], camera.resolution[1], 3))

        if render_source.zbuffer is not None:
            image = render_source.zbuffer.rgba
        else:
            image = self.new_image(camera)
        dummy = np.ones(3, dtype='float64')
        image.shape = (camera.resolution[0]*camera.resolution[1],1,4)
        vectors.shape = (camera.resolution[0]*camera.resolution[1],1,3)
        positions.shape = (camera.resolution[0]*camera.resolution[1],1,3)
        sampler_params = dict(
                vp_pos = positions,
                vp_dir = vectors,
                center = self.back_center.d,
                bounds = (0.0, 1.0, 0.0, 1.0),
                x_vec = dummy,
                y_vec = dummy,
                width = np.zeros(3, dtype="float64"),
                image = image)
        return sampler_params

    def set_viewpoint(self, camera):
        """
        For a PerspectiveLens, the viewpoint is the front center.
        """
        self.viewpoint = camera.position

    def project_to_plane(self, camera, pos, res=None):
        if res is None:
            res = camera.resolution
        # Much of our setup here is the same as in the fisheye, except for the
        # actual conversion back to the px, py values.
        lpos = camera.position - pos
        #inv_mat = np.linalg.inv(self.rotation_matrix)
        #lpos = lpos.dot(self.rotation_matrix)
        mag = (lpos * lpos).sum(axis=1)**0.5
        lpos /= mag[:,None]
        # originally:
        #  the x vector is cos(px) * cos(py)
        #  the y vector is sin(px) * cos(py)
        #  the z vector is sin(py)
        # y / x = tan(px), so arctan2(lpos[:,1], lpos[:,0]) => px
        # z = sin(py) so arcsin(z) = py
        # px runs from -pi to pi
        # py runs from -pi/2 to pi/2
        px = np.arctan2(lpos[:,1], lpos[:,0])
        py = np.arcsin(lpos[:,2])
        dz = mag / self.radius
        u = camera.focus.uq
        # dz is distance the ray would travel
        px = ((-px + np.pi) / (2.0*np.pi)) * res[0]
        py = ((-py + np.pi/2.0) / np.pi) * res[1]
        px = (u * np.rint(px)).astype("int64")
        py = (u * np.rint(py)).astype("int64")
        return px, py, dz


lenses = {'plane-parallel': PlaneParallelLens,
          'perspective': PerspectiveLens,
          'fisheye': FisheyeLens,
          'spherical': SphericalLens}
