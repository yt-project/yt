import numpy as na

from yt.funcs import *
from yt.utilities.math_utils import get_rotation_matrix

class Orientation:
    def __init__(self, normal_vector, north_vector=None, steady_north=False):
        r"""An object that returns a set of basis vectors for orienting
        cameras a data containers.

        Parameters
        ----------
        center        : array_like
           The current "center" of the view port -- the normal_vector connects
           the center and the origin
        normal_vector : array_like
           A vector normal to the image plane
        north_vector  : array_like, optional
           The 'up' direction to orient the image plane.  
           If not specified, gets calculated automatically
        steady_north  : bool, optional
           Boolean to control whether to normalize the north_vector
           by subtracting off the dot product of it and the normal 
           vector.  Makes it easier to do rotations along a single
           axis.  If north_vector is specified, is switched to
           True.  Default: False
           
        """
        self.steady_north = steady_north
        if na.all(north_vector == normal_vector):
            mylog.error("North vector and normal vector are the same.  Disregarding north vector.")
            north_vector = None
        if north_vector is not None: self.steady_north = True
        self._setup_normalized_vectors(normal_vector, north_vector)

    def _setup_normalized_vectors(self, normal_vector, north_vector):
        # Now we set up our various vectors
        normal_vector /= na.sqrt( na.dot(normal_vector, normal_vector))
        if north_vector is None:
            vecs = na.identity(3)
            t = na.cross(normal_vector, vecs).sum(axis=1)
            ax = t.argmax()
            north_vector = na.cross(vecs[ax,:], normal_vector).ravel()
        else:
            if self.steady_north:
                north_vector = north_vector - na.dot(north_vector,normal_vector)*normal_vector
        north_vector /= na.sqrt(na.dot(north_vector, north_vector))
        east_vector = -na.cross(north_vector, normal_vector).ravel()
        east_vector /= na.sqrt(na.dot(east_vector, east_vector))
        self.normal_vector = normal_vector
        self.north_vector = north_vector
        self.unit_vectors = [east_vector, north_vector, normal_vector]
        self.inv_mat = na.linalg.pinv(self.unit_vectors)

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
        self._setup_normalized_vectors(normal_vector, north_vector)


    def switch_orientation(self, normal_vector=None, center=None, north_vector=None):
        r"""Change the view direction based on any of the orientation parameters.

        This will recalculate all the necessary vectors and vector planes related
        to a camera with new normal vectors, widths, centers, or north vectors.

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
        if north_vector is None:
            north_vector = self.north_vector
        if normal_vector is None:
            normal_vector = self.front_center-center
        self._setup_normalized_vectors(normal_vector, north_vector)

        
