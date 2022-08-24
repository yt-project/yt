import numpy as np

from yt.units.yt_array import YTArray
from yt.utilities.exceptions import YTException


def _aligned(a, b):
    aligned_component = np.abs(np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b))
    return np.isclose(aligned_component, 1.0, 1.0e-13)


def _validate_unit_vectors(normal_vector, north_vector):

    # Make sure vectors are unitless
    if north_vector is not None:
        north_vector = YTArray(north_vector, "", dtype="float64")
    if normal_vector is not None:
        normal_vector = YTArray(normal_vector, "", dtype="float64")

    if not np.dot(normal_vector, normal_vector) > 0:
        raise YTException("normal_vector cannot be the zero vector.")
    if north_vector is not None and _aligned(north_vector, normal_vector):
        raise YTException("normal_vector and north_vector cannot be aligned.")

    return normal_vector, north_vector


class Orientation:
    def __init__(self, normal_vector, north_vector=None, steady_north=False):
        r"""An object that returns a set of basis vectors for orienting
        cameras a data containers.

        Parameters
        ----------
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

        normal_vector, north_vector = _validate_unit_vectors(
            normal_vector, north_vector
        )
        self.steady_north = steady_north
        if north_vector is not None:
            self.steady_north = True
        self.north_vector = north_vector
        self._setup_normalized_vectors(normal_vector, north_vector)
        if self.north_vector is None:
            self.north_vector = self.unit_vectors[1]

    def _setup_normalized_vectors(self, normal_vector, north_vector):
        normal_vector, north_vector = _validate_unit_vectors(
            normal_vector, north_vector
        )
        # Now we set up our various vectors
        normal_vector /= np.sqrt(np.dot(normal_vector, normal_vector))
        if north_vector is None:
            vecs = np.identity(3)
            t = np.cross(normal_vector, vecs).sum(axis=1)
            ax = t.argmax()
            east_vector = np.cross(vecs[ax, :], normal_vector).ravel()
            # self.north_vector must remain None otherwise rotations about a fixed axis
            # will break. The north_vector calculated here will still be included
            # in self.unit_vectors.
            north_vector = np.cross(normal_vector, east_vector).ravel()
        else:
            if self.steady_north or (np.dot(north_vector, normal_vector) != 0.0):
                north_vector = (
                    north_vector - np.dot(north_vector, normal_vector) * normal_vector
                )
            east_vector = np.cross(north_vector, normal_vector).ravel()
        north_vector /= np.sqrt(np.dot(north_vector, north_vector))
        east_vector /= np.sqrt(np.dot(east_vector, east_vector))
        self.normal_vector = normal_vector
        self.north_vector = north_vector
        self.unit_vectors = YTArray([east_vector, north_vector, normal_vector], "")
        self.inv_mat = np.linalg.pinv(self.unit_vectors)
