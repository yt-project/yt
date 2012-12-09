"""
A class that manages the coordinate system for orientable data
containers and cameras.

Author: Nathan Goldbaum <goldbaum@ucolick.org>
Affiliation: UCSC Astronomy
License:
  Copyright (C) 2008-2011 Matthew Turk.  All Rights Reserved.

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

import numpy as np

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
        if np.all(north_vector == normal_vector):
            mylog.error("North vector and normal vector are the same.  Disregarding north vector.")
            north_vector = None
        if north_vector is not None: self.steady_north = True
        self.north_vector = north_vector
        self._setup_normalized_vectors(normal_vector, north_vector)

    def _setup_normalized_vectors(self, normal_vector, north_vector):
        # Now we set up our various vectors
        normal_vector /= np.sqrt( np.dot(normal_vector, normal_vector))
        if north_vector is None:
            vecs = np.identity(3)
            t = np.cross(normal_vector, vecs).sum(axis=1)
            ax = t.argmax()
            east_vector = np.cross(vecs[ax,:], normal_vector).ravel()
            # self.north_vector must remain None otherwise rotations about a fixed axis will break.  
            # The north_vector calculated here will still be included in self.unit_vectors.
            north_vector = np.cross(normal_vector, east_vector).ravel()
        else:
            if self.steady_north:
                north_vector = north_vector - np.dot(north_vector,normal_vector)*normal_vector
            east_vector = np.cross(north_vector, normal_vector).ravel()
        north_vector /= np.sqrt(np.dot(north_vector, north_vector))
        east_vector /= np.sqrt(np.dot(east_vector, east_vector))
        self.normal_vector = normal_vector
        self.unit_vectors = [east_vector, north_vector, normal_vector]
        self.inv_mat = np.linalg.pinv(self.unit_vectors)
        
    def switch_orientation(self, normal_vector=None, north_vector=None):
        r"""Change the view direction based on any of the orientation parameters.

        This will recalculate all the necessary vectors and vector planes related
        to an orientable object.

        Parameters
        ----------
        normal_vector: array_like, optional
            The new looking vector.
        north_vector : array_like, optional
            The 'up' direction for the plane of rays.  If not specific,
            calculated automatically.
        """
        if north_vector is None:
            north_vector = self.north_vector
        if normal_vector is None:
            normal_vector = self.normal_vector
        self._setup_normalized_vectors(normal_vector, north_vector)

        
