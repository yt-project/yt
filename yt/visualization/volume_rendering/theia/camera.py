#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
  
#Matix helpers for camera control
import yt.visualization.volume_rendering.theia.helpers.matrix_transformation as mth

import numpy as np

class Camera:
    r"""
    This is a basic camera intended to be a base class. 
    The camera provies utilities for controling the view point onto the data.


    Parameters
    ----------
    pos   : numpy array
        The 3d position vector of the camera
    rot   : numpy array
        The 3d rotation vector of the camera
    scale : float
        The scalar acting on camera position to zoom

    Example:

    from yt.visualization.volume_rendering.cameras.camera import Camera

    cam = Camera()
    
    cam.zoom(10.0)
    cam.rotateX(np.pi)
    cam.rotateY(np.pi)
    cam.rotateZ(np.pi)
    cam.translateX(-0.5)
    cam.translateY(-0.5)
    cam.translateZ(-0.5)

    view  = cam.get_matrix()

    """
    def __init__(self, rot = np.array([0.0, 0.0, 0.0]), scale = 1.0,
                       pos = np.array([0.0, 0.0, 0.0]), ):
        #Values to control camera
        self.rot   = rot
        self.scale = scale
        self.pos   = pos

        #TODO : users should have control over perspective frustrum
        self.stdmatrix  = mth.perspective( [0.0, 0.0, 0.25] )


    def zoom(self, percentage = 10.0):
        r"""This will zoom the camera by percentage specified.

        Parameters
        ----------
        percentage   : float
            If percentage is postive the camera will zoom in, if negative
            then the camera will zoom out. Cannot exceed 100%.
        """
        if (percentage > 100.0) :
            percentage = 100.0
        self.scale = ((100.0 - percentage)/100.0)*self.scale
      
    def rotateX(self, degrees = 0.1):
        r"""This will rotate the camera around its X axis.

        Parameters
        ----------
        degrees   : float
            The amount of rotation in specified in radians.
        """
        self.rot[0] += degrees

    def rotateY(self, degrees = 0.1):
        r"""This will rotate the camera around its Y axis.

        Parameters
        ----------
        degrees   : float
            The amount of rotation in specified in radians.
        """
        self.rot[1] += degrees

    def rotateZ(self, degrees = 0.1):
        r"""This will rotate the camera around its Z axis.

        Parameters
        ----------
        degrees   : float
            The amount of rotation in specified in radians.
        """
        self.rot[2] += degrees

    def translateX(self, dist = 0.1):
        r"""This will move the camera's position along its X axis.

        Parameters
        ----------
        dist   : float
            The amount of movement. **This unit is unknown!**
        """
        self.pos[0] += dist

    def translateY(self, dist = 0.1):
        r"""This will move the camera's position along its Y axis.

        Parameters
        ----------
        dist   : float
            The amount of movement. **This unit is unknown!**
        """
        self.pos[1] += dist

    def translateZ(self, dist = 0.1):
        r"""This will move the camera's position along its Z axis.

        Parameters
        ----------
        dist   : float
            The amount of movement. **This unit is unknown!**
        """
        self.pos[2] += dist

    def get_matrix(self):
        r"""This will calculate the curent modelview matrix
	    from the camera's position, rotation, and zoom scale.  

        Parameters
        ----------
        dist   : float
            The amount of movement. **This unit is unknown!**
        """
        return mth.rotate_x(self.rot[0]) * mth.rotate_y(self.rot[1]) * mth.rotate_z(self.rot[2]) * mth.scale((self.scale, self.scale, self.scale)) * mth.perspective([0.0, 0.0, 0.25]) *mth.translate((self.pos[0], self.pos[1], self.pos[2])) 
