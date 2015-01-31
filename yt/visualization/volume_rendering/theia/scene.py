#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.visualization.volume_rendering.theia.source import TheiaSource
from yt.visualization.volume_rendering.theia.camera import Camera

class TheiaScene:
    r"""
        TheiaScene

        This is the highest level entry point into Theia.
        This class is a container for a Theia source and 
        a camera for controlling a view onto the source.

        Call :
  	    update()
        To perform the volume rendering.
        Call :
            get_results()
        To get the resulting image array.
 

        Parameters
        ----------
        ts : TheiaSource
            The container for the raycasting algorithm and the volumetric data    
        volume : 
            If a TheiaSource is not specified a volme can be used to initialize
            the scene. However, rendering will only occur when a raycaster
            is also provided.

        camera : 
            A camera for controling view point. Must be able to provide a 4x4
            matrix to the raycaster via get_matrix()
            If no camera is provided the default Theia Camera is used.

        Example:

        from yt.visualization.volume_rendering.theia.scene import TheiaScene
        from yt.visualization.volume_rendering.algorithms.front_to_back import FrontToBackRaycaster
        import numpy as np

        #load 3D numpy array of float32
        volume = np.load("/home/bogert/log_densities_1024.npy")

        scene = TheiaScene( volume = volume, raycaster = FrontToBackRaycaster() )

        scene.camera.rotateX(1.0)
        scene.update()

        surface = scene.get_results()
        #now do something with surface
    
    """
    def __init__(self, ts = None, camera = None, volume = None, raycaster = None) :
        #Values to control screen resolution
        if ts != None :
            self.source = ts
        else : 
            self.source = TheiaSource(volume  = volume, raycaster = raycaster)

        if camera == None :
              self.camera = Camera()
        else :
              self.camera = camera

    def get_raycaster(self):
        """
            This is for programmers who would prefer calling :
                scene.get_raycaster() 
            instead of :
                scene.source.raycaster
	"""
        return self.source.raycaster

    def update(self):
        """
            This will send the camera's modelview matrix to
            the raycaster and call the raycaster's cast function
            to perform the volume rendering. 
        """
        self.source.update(self.camera.get_matrix())

    def get_results(self):
        """
            This will cause the results of the volume rendering to be 
            copied back from the GPU and then returned as a 2d numpy array. 
        """
        return self.source.get_results()

