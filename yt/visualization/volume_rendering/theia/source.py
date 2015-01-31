#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

class TheiaSource:
    r"""
        TheiaSource

        This is a wrapper class to create an interface for consistent access to 
        raycasters and the volume they act on, as well as the returned results.


	The raycaster shold be a GPU based raycasting algorithm to act on
        the volume. The raycaster class must provide the following functions:
            send_volume_to_gpu(volume)
                Accepts a volume and moves the data onto the graphics card 
                memory.
            set_matrix(camera.get_matrix()) 
                Accepts a modelview matrix (4x4 numpy matrix) from the 
                Theia Camera class and sets the view on the volume
            cast()
                This must cause the ray caster to act on the volume
                in the graphics cards memory, sampling the volume and 
                filling the surface with results

            get_surface()
        	This must return a 2d numpy array with the results of
                the last raycast


        Any UI's using a TheiaSource can call the functions:

        update()
        	This will call the raycaster's cast function

        get_results()
        	This returns a 2d numpy array with the results of the latest update

        Parameters
        ----------
        volume   : 
            The data the camera will volume render. 
        raycaster : 
            A ray casting algorithm conforming to the requirements
            detailed above

        Example:

        from yt.visualization.volume_rendering.theia.source import TheiaSource
        from yt.visualization.volume_rendering.algorithms.front_to_back import FrontToBackRaycaster
        import numpy as np

        #load a 3d numpy array of float32 
        volume = np.load("/home/bogert/log_densities_1024.npy")
        raycaster = FrontToBackRaycaster()

        ts = TheiaSource(volume = volume, raycaster = raycaster)

        ts.update()

        image_array = ts.get_results()
    
    """
    def __init__(self, volume = None, raycaster = None):
        #Values to control screen resolution
        self.raycaster = None
        self.volume    = None

        self.set_volume(volume)
        self.set_raycaster(raycaster)

    def set_volume(self, volume = None):
        r"""This will set the TheiaSource volume to
            the volume supplied. If the TheiaSource 
            has a raycaster object it will send the
            volume to the gpu. 

        Parameters
        ----------
        volume   : 
            The volumetric data to send to the raycaster
        """
        self.volume = volume
        if (self.volume != None  and self.raycaster != None) :
            self.raycaster.send_volume_to_gpu(self.volume)

    def set_raycaster(self, raycaster = None):
        r"""This will set the TheiaSource raycaster to
            the raycaster supplied. If the TheiaSource 
            has a volume set previously it will send the
            volume to the gpu via the raycaster.

        Parameters
        ----------
        raycaster   : 
            The raycaster used to render the volume.
        """
        self.raycaster = raycaster
        if (self.volume != None  and self.raycaster != None) :
            self.raycaster.send_volume_to_gpu(self.volume)

    def update(self, matrix = None):
        r"""This will update the raycaster with any 
            matrix information if given and then 
            call the raycaster's cast() function to
            render the volume.

        Parameters
        ----------
        matrix   : 4x4 numpy matrix
            The modelview matrix used to define the
            view onto the volume
        """
        if (matrix != None) :
            self.raycaster.set_matrix(matrix)
        self.raycaster.cast()

    def get_results(self):
        r"""This will return the results of the 
            raycasting after the update function
            is called as a 2d numpy  array.

        """
        return self.raycaster.get_surface()
    
