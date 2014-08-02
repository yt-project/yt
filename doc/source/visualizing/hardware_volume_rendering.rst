.. _hardware_volume_rendering:

Hardware Volume Rendering on NVidia Graphics cards
--------------------------------------------------

Theia is a hardware volume renderer that takes advantage of NVidias CUDA language
to peform ray casting with GPUs instead of the CPU. 

Only unigrid rendering is supported, but yt provides a grid mapping function
to get unigrid data from amr or sph formats, see :ref:`extract_frb`.

System Requirements
+++++++++++++++++++

Nvidia graphics card - The memory limit of the graphics card sets the limit
                       on the size of the data source.

CUDA 5 or later and

The environment variable CUDA_SAMPLES must be set pointing to
the common/inc samples shipped with CUDA. The following shows an example
in bash with CUDA 5.5 installed in /usr/local :

    export CUDA_SAMPLES=/usr/local/cuda-5.5/samples/common/inc

PyCUDA must also be installed to use Theia. 

PyCUDA can be installed following these instructions :

    git clone --recursive http://git.tiker.net/trees/pycuda.git

    python configure.py
    python setup.py install


Tutorial
++++++++

Currently rendering only works on uniform grids. Here is an example
on a 1024 cube of float32 scalars.

.. code-block:: python

   from yt.visualization.volume_rendering.theia.scene import TheiaScene
   from yt.visualization.volume_rendering.algorithms.front_to_back import FrontToBackRaycaster
   import numpy as np

   #load 3D numpy array of float32
   volume = np.load("/home/bogert/log_densities_1024.npy")

   scene = TheiaScene( volume = volume, raycaster = FrontToBackRaycaster() )

   scene.camera.rotateX(1.0)
   scene.update()

   surface = scene.get_results()
   #surface now contains an image array 2x2 int32 rbga values

.. _the-theiascene-interface:

The TheiaScene Interface
++++++++++++++++++++++++

A TheiaScene object has been created to provide a high level entry point for
controlling the raycaster's view onto the data. The class
:class:`~yt.visualization.volume_rendering.theia.TheiaScene` encapsulates a
Camera object and a TheiaSource that intern encapsulates a volume. The
:class:`~yt.visualization.volume_rendering.theia.Camera` provides controls for
rotating, translating, and zooming into the volume.  Using the
:class:`~yt.visualization.volume_rendering.theia.TheiaSource` automatically
transfers the volume to the graphic's card texture memory.

Example Cookbooks
+++++++++++++++++

OpenGL Example for interactive volume rendering:

.. literalinclude:: opengl_volume_rendering.py

.. warning::  Frame rate will suffer significantly from stereoscopic rendering.
              ~2x slower since the volume must be rendered twice.

OpenGL Stereoscopic Example: 

.. literalinclude:: opengl_stereo_volume_rendering.py

Pseudo-Realtime video rendering with ffmpeg:

.. literalinclude:: ffmpeg_volume_rendering.py
