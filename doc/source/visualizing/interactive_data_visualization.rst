.. _interactive_data_visualization:

Interactive Data Visualization
==============================

In version 3.3 we introduced an experimental, hardware accelerated, interactive
volume renderer based on OpenGL. It allows to load a grid dataset into a GPU
memory as a series of 3D textures and subsequently view it using an interactive
window, which represents "camera", that can be rotated or zoomed. Color of the each
pixel in the window is evaluated by applying a small gpu program called
"fragment shader" to a set of rays that traverse the volume parallel to the
camera's line of sight. The currently implemented shaders allow to compute a
maximum value along the line of sight (also know as a `Maximum Intensity
Projection <https://en.wikipedia.org/wiki/Maximum_intensity_projection>`_) and
an unweighted integration of the grid values for a given field along the line of
sight. See :ref:`projection-types` for more information.

A comprehensive description of the OpenGL volume rendering is beyond the scope
of this document. However, a more detailed explanation can be found in `this
guide <https://open.gl/>`_.

Installation
^^^^^^^^^^^^

In order to use Interactive Data Visualization (IDV) you need to install
`PyOpenGL <https://pypi.python.org/pypi/PyOpenGL>`_ and `cyglfw3
<https://pypi.python.org/pypi/cyglfw3/>`_ along with their respective
dependencies, e.g. `glfw3 <http://www.glfw.org/>`_ is required to be installed
before you can ``pip install cyglfw3``. Please carefully read installation
instructions provided on pypi pages of both packages. 

If you are using conda, ``cyglfw3`` is provided in our conda channel
(``pyopengl`` is shipped by Continuum already) and can be installed via:

.. code-block:: bash

    conda install -c http://use.yt/with_conda/ cyglfw3 pyopengl

Using the interactive renderer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can simply pass dataset to :meth:`~yt.interactive_render`. By default
it will load all data and render gas density:

.. code-block:: python

    import yt
    
    ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    yt.interactive_render(ds)

Alternatively you can provide a data object as a first argument to
:meth:`~yt.interactive_render` if your dataset is too big to fit GPU memory:

.. code-block:: python

    import yt

    ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    sp = ds.sphere("max", (0.1, "Mpc"))

    cam_pos = ds.arr([0.1, 0.1, 0.1], "Mpc").in_units("code_length")
    yt.interactive_render(sp, field="pressure", cam_position=cam_pos,
                          window_size=(512, 512))

A successful call to :meth:`~yt.interactive_render` should create a new window
called *vol_render*. 

.. image:: _images/idv.jpg
   :width: 1000

By default it renders a Maximum Intensity Projection of the density field. The
rendering can be dynamically modified using the following keybindings:

1
   Switch to MIP fragment shader
2
   Switch to integration fragement shader
L
   Switch between linear and logarithmic scales
W
   Zoom in the camera
S
   Zoom out the camera
C
   Change the colormap

Pressing the *h* key will print all the available key bindings in a terminal window.
The camera can be moved around by holding a left mouse button while moving the mouse.

More advanced initialization of interactive volume renderer can be found in
:ref:`cookbook-opengl_vr`.
