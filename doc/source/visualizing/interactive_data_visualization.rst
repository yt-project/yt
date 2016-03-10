.. _interactive_data_visualization:

Interactive Data Visualization
==============================

Installation
^^^^^^^^^^^^

In version 3.3 we introduced experimental, hardware accelerated, interactive
volume renderer based on OpenGL. In order to use it you need to install
`PyOpenGL <https://pypi.python.org/pypi/PyOpenGL>`_ and `cyglfw3
<https://pypi.python.org/pypi/cyglfw3/>`_ along with their respective
dependencies. Both packages are available in our conda channel and can be
installed via:

.. code-block:: bash

    conda install -c http://use.yt/with_conda/ cyglfw3 pyopengl

Using interactive renderer
^^^^^^^^^^^^^^^^^^^^^^^^^^

You can simply pass dataset to :meth:`~yt.interactive_render`. By default
it will load all data and render gas density:

.. code-block:: python

    import yt
    
    ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    yt.interactive_render(ds)

Alternatively you can provide data object as a main argument to
:meth:`~yt.interactive_render` if your dataset is too big to fit GPU memory.

.. code-block:: python

    import yt

    ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    sp = ds.sphere("max", (0.1, "Mpc"))

    cam_pos = ds.arr([0.1, 0.1, 0.1], "Mpc").in_units("code_length")
    yt.interactive_render(sp, field="pressure", cam_position=cam_pos,
                          window_size=(512, 512))

Successful call to :meth:`~yt.interactive_render` should create new window
called *vol_render*. 

.. image:: _images/idv.jpg
   :width: 1000

By default it renders maximum intensity of your data.
Camera can be moved around by holding left mouse button while moving the mouse.
Apart from doing maximum intensity rendering, it's possible to create projection
along the line of sight (press *2*) which is equivalent to
:meth:`~yt.visualization.volume_rendering.off_axis_projection.off_axis_projection`.
Pressing *h* key will print all available key bindings in the terminal window.
More advanced initialization of interactive volume renderer can be found in
:ref:`cookbook-opengl-vr`.
