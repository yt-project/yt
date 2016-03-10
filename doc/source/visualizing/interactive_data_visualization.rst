.. _interactive_data_visualization:

Interactive Data Visualization
==============================

Installation
^^^^^^^^^^^^

In version 3.3 we introduced experimental, hardware accelerated, interactive
volume renderer based on OpenGL. In order to use it you need to install
`PyOpenGL <https://pypi.python.org/pypi/PyOpenGL>`_ and `cyglfw3
<https://pypi.python.org/pypi/cyglfw3/`_ along with their respective
dependencies.

Using interactive renderer
^^^^^^^^^^^^^^^^^^^^^^^^^^

You can simply pass dataset to :meth:`~yt.interactive_render`. By default
it will load all data and render gas density:

.. python-script::

    import yt
    
    ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    yt.interactive_render(ds)

Alternatively you can provide data object as a main argument to
:meth:`~yt.interactive_render` if your dataset is too big to fit GPU memory.

.. python-script::

    import yt

    ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    sp = ds.sphere("max", (0.1, "Mpc"))

    cam_pos = ds.arr([0.1, 0.1, 0.1], "Mpc").in_units("code_length")
    yt.interactive_render(sp, field="pressure", cam_position=cam_pos)

Successful call to :meth:`~yt.interactive_render` should create new
window called *vol_render*. By default it renders maximum intensity projection
of your data. You can move the camera around by holding left mouse button and
moving it. Pressing *h* key will print all available key binding in the
terminal window.
