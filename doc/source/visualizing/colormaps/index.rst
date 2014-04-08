.. _colormaps:

Colormaps
=========

There are several colormaps available for yt.  yt includes all of the 
matplotlib colormaps as well for nearly all functions.  Individual visualization
functions usually allow you to specify a colormap with the ``cmap`` flag.
There are a small number of functions (mostly contained in the image_writer 
module; e.g. write_bitmap, write_image, write_projection, etc.), which do 
not load the matplotlib infrastructure and can only access the colormaps 
native to yt.  

Here is a chart of all of the colormaps available.  In addition to each 
colormap displayed here, you can access its "reverse" by simply appending a 
``"_r"`` to the end of the colormap name.

All Colormaps (including matplotlib)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. image:: ../_images/all_colormaps.png
   :width: 512

Native yt Colormaps
~~~~~~~~~~~~~~~~~~~

.. image:: ../_images/native_yt_colormaps.png
   :width: 512

Displaying Colormaps Locally
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To display the most up to date colormaps locally, you can run:

.. code-block:: python

    from yt.mods import *
    show_colormaps()

or to output just the colormaps native to yt to an image file, try:

.. code-block:: python

    from yt.mods import *
    show_colormaps(subset = "yt_native", filename = "yt_native.png")

Applying a Colormap to your Rendering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All of the visualization functions in yt have a keyword allowing you to
manually specify a specific colormap.  For example:

.. code-block:: python

    write_image(im, "output.png", cmap_name = 'jet')

If you're using the Plot Window interface (e.g. SlicePlot, ProjectionPlot, 
etc.), it's even easier than that.  Simply create your rendering, and you
can quickly swap the colormap on the fly after the fact with the ``set_cmap``
callback:

.. code-block:: python

    pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")
    p = ProjectionPlot(pf, "z", "density")

    p.set_cmap(field="density", cmap='jet')
    p.save('proj_with_jet_cmap.png')

    p.set_cmap(field="density", cmap='hot')
    p.save('proj_with_hot_cmap.png')

For more information about the callbacks available to Plot Window objects, 
see :ref:`callbacks`.

Examples of Each Colormap
~~~~~~~~~~~~~~~~~~~~~~~~~

To give the reader a better feel for how a colormap appears once it is applied
to a dataset, below we provide a library of identical projections of an 
isolated galaxy where only the colormap has changed.  They use the sample 
dataset "IsolatedGalaxy" available at 
`http://yt-project.org/data <http://yt-project.org/data>`_.

.. yt_colormaps:: cmap_images.py
