.. _geographic_projections_and_transforms:

Geographic Projections and Transforms
=====================================

Geographic data that is on a sphere can be visualized by projecting that data
onto a representation of that sphere flattened into 2d space. There exist a
number of projection types, which can be found in the `the cartopy
documentation <https://scitools.org.uk/cartopy/docs/v0.15/crs/projections.html>`_.
yt now supports these projection types for geographically loaded data.
Underlying data is assumed to have an underlying projection type of `PlateCarree`. If
your data is not of this form, feel free to open an issue or file a pull
request on the yt github page for this feature.

It should be noted that
these projections are not the same as yt's ProjectionPlot. For more information
on yt's projection plots, see :ref:`projection-types`.

Installation
^^^^^^^^^^^^

In order to access the geographic projection functionality, you will need to have an
installation of cartopy available on your machine.

.. code-block:: bash

    conda install cartopy

Using Basic Projections
^^^^^^^^^^^^^^^^^^^^^^^

All of the transforms available in Cartopy v0.15 are accessible with this
functionality.

The next few examples will use the GEOS dataset. For details about loading this
data, please see the geographic projections cookbook.

If a geographic dataset is loaded without any defined projection the default
option of `PlateCarree` will be displayed.

.. code-block:: python

    ds = yt.load_uniform_grid(data, sizes, 1.0, geometry=("geographic", dims),
    bbox=bbox)
    p = yt.SlicePlot(ds, "altitude", 'AIRDENS')

If an option other than `PlateCarree` is desired, the plot projection type can
be set manually with the `set_mpl_projection` function. This next code block
sets the projection to a `Robinson` projection from the default `PlateCarree`.

.. code-block:: python

    p.set_mpl_projection('Robinson')
    p.show()

The axes attributes of the plot can be accessed to add in annotations, such as
coastlines. The axes are matplotlib `GeoAxes` so any of the annotations
available with matplotlib should be available for customization. Here a
`Robinson` plot is made with coastline annotations.

.. code-block:: python

    p.set_mpl_projection('Robinson')
    p._setup_plots()
    p.plots['AIRDENS'].axes.set_global()
    p.plots['AIRDENS'].axes.coastlines()
    p.show()

`p._setup_plots()` is required here to access the plot axes. When a new
projection is called the plot axes are reset and are not available unless set
up again.

Additional arguments can be passed to the projection function for further
customization. If additional arguments are desired, then rather than passing a
string of the projection name, one would pass a 2 or 3-item tuple.

The function set_mpl_projection can take one of three input types:

.. code-block:: python

    set_mpl_projection('ProjectionType')
    set_mpl_projection(('ProjectionType', (args)))
    set_mpl_projection(('ProjectionType', (args), {kwargs}))

Further examples of using the geographic transforms with this dataset
can be found in :ref:`cookbook-geographic_projections`.
