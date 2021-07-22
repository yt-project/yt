.. _geographic_projections_and_transforms:

Geographic Projections and Transforms
=====================================

Geographic data that is on a sphere can be visualized by projecting that data
onto a representation of that sphere flattened into 2d space. There exist a
number of projection types, which can be found in the `the cartopy
documentation <https://scitools.org.uk/cartopy/docs/latest/crs/projections.html>`_.
With support from `cartopy <https://scitools.org.uk/cartopy/docs/latest/>`_,
``yt`` now supports these projection
types for geographically loaded data.
Underlying data is assumed to have a transform of `PlateCarree
<https://scitools.org.uk/cartopy/docs/latest/crs/projections.html#platecarree>`__,
which is data on a flattened, rectangular, latitude/longitude grid. This is a
a typical format for geographic data.

The distinction between the data transform and projection is worth noting. The data
transform is what system your data is defined with and the data projection is
what the resulting plot will display. For more information on this difference,
refer to `the cartopy documentation on these differences
<https://scitools.org.uk/cartopy/docs/latest/tutorials/understanding_transform.html>`_.
If your data is not of this form, feel free to open an issue or file a pull
request on the ``yt`` github page for this feature.

It should be noted that
these projections are not the same as yt's ProjectionPlot. For more information
on yt's projection plots, see :ref:`projection-types`.

.. _install-cartopy:

Installing Cartopy
^^^^^^^^^^^^^^^^^^

In order to access the geographic projection functionality, you will need to have an
installation of ``cartopy`` available on your machine. If you're using conda as
your package manager, this should be sufficient:

.. code-block:: bash

    conda install cartopy

If you're on a mac and are using pip, there can be conflicts with the GEOS
library and cartopy / cartopy dependencies. To avoid these issues use your
package manager of choice to install proj4 and geos (``proj`` and ``geos``).
Following that, build cartopy and shapely from source.
For example, a user using homebrew and pip
would execute the following commands:

.. code-block:: bash

    $ brew install proj geos
    $ brew upgrade proj geos
    $ python -m pip install --no-binary :all: shapely cartopy


On ubuntu you'll need to install the following packages: ``libproj-dev``,
``proj-data``, ``proj-bin``, and ``libgeos-dev``.

Using Basic Transforms
^^^^^^^^^^^^^^^^^^^^^^^

As mentioned above, the default data transform is assumed to be of `PlateCarree
<https://scitools.org.uk/cartopy/docs/latest/crs/projections.html#platecarree>`__,
which is data on a flattened, rectangular, latitude/longitude grid. To set
something other than ``PlateCarree``, the user can access the dictionary in the coordinate
handler that defines the coordinate transform to change the default transform
type. Because the transform
describes the underlying data coordinate system, the loaded dataset will carry
this newly set attribute and all future plots will have the user-defined data
transform. Also note that the dictionary is ordered by axis type. Because
slicing along the altitude may differ from, say, the latitude axis, we may
choose to have different transforms for each axis.

.. code-block:: python

    ds = yt.load_uniform_grid(data, sizes, 1.0, geometry=("geographic", dims), bbox=bbox)
    ds.coordinates.data_transform["altitude"] = "Miller"
    p = yt.SlicePlot(ds, "altitude", "AIRDENS")

In this example, the ``data_transform`` kwarg has been changed from its default
of ``PlateCarree`` to ``Miller``. You can check that you have successfully changed
the defaults by inspecting the ``data_transform`` and ``data_projection`` dictionaries
in the coordinate
handler. For this dataset, that would be accessed by:

.. code-block:: python

    print(ds.coordinates.data_transform["altitude"])
    print(ds.coordinates.data_projection["altitude"])



Using Basic Projections
^^^^^^^^^^^^^^^^^^^^^^^

All of the transforms available in ``Cartopy`` v0.15 and above are accessible
with this functionality.

The next few examples will use a GEOS dataset accessible from the ``yt`` data
downloads page. For details about loading this data, please
see :ref:`cookbook-geographic_projections`.

If a geographic dataset is loaded without any defined projection the default
option of ``Mollweide`` will be displayed.

.. code-block:: python

    ds = yt.load_uniform_grid(data, sizes, 1.0, geometry=("geographic", dims), bbox=bbox)
    p = yt.SlicePlot(ds, "altitude", "AIRDENS")

If an option other than ``Mollweide`` is desired, the plot projection type can
be set with the ``set_mpl_projection`` function. The next code block illustrates how to
set the projection to a ``Robinson`` projection from the default ``PlateCarree``.

.. code-block:: python

    ds = yt.load_uniform_grid(data, sizes, 1.0, geometry=("geographic", dims), bbox=bbox)
    p = yt.SlicePlot(ds, "altitude", "AIRDENS")
    p.set_mpl_projection("Robinson")
    p.show()

The axes attributes of the plot can be accessed to add in annotations, such as
coastlines. The axes are matplotlib ``GeoAxes`` so any of the annotations
available with matplotlib should be available for customization. Here a
``Robinson`` plot is made with coastline annotations.

.. code-block:: python

    p.set_mpl_projection("Robinson")
    p._setup_plots()
    p.plots["AIRDENS"].axes.set_global()
    p.plots["AIRDENS"].axes.coastlines()
    p.show()

``p._setup_plots()`` is required here to access the plot axes. When a new
projection is called the plot axes are reset and are not available unless set
up again.

Additional arguments can be passed to the projection function for further
customization. If additional arguments are desired, then rather than passing a
string of the projection name, one would pass a 2 or 3-item tuple, the first
item of the tuple corresponding to a string of the transform name, and the
second and third items corresponding to the args and kwargs of the transform,
respectively.

Alternatively, a user can pass a transform object rather than a string or tuple.
This allows for users to
create and define their own transforms, beyond what is available in cartopy.
The type must be a cartopy GeoAxes object or a matplotlib transform object. For
creating custom transforms, see `the matplotlib example
<https://matplotlib.org/examples/api/custom_projection_example.html>`_.

The function ``set_mpl_projection()`` accepts several input types for varying
levels of customization:

.. code-block:: python

    set_mpl_projection("ProjectionType")
    set_mpl_projection(("ProjectionType", (args)))
    set_mpl_projection(("ProjectionType", (args), {kwargs}))
    set_mpl_projection(cartopy.crs.PlateCarree())

Further examples of using the geographic transforms with this dataset
can be found in :ref:`cookbook-geographic_projections`.
