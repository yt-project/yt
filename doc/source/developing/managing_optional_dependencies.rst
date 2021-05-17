How to manage soft (optional) dependencies
------------------------------------------

We might sometimes rely on heavy external libraries to support some features
outside of yt's core. Typically, many frontends require HDF5 support, provided
by the ``h5py`` library, but many users do not need it and shouldn't have to
install such a large library to get yt.

A mechanism to support soft-dependencies is implemented in the
``yt/utilities/on_demand_imports.py`` module. Existing soft dependencies are
listed in a ``optional_dependencies`` dictionary using package names for keys.
The dictionary is then processed to generate drop-in replacement for the actual
packages, which behave as normal instances of the packages they mirror:

.. code-block:: python

    from yt.utilities.on_demand_imports import h5py


    def process(file):
        with h5py.File(file, mode="r") as fh:
            # perform interesting operations
            ...

In case the package is missing, an ``ImportError`` will be raised at runtime if
and only if the code using it is run. In the example above, this means that as
long as the ``process`` function is not run, no error will be raised. An
implication is that soft dependencies cannot be used inconditionally in the
global scope of a module, because it is run at yt's importtime, but should
rather be restricted to lower level scopes, inside of functions or properly
guarded conditional blocks.

Note that, mocking an external package using the ``OnDemandImport`` class also
delays any actual import to the first time it is used (this is how we can afford
to fail only at the conditions we have to). A direct implication is that objects
defined as such can safely be imported at the top level of a module, where the
bulk of import statements live, without affecting import durations in any
significant way.
