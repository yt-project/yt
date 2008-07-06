The Plugin File
===============

The plugin file is a means of modifying the available fields, quantities, data
objects and so on without modifying the source code of yt.

The following configuration parameters need to be set in the ``~/.yt/config``
file in order to enable the usage of a plugin file:

.. code-block:: python

   [lagos] 

   loadfieldplugins: True
   pluginfilename: my_plugins.py

You can call your plugin file whatever you like, and after the imports inside
the :mod:`lagos` module are completed, it is executed in that namespace.

The code in this file can thus add fields, add derived quantities, add
datatypes, and on and on.  For example, if I created a plugin file containing:

.. code-block:: python

   def _myfunc(field, data):
       return na.random.random(data["Density"].shape)
   add_field("SomeQuantity", function=_myfunc)

then all of my data objects would have access to the field "SomeQuantity"
despite its lack of use.
