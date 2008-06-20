Extracting Data
===============

Oftentimes you will want to export data from yt to plot in another plotting
package, or to compare against another set of data elsewhere.  (Although, note
should be made that ever effort has been made to make yt into a
publication-quality plotter.)

Data objects all have the :meth:`write_out` method, which accepts a field
list and a filename.  This will dump to text.

However, if you have `PyTables <http://www.pytables.org/>`_ installed, you can
also export things in a different way.  

The following code snippet (``cookbook_extract_data.py``) will generate an
HDF5 file with the "Density" array from our smoothed cube at the root node, and
called "my_cube_density".  With PyTables you can nest these, set attributes,
and access the data through any application that can read HDF5.

.. literalinclude:: ../../../examples/cookbook_extract_data.py
   :language: python
   :linenos:
