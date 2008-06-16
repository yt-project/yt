Extracting Data
===============

Oftentimes you will want to export data from yt to plot in another plotting
package, or to compare against another set of data elsewhere.  (Although, note
should be made that ever effort has been made to make yt into a
publication-quality plotter.)

Data objects all have the :meth:`write_out` method, which accepts a field
list and a filename.  This will dump to text.

However, if you have `PyTables <http://www.pytables.org/>`_ installed, you can
also export things in a different way: ::

   >>> cube = pf.h.smoothed_covering_grid(3, left_edge=[0.0,0.0,0.0],
   ...                                       right_edge=[1.0,1.0,1.0],
   ...                                       dims=[512,512,512],
   ...                                       fields=["Density"])
   >>> f = tables.openFile("my_cube.h5","w")
   >>> f.createArray("/","my_cube_density", cube["Density"])
   >>> f.close()

This will generate an HDF5 file with the "Density" array from our smoothed cube
at the root node, and called "my_cube_density".  With PyTables you can nest
these, set attributes, and access the data through any application that can
read HDF5.
