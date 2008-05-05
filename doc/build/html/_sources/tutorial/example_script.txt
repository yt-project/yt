A Demonstrative Example Script
==============================

During the course of my work, I often have to run a selection of data analysis
tasks in order to keep tabs on my simulation as it runs.  I've included below a
script that I run on each output as it is created.

.. code-block:: python

   import sys
   import yt.lagos as lagos
   import yt.raven as raven
   
   fn = sys.argv[-1]
   pf = lagos.EnzoStaticOutput(fn)
   pf.h.print_stats()

   fields_to_plot = ["NumberDensity", "Density", "Temperature", "H2I_Fraction"]
   
   my_vals = [1000.0, 100.0, 10.0, 1.0]
   my_units = ['au','rsun']
   
   my_widths = []
   
   for unit in my_units:
       my_widths += [(v,units) for v in my_vals]
   
   pc = raven.PlotCollection(pf)
   for i in range(3): pc.add_slice("NumberDensity", i)
   
   for field in fields_to_plot:
       pc.switch_field(field)
       for w, u in my_widths:
           pc.set_width(w,u)
           pc.save("%s" % pf)

And that's it!
