Making Simple Plots
===================

Single Width Plots
------------------

Making plots can be considered one of the simplest things to do in yt, and at
its most primitive you can simply declare which plots to make.

The following recipe opens the parameter file, creates a 'collection' of plots,
adds a few slices, and then saves it at a fixed width.

.. code-block:: python

   pf = lagos.EnzoStaticOutput("my_data/my_data0001")
   pc = raven.PlotCollection(pf)
   pc.add_slice("Density",0)
   pc.add_slice("Density",1)
   pc.add_slice("Density",2)
   pc.set_width(100.0,'kpc')
   pc.save("my_data0001_100kpc")

Multiple Fields Single Width
----------------------------

The 'save' command encodes the name of the field into the filename.  It is more
efficient to 'switch' a field than to add a new slice with a different field.

.. code-block:: python

   fields = ["Density", "Temperature", "x-velocity"]
   pf = lagos.EnzoStaticOutput("my_data/my_data0001")
   pc = raven.PlotCollection(pf)
   pc.add_slice(fields[0],0)
   pc.add_slice(fields[0],1)
   pc.add_slice(fields[0],2)
   pc.set_width(100.0,'kpc')
   for field in fields:
       pc.switch_field(field)
       pc.save("my_data0001_100kpc")

Single Field Multiple Widths
----------------------------

We can zoom in on our slice very easily, and we can define it to do that and
vary units, too, thus ensuring we have a consistent set of plots.

.. note::
   We do some fancy footwork here with the creation of *my_pairs* but it is an
   idiom that can be applied elsewhere.

.. code-block:: python

   widths = [1000.0, 100.0, 10.0, 1.0]
   units = ['mpc','kpc','pc','au']
   my_pairs = [ (w,u) for u in units for w in widths ]

   pf = lagos.EnzoStaticOutput("my_data/my_data0001")
   pc = raven.PlotCollection(pf)
   pc.add_slice("Density",0)
   pc.add_slice("Density",1)
   pc.add_slice("Density",2)
   for w, u in my_pairs:
        pc.set_width(w,u)
        pc.save("my_data0001_%05i%s" % (w, u))

Multiple Fields Multiple Widths
-------------------------------

Because of the way the slices are created and the fields handled, we set our
outer loop to be over the field and the inner loop to be over the widths.

.. code-block:: python

   fields = ["Density", "Temperature", "x-velocity"]
   widths = [1000.0, 100.0, 10.0, 1.0]
   units = ['mpc','kpc','pc','au']
   my_pairs = [ (w,u) for u in units for w in widths ]

   pf = lagos.EnzoStaticOutput("my_data/my_data0001")
   pc = raven.PlotCollection(pf)
   pc.add_slice(fields[0],0)
   pc.add_slice(fields[0],1)
   pc.add_slice(fields[0],2)
   for field in fields:
       pc.switch_field(field)
       for w, u in my_pairs:
            pc.set_width(w,u)
            pc.save("my_data0001_%05i%s" % (w, u))


