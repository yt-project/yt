Making Simple Plots
===================

Single Width Plots
------------------

Making plots can be considered one of the simplest things to do in yt, and at
its most primitive you can simply declare which plots to make.

The following recipe opens the parameter file, creates a 'collection' of plots,
adds a few slices, and then saves it at a fixed width.
(``cookbook_single_width_plot.py``)

.. literalinclude:: ../../../examples/cookbook_single_width_plot.py
   :language: python
   :linenos:

Multiple Fields Single Width
----------------------------

The 'save' command encodes the name of the field into the filename.  It is more
efficient to 'switch' a field than to add a new slice with a different field.
(``cookbook_multiple_fields_single_width_plot.py``)


.. literalinclude:: ../../../examples/cookbook_multiple_fields_single_width_plot.py
   :language: python
   :linenos:

Single Field Multiple Widths
----------------------------

We can zoom in on our slice very easily, and we can define it to do that and
vary units, too, thus ensuring we have a consistent set of plots.
(``cookbook_multiple_widths_plot.py``)

.. note::
   We do some fancy footwork here with the creation of *my_pairs* but it is an
   idiom that can be applied elsewhere.

.. literalinclude:: ../../../examples/cookbook_multiple_widths_plot.py
   :language: python
   :linenos:

Multiple Fields Multiple Widths
-------------------------------

Because of the way the slices are created and the fields handled, we set our
outer loop to be over the field and the inner loop to be over the widths.
(``cookbook_multiple_fields_multiple_widths_plot.py``)

.. literalinclude:: ../../../examples/cookbook_multiple_fields_multiple_widths_plot.py
   :language: python
   :linenos:

Awesome Linked Plots Saving
---------------------------

This is an idiom with which one can link several plots in widths, naming, and
colorbars.
(``cookbook_linked_plot_save.py``)

.. literalinclude:: ../../../examples/cookbook_linked_plot_save.py
   :language: python
   :linenos:

