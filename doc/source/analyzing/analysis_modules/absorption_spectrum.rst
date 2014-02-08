.. _absorption_spectrum:

Making an Absorption Spectrum
=============================
.. sectionauthor:: Britton Smith <brittonsmith@gmail.com>

Absorption line spectra, such as shown below, can be made with data created by the 
(:ref:`light-ray-generator`).  For each element of the ray, column densities are 
calculated multiplying the number density within a grid cell with the path length 
of the ray through the cell.  Line profiles are generated using a voigt profile based 
on the temperature field.  The lines are then shifted according to the redshift 
recorded by the light ray tool and (optionally) the line of sight peculiar velocity.  
Inclusion of the peculiar velocity requires setting **get_los_velocity** to True in 
the call to :meth:`make_light_ray`.

The spectrum generator will output a file containing the wavelength and normalized flux.  
It will also output a text file listing all important lines.

.. image:: _images/spectrum_full.png
   :width: 500

An absorption spectrum for the wavelength range from 900 to 1800 Angstroms made with 
a light ray extending from z = 0 to z = 0.4.

.. image:: _images/spectrum_zoom.png
   :width: 500

A zoom-in of the above spectrum.

Creating an Absorption Spectrum
-------------------------------

To instantiate an AbsorptionSpectrum object, the arguments required are the minimum and 
maximum wavelengths, and the number of wavelength bins.

.. code-block:: python

  from yt.analysis_modules.api import AbsorptionSpectrum

  sp = AbsorptionSpectrum(900.0, 1800.0, 10000)

Adding Features to the Spectrum
-------------------------------

Absorption lines and continuum features can then be added to the spectrum.  To add a 
line, you must know some properties of the line: the rest wavelength, f-value, gamma value, 
and the atomic mass in amu of the atom.  Below, we will add the H Lyman-alpha line.

.. code-block:: python
  
  my_label = 'HI Lya'
  field = 'HI_NumberDensity'
  wavelength = 1215.6700 # Angstroms
  f_value = 4.164E-01
  gamma = 6.265e+08
  mass = 1.00794
  
  sp.add_line(my_label, field, wavelength, f_value, gamma, mass, label_threshold=1.e10)

In the above example, the *field* argument tells the spectrum generator which field from the 
ray data to use to calculate the column density.  The **label_threshold** keyword tells the 
spectrum generator to add all lines above a column density of 10 :superscript:`10` 
cm :superscript:`-2` to the text line list.  If None is provided, as is the default, no 
lines of this type will be added to the text list.

Continuum features who optical depths follow a power law can also be added.  Below, we will add 
H Lyman continuum.

.. code-block:: python

  my_label = 'HI Lya'
  field = 'HI_NumberDensity'
  wavelength = 912.323660 # Angstroms
  normalization = 1.6e17
  index = 3.0
  
  sp.add_continuum(my_label, field, wavelength, normalization, index)

Making the Spectrum
-------------------

Once all the lines and continuum are added, it is time to make a spectrum out of 
some light ray data.

.. code-block:: python

  wavelength, flux = sp.make_spectrum('lightray.h5', output_file='spectrum.fits', 
                                      line_list_file='lines.txt',
                                      use_peculiar_velocity=True)

A spectrum will be made using the specified ray data and the wavelength and flux arrays 
will also be returned.  If **use_peculiar_velocity** is set to False, the lines will only 
be shifted according to the redshift.

Three output file formats are supported for writing out the spectrum: fits, hdf5, and ascii.  
The file format used is based on the extension provided in the **output_file** keyword: '.fits' 
for a fits file, '.h5' for an hdf5 file, and anything else for an ascii file.

.. note:: To write out a fits file, you must install the `pyfits <http://www.stsci.edu/resources/software_hardware/pyfits>`_ module.

What can I do with this?
------------------------

Try :ref:`quick_start_fitting`
