.. _xray_emission_fields:

X-ray Emission Fields
=====================
.. sectionauthor:: Britton Smith <brittonsmith@gmail.com>, John ZuHone <jzuhone@gmail.com>

.. note::

  If you came here trying to figure out how to create simulated X-ray photons and observations,
  you should go `here <photon_simulator.html>`_ instead.

This functionality provides the ability to create metallicity-dependent 
X-ray luminosity, emissivity, and photon emissivity fields for a given
photon energy range.  This works by interpolating from emission tables 
created from the photoionization code `Cloudy <http://nublado.org/>`_ or
the collisional ionization database `AtomDB <http://www.atomdb.org>`_. If
you installed yt with the install script, these data files should be located in
the *data* directory inside the installation directory, or can be downloaded
from `<http://yt-project.org/data>`_. Emission fields can be made for any
interval between 0.1 keV and 100 keV.

Adding Emission Fields
----------------------

Fields will be created for luminosity :math:`{\rm (erg~s^{-1})}`, emissivity :math:`{\rm (erg~s^{-1}~cm^{-3})}`,
and photon emissivity :math:`{\rm (photons~s^{-1}~cm^{-3})}`.  The only required arguments are the
dataset object, and the minimum and maximum energies of the energy band.

.. code-block:: python

  import yt
  from yt.analysis_modules.spectral_integrator.api import \
       add_xray_emissivity_field

  xray_fields = add_xray_emissivity_field(0.5, 7.0)

Additional keyword arguments are:

 * **filename** (*string*): Path to data file containing emissivity values. If None,
   a file called "cloudy_emissivity.h5" is used, for photoionized plasmas. A second
   option, for collisionally ionized plasmas, is in the file "apec_emissivity.h5",
   available at http://yt-project.org/data. These files contain emissivity tables
   for primordial elements and for metals at solar metallicity for the energy range
   0.1 to 100 keV. Default: None.

 * **with_metals** (*bool*): If True, use the metallicity field to add the 
   contribution from metals.  If False, only the emission from H/He is 
   considered.  Default: True.

 * **constant_metallicity** (*float*): If specified, assume a constant 
   metallicity for the emission from metals.  The *with_metals* keyword 
   must be set to False to use this.  Default: None.

The resulting fields can be used like all normal fields. The function will return the names of
the created fields in a Python list.

.. code-block:: python

  import yt
  from yt.analysis_modules.spectral_integrator.api import \
       add_xray_emissivity_field

  xray_fields = add_xray_emissivity_field(0.5, 7.0, filename="apec_emissivity.h5")

  ds = yt.load("enzo_tiny_cosmology/DD0046/DD0046")
  plot = yt.SlicePlot(ds, 'x', 'xray_luminosity_0.5_7.0_keV')
  plot.save()
  plot = yt.ProjectionPlot(ds, 'x', 'xray_emissivity_0.5_7.0_keV')
  plot.save()
  plot = yt.ProjectionPlot(ds, 'x', 'xray_photon_emissivity_0.5_7.0_keV')
  plot.save()

.. warning::

  The X-ray fields depend on the number density of hydrogen atoms, in the yt field
  ``H_number_density``. If this field is not defined (either in the dataset or by the user),
  the primordial hydrogen mass fraction (X = 0.76) will be used to construct it.
