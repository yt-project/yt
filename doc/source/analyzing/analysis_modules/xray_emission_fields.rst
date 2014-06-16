.. _xray_emission_fields:

X-ray Emission Fields
=====================
.. sectionauthor:: Britton Smith <brittonsmith@gmail.com>

This functionality provides the ability to create metallicity-dependent 
X-ray luminosity, emissivity, and photo emissivity fields for a given 
photon energy range.  This works by interpolating from emission tables 
created with the photoionization code, `Cloudy <http://nublado.org/>`_.  
If you installed yt with the install script, the data should be located in 
the *data* directory inside the installation directory.  Emission fields can 
be made for any interval between 0.1 keV and 100 keV.

Adding Emission Fields
----------------------

Fields can be created for luminosity (erg/s), emissivity (erg/s/cm^3), 
and photon emissivity (photons/s/cm^3).  The only required arguments are 
the minimum and maximum energies.

.. code-block:: python

  from yt.mods import *
  from yt.analysis_modules.spectral_integrator.api import \
       add_xray_luminosity_field, \
       add_xray_emissivity_field, \
       add_xray_photon_emissivity_field

  add_xray_luminosity_field(0.5, 7)
  add_xray_emissivity_field(0.5, 7)
  add_xray_photon_emissivity_field(0.5, 7)

Additional keyword arguments are:

 * **filename**  (*string*): Path to data file containing emissivity 
   values.  If None, a file called xray_emissivity.h5 is used.  This file 
   contains emissivity tables for primordial elements and for metals at 
   solar metallicity for the energy range 0.1 to 100 keV.  Default: None.

 * **with_metals** (*bool*): If True, use the metallicity field to add the 
   contribution from metals.  If False, only the emission from H/He is 
   considered.  Default: True.

 * **constant_metallicity** (*float*): If specified, assume a constant 
   metallicity for the emission from metals.  The *with_metals* keyword 
   must be set to False to use this.  Default: None.

The resulting fields can be used like all normal fields.

.. python-script::

  from yt.mods import *
  from yt.analysis_modules.spectral_integrator.api import \
       add_xray_luminosity_field, \
       add_xray_emissivity_field, \
       add_xray_photon_emissivity_field

  add_xray_luminosity_field(0.5, 7)
  add_xray_emissivity_field(0.5, 7)
  add_xray_photon_emissivity_field(0.5, 7)

  ds = load("enzo_tiny_cosmology/DD0046/DD0046")
  plot = SlicePlot(ds, 'x', 'Xray_Luminosity_0.5_7keV')
  plot.save()
  plot = ProjectionPlot(ds, 'x', 'Xray_Emissivity_0.5_7keV')
  plot.save()
  plot = ProjectionPlot(ds, 'x', 'Xray_Photon_Emissivity_0.5_7keV')
  plot.save()
