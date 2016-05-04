"""
Integrator classes to deal with interpolation and integration of input spectral
bins.  Currently only supports Cloudy and APEC-style data.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.on_demand_imports import _h5py as h5py
import numpy as np
import os

from yt.funcs import \
     download_file, \
     mylog, \
     only_on_root

from yt.utilities.exceptions import YTFieldNotFound
from yt.utilities.exceptions import YTException
from yt.utilities.linear_interpolators import \
    UnilinearFieldInterpolator, BilinearFieldInterpolator
from yt.utilities.physical_constants import \
    hcgs, mp
from yt.units.yt_array import YTArray, YTQuantity
from yt.utilities.physical_ratios import \
    primordial_H_mass_fraction, erg_per_keV

xray_data_version = 1

def _get_data_file(data_file=None):
    if data_file is None:
        data_file = "cloudy_emissivity.h5"
    data_url = "http://yt-project.org/data"
    if "YT_DEST" in os.environ and \
      os.path.isdir(os.path.join(os.environ["YT_DEST"], "data")):
        data_dir = os.path.join(os.environ["YT_DEST"], "data")
    else:
        data_dir = "."
    data_path = os.path.join(data_dir, data_file)
    if not os.path.exists(data_path):
        mylog.info("Attempting to download supplementary data from %s to %s." % 
                   (data_url, data_dir))
        fn = download_file(os.path.join(data_url, data_file), data_path)
        if fn != data_path:
            raise RuntimeError("Failed to download supplementary data.")
    return data_path

class EnergyBoundsException(YTException):
    def __init__(self, lower, upper):
        self.lower = lower
        self.upper = upper

    def __str__(self):
        return "Energy bounds are %e to %e keV." % \
          (self.lower, self.upper)

class ObsoleteDataException(YTException):
    def __str__(self):
        return "X-ray emissivity data is out of date.\n" + \
               "Download the latest data from http://yt-project.org/data/cloudy_emissivity.h5 and move it to %s." % \
          os.path.join(os.environ["YT_DEST"], "data", "cloudy_emissivity.h5")
          
class EmissivityIntegrator(object):
    r"""Class for making X-ray emissivity fields with hdf5 data tables 
    from Cloudy.
    
    Initialize an EmissivityIntegrator object.

    Parameters
    ----------
    filename: string, default None
        Path to data file containing emissivity values.  If None,
        a file called "cloudy_emissivity.h5" is used, for photoionized
        plasmas. A second option, for collisionally ionized plasmas, is
        in the file "apec_emissivity.h5", available at http://yt-project.org/data.
        These files contain emissivity tables for primordial elements and
        for metals at solar metallicity for the energy range 0.1 to 100 keV.
        Default: None.
        
    """
    def __init__(self, filename=None):

        default_filename = False
        if filename is None:
            filename = _get_data_file()
            default_filename = True

        if not os.path.exists(filename):
            mylog.warning("File %s does not exist, will attempt to find it." % filename)
            filename = _get_data_file(data_file=filename)
        only_on_root(mylog.info, "Loading emissivity data from %s." % filename)
        in_file = h5py.File(filename, "r")
        if "info" in in_file.attrs:
            only_on_root(mylog.info, in_file.attrs["info"])
        if default_filename and \
          in_file.attrs["version"] < xray_data_version:
            raise ObsoleteDataException()
        else:
            only_on_root(mylog.info, "X-ray emissivity data version: %s." % \
                         in_file.attrs["version"])

        for field in ["emissivity_primordial", "emissivity_metals",
                      "log_nH", "log_T", "log_E"]:
            if field in in_file:
                setattr(self, field, in_file[field][:])
        in_file.close()

        E_diff = np.diff(self.log_E)
        self.E_bins = \
                  YTArray(np.power(10, np.concatenate([self.log_E[:-1] - 0.5 * E_diff,
                                                      [self.log_E[-1] - 0.5 * E_diff[-1],
                                                       self.log_E[-1] + 0.5 * E_diff[-1]]])),
                          "keV")
        self.dnu = (np.diff(self.E_bins)/hcgs).in_units("Hz")

    def get_interpolator(self, data, e_min, e_max):
        e_min = YTQuantity(e_min, "keV")
        e_max = YTQuantity(e_max, "keV")
        if (e_min - self.E_bins[0]) / e_min < -1e-3 or \
          (e_max - self.E_bins[-1]) / e_max > 1e-3:
            raise EnergyBoundsException(self.E_bins[0], self.E_bins[-1])
        e_is, e_ie = np.digitize([e_min, e_max], self.E_bins)
        e_is = np.clip(e_is - 1, 0, self.E_bins.size - 1)
        e_ie = np.clip(e_ie, 0, self.E_bins.size - 1)

        my_dnu = self.dnu[e_is: e_ie].copy()
        # clip edge bins if the requested range is smaller
        my_dnu[0] -= ((e_min - self.E_bins[e_is])/hcgs).in_units("Hz")
        my_dnu[-1] -= ((self.E_bins[e_ie] - e_max)/hcgs).in_units("Hz")

        interp_data = (data[..., e_is:e_ie] * my_dnu).sum(axis=-1)
        if len(data.shape) == 2:
            emiss = UnilinearFieldInterpolator(np.log10(interp_data),
                                               [self.log_T[0],  self.log_T[-1]],
                                               "log_T", truncate=True)
        else:
            emiss = BilinearFieldInterpolator(np.log10(interp_data),
                                              [self.log_nH[0], self.log_nH[-1],
                                               self.log_T[0],  self.log_T[-1]],
                                              ["log_nH", "log_T"], truncate=True)

        return emiss

def add_xray_emissivity_field(ds, e_min, e_max,
                              filename=None,
                              with_metals=True,
                              constant_metallicity=None):
    r"""Create X-ray emissivity fields for a given energy range.

    Parameters
    ----------
    e_min: float
        the minimum energy in keV for the energy band.
    e_min: float
        the maximum energy in keV for the energy band.
    filename: string, optional
        Path to data file containing emissivity values.  If None,
        a file called "cloudy_emissivity.h5" is used, for photoionized
        plasmas. A second option, for collisionally ionized plasmas, is
        in the file "apec_emissivity.h5", available at http://yt-project.org/data.
        These files contain emissivity tables for primordial elements and
        for metals at solar metallicity for the energy range 0.1 to 100 keV.
        Default: None.
    with_metals: bool, optional
        If True, use the metallicity field to add the contribution from 
        metals.  If False, only the emission from H/He is considered.
        Default: True.
    constant_metallicity: float, optional
        If specified, assume a constant metallicity for the emission 
        from metals.  The *with_metals* keyword must be set to False 
        to use this. It should be given in unit of solar metallicity.
        Default: None.

    This will create three fields:

    "xray_emissivity_{e_min}_{e_max}_keV" (erg s^-1 cm^-3)
    "xray_luminosity_{e_min}_{e_max}_keV" (erg s^-1)
    "xray_photon_emissivity_{e_min}_{e_max}_keV" (photons s^-1 cm^-3)

    Examples
    --------

    >>> from yt.mods import *
    >>> from yt.analysis_modules.spectral_integrator.api import *
    >>> ds = load(dataset)
    >>> add_xray_emissivity_field(ds, 0.5, 2)
    >>> p = ProjectionPlot(ds, 'x', "xray_emissivity_0.5_2_keV")
    >>> p.save()

    """

    if with_metals:
        try:
            ds._get_field_info("metal_density")
        except YTFieldNotFound:
            raise RuntimeError("Your dataset does not have a \"metal_density\" field! " +
                               "Perhaps you should specify a constant metallicity?")

    my_si = EmissivityIntegrator(filename=filename)

    em_0 = my_si.get_interpolator(my_si.emissivity_primordial, e_min, e_max)
    em_Z = None
    if with_metals or constant_metallicity is not None:
        em_Z = my_si.get_interpolator(my_si.emissivity_metals, e_min, e_max)

    energy_erg = np.power(10, my_si.log_E) * erg_per_keV
    emp_0 = my_si.get_interpolator((my_si.emissivity_primordial[..., :] / energy_erg),
                                   e_min, e_max)
    emp_Z = None
    if with_metals or constant_metallicity is not None:
        emp_Z = my_si.get_interpolator((my_si.emissivity_metals[..., :] / energy_erg),
                                       e_min, e_max)

    try:
        ds._get_field_info("H_number_density")
    except YTFieldNotFound:
        mylog.warning("Could not find a field for \"H_number_density\". Assuming primordial H " +
                      "mass fraction.")
        def _nh(field, data):
            return primordial_H_mass_fraction*data["gas","density"]/mp
        ds.add_field(("gas", "H_number_density"), function=_nh, units="cm**-3")

    def _emissivity_field(field, data):
        dd = {"log_nH" : np.log10(data["gas","H_number_density"]),
              "log_T"   : np.log10(data["gas","temperature"])}

        my_emissivity = np.power(10, em_0(dd))
        if em_Z is not None:
            if with_metals:
                my_Z = data["gas","metallicity"]
            elif constant_metallicity is not None:
                my_Z = constant_metallicity
            my_emissivity += my_Z * np.power(10, em_Z(dd))

        return data["gas","H_number_density"]**2 * \
            YTArray(my_emissivity, "erg*cm**3/s")

    emiss_name = "xray_emissivity_%s_%s_keV" % (e_min, e_max)
    ds.add_field(("gas", emiss_name), function=_emissivity_field,
                 display_name=r"\epsilon_{X} (%s-%s keV)" % (e_min, e_max),
                 units="erg/cm**3/s")

    def _luminosity_field(field, data):
        return data[emiss_name] * data["cell_volume"]

    lum_name = "xray_luminosity_%s_%s_keV" % (e_min, e_max)
    ds.add_field(("gas", lum_name), function=_luminosity_field,
                 display_name=r"\rm{L}_{X} (%s-%s keV)" % (e_min, e_max),
                 units="erg/s")

    def _photon_emissivity_field(field, data):
        dd = {"log_nH" : np.log10(data["gas","H_number_density"]),
              "log_T"   : np.log10(data["gas","temperature"])}

        my_emissivity = np.power(10, emp_0(dd))
        if emp_Z is not None:
            if with_metals:
                my_Z = data["gas","metallicity"]
            elif constant_metallicity is not None:
                my_Z = constant_metallicity
            my_emissivity += my_Z * np.power(10, emp_Z(dd))

        return data["gas","H_number_density"]**2 * \
            YTArray(my_emissivity, "photons*cm**3/s")

    phot_name = "xray_photon_emissivity_%s_%s_keV" % (e_min, e_max)
    ds.add_field(("gas", phot_name), function=_photon_emissivity_field,
                 display_name=r"\epsilon_{X} (%s-%s keV)" % (e_min, e_max),
                 units="photons/cm**3/s")

    return emiss_name, lum_name, phot_name
