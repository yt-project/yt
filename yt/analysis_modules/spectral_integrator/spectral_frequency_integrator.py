"""
Integrator classes to deal with interpolation and integration of input spectral
bins.  Currently only supports Cloudy-style data.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Author: Britton Smith <brittons@origins.colorado.edu>
Affiliation: Michigan State University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2007-2012 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from exceptions import IOError
import h5py
import numpy as np
import os

from yt.funcs import *

from yt.data_objects.field_info_container import add_field
from yt.utilities.exceptions import YTException
from yt.utilities.linear_interpolators import \
    BilinearFieldInterpolator

xray_data_version = 1
    
class SpectralFrequencyIntegrator(object):
    """
    From a table, interpolate over field_names to get resultant luminosity.
    Table must be of the style such that it is ordered by
    ``[field_names[0], field_names[1], ev]``
    """
    def __init__(self, table, field_names,
                 bounds, ev_bounds):
        self.table = table
        self.field_names = field_names

        self.bounds = bounds
        self.ev_bounds = ev_bounds
        self.ev_vals = np.logspace(ev_bounds[0], ev_bounds[1], table.shape[-1])
        
    def _get_interpolator(self, ev_min, ev_max):
        """
        Integrates from ev_min to ev_max and returns an interpolator.
        """
        e_is, e_ie = np.digitize([ev_min, ev_max], self.ev_vals)
        bin_table = np.trapz(self.table[...,e_is-1:e_ie],
                             2.41799e17*
            (self.ev_vals[e_is:e_ie+1]-self.ev_vals[e_is-1:e_is]),
                             axis=-1)
        bin_table = np.log10(bin_table.clip(1e-80,bin_table.max()))
        return BilinearFieldInterpolator(
            bin_table, self.bounds, self.field_names[:],
            truncate=True)



    def add_frequency_bin_field(self, ev_min, ev_max):
        """
        Add a new field to the FieldInfoContainer, which is an
        integrated bin from *ev_min* to *ev_max*.
        
        Returns the name of the new field.
        """
        interp = self._get_interpolator(ev_min, ev_max)
        name = "XRay_%s_%s" % (ev_min, ev_max)
        def frequency_bin_field(field, data):
            dd = {'H_NumberDensity' : np.log10(data["H_NumberDensity"]),
                  'Temperature'   : np.log10(data["Temperature"])}
            return 10**interp(dd)
        add_field(name, function=frequency_bin_field,
                        projection_conversion="cm",
                        units=r"\rm{ergs}\ \rm{cm}^{-3}\ \rm{s}^{-1}",
                        projected_units=r"\rm{ergs}\ \rm{cm}^{-2}\ \rm{s}^{-1}")
        return name

def create_table_from_textfiles(pattern, rho_spec, e_spec, T_spec):
    """
    This accepts a CLOUDY text file of emissivities and constructs an
    interpolation table for spectral integration.
    """
    rho_n_bins, rho_min, rho_max = rho_spec
    e_n_bins, e_min, e_max = e_spec
    T_n_bins, T_min, T_max = T_spec
    # The second one is the fast-varying one
    rho_is, e_is = np.mgrid[0:rho_n_bins,0:e_n_bins]
    table = np.zeros((rho_n_bins, T_n_bins, e_n_bins), dtype='float64')
    mylog.info("Parsing Cloudy files")
    for i,ri,ei in zip(range(rho_n_bins*e_n_bins), rho_is.ravel(), e_is.ravel()):
        table[ri,:,ei] = [float(l.split()[-1]) for l in open(pattern%(i+1)) if l[0] != "#"]
    return table

class EnergyBoundsException(YTException):
    def __init__(self, lower, upper):
        self.lower = lower
        self.upper = upper

    def __str__(self):
        return "Energy bounds are %e to %e keV." % \
          (self.lower, self.upper)

class ObsoleteDataException(YTException):
    def __str__(self):
        return "X-ray emissivity data is out of data.\nDownload the latest data from http://yt-project.org/data/xray_emissivity.h5 and move it to %s." % \
          os.path.join(os.environ["YT_DEST"], "data", "xray_emissivity.h5")
          
class EmissivityIntegrator(object):
    r"""Initialize an EmissivityIntegrator object.

    Keyword Parameters
    ------------------
    filename: string
        Path to data file containing emissivity values.  If None,
        a file called xray_emissivity.h5 is used.  This file contains 
        emissivity tables for primordial elements and for metals at 
        solar metallicity for the energy range 0.1 to 100 keV.
        Default: None.
        
    """
    def __init__(self, filename=None):

        default_filename = False
        if filename is None:
            filename = os.path.join(os.environ["YT_DEST"], 
                                    "data", "xray_emissivity.h5")
            default_filename = True

        if not os.path.exists(filename):
            raise IOError("File does not exist: %s." % filename)
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
            setattr(self, field, in_file[field][:])
        in_file.close()

        E_diff = np.diff(self.log_E)
        self.E_bins = \
                  np.power(10, np.concatenate([self.log_E[:-1] - 0.5 * E_diff,
                                               [self.log_E[-1] - 0.5 * E_diff[-1],
                                                self.log_E[-1] + 0.5 * E_diff[-1]]]))
        self.dnu = 2.41799e17 * np.diff(self.E_bins)

    def _get_interpolator(self, data, e_min, e_max):
        r"""Create an interpolator for total emissivity in a 
        given energy range.

        Parameters
        ----------
        e_min: float
            the minimum energy in keV for the energy band.
        e_min: float
            the maximum energy in keV for the energy band.

        """
        if (e_min - self.E_bins[0]) / e_min < -1e-3 or \
          (e_max - self.E_bins[-1]) / e_max > 1e-3:
            raise EnergyBoundsException(np.power(10, self.E_bins[0]),
                                        np.power(10, self.E_bins[-1]))
        e_is, e_ie = np.digitize([e_min, e_max], self.E_bins)
        e_is = np.clip(e_is - 1, 0, self.E_bins.size - 1)
        e_ie = np.clip(e_ie, 0, self.E_bins.size - 1)

        my_dnu = np.copy(self.dnu[e_is: e_ie])
        # clip edge bins if the requested range is smaller
        my_dnu[0] -= e_min - self.E_bins[e_is]
        my_dnu[-1] -= self.E_bins[e_ie] - e_max

        interp_data = (data[..., e_is:e_ie] * my_dnu).sum(axis=-1)
        return BilinearFieldInterpolator(np.log10(interp_data),
                                         [self.log_nH[0], self.log_nH[-1],
                                          self.log_T[0],  self.log_T[-1]],
                                         ["log_nH", "log_T"], truncate=True)

def add_xray_emissivity_field(e_min, e_max, filename=None,
                              with_metals=True,
                              constant_metallicity=None):
    r"""Create an X-ray emissivity field for a given energy range.

    Parameters
    ----------
    e_min: float
        the minimum energy in keV for the energy band.
    e_min: float
        the maximum energy in keV for the energy band.

    Keyword Parameters
    ------------------
    filename: string
        Path to data file containing emissivity values.  If None,
        a file called xray_emissivity.h5 is used.  This file contains 
        emissivity tables for primordial elements and for metals at 
        solar metallicity for the energy range 0.1 to 100 keV.
        Default: None.
    with_metals: bool
        If True, use the metallicity field to add the contribution from 
        metals.  If False, only the emission from H/He is considered.
        Default: True.
    constant_metallicity: float
        If specified, assume a constant metallicity for the emission 
        from metals.  The *with_metals* keyword must be set to False 
        to use this.
        Default: None.

    This will create a field named "Xray_Emissivity_{e_min}_{e_max}keV".
    The units of the field are erg s^-1 cm^-3.

    Examples
    --------

    >>> from yt.mods import *
    >>> from yt.analysis_modules.spectral_integrator.api import *
    >>> add_xray_emissivity_field(0.5, 2)
    >>> pf = load(dataset)
    >>> p = ProjectionPlot(pf, 'x', "Xray_Emissivity_0.5_2keV")
    >>> p.save()

    """

    my_si = EmissivityIntegrator(filename=filename)

    em_0 = my_si._get_interpolator(my_si.emissivity_primordial, e_min, e_max)
    em_Z = None
    if with_metals or constant_metallicity is not None:
        em_Z = my_si._get_interpolator(my_si.emissivity_metals, e_min, e_max)

    def _emissivity_field(field, data):
        dd = {"log_nH" : np.log10(data["H_NumberDensity"]),
              "log_T"   : np.log10(data["Temperature"])}

        my_emissivity = np.power(10, em_0(dd))
        if em_Z is not None:
            if with_metals:
                my_Z = data["Metallicity"]
            elif constant_metallicity is not None:
                my_Z = constant_metallicity
            my_emissivity += my_Z * np.power(10, em_Z(dd))

        return data["H_NumberDensity"]**2 * my_emissivity

    field_name = "Xray_Emissivity_%s_%skeV" % (e_min, e_max)
    add_field(field_name, function=_emissivity_field,
              projection_conversion="cm",
              display_name=r"\epsilon_{X}\/(%s-%s\/keV)" % (e_min, e_max),
              units=r"\rm{erg}\/\rm{cm}^{-3}\/\rm{s}^{-1}",
              projected_units=r"\rm{erg}\/\rm{cm}^{-2}\/\rm{s}^{-1}")
    return field_name

def add_xray_luminosity_field(e_min, e_max, filename=None,
                              with_metals=True,
                              constant_metallicity=None):
    r"""Create an X-ray luminosity field for a given energy range.

    Parameters
    ----------
    e_min: float
        the minimum energy in keV for the energy band.
    e_min: float
        the maximum energy in keV for the energy band.

    Keyword Parameters
    ------------------
    filename: string
        Path to data file containing emissivity values.  If None,
        a file called xray_emissivity.h5 is used.  This file contains 
        emissivity tables for primordial elements and for metals at 
        solar metallicity for the energy range 0.1 to 100 keV.
        Default: None.
    with_metals: bool
        If True, use the metallicity field to add the contribution from 
        metals.  If False, only the emission from H/He is considered.
        Default: True.
    constant_metallicity: float
        If specified, assume a constant metallicity for the emission 
        from metals.  The *with_metals* keyword must be set to False 
        to use this.
        Default: None.

    This will create a field named "Xray_Luminosity_{e_min}_{e_max}keV".
    The units of the field are erg s^-1.

    Examples
    --------

    >>> from yt.mods import *
    >>> from yt.analysis_modules.spectral_integrator.api import *
    >>> add_xray_luminosity_field(0.5, 2)
    >>> pf = load(dataset)
    >>> sp = pf.h.sphere('max', (2., 'mpc'))
    >>> print sp.quantities['TotalQuantity']('Xray_Luminosity_0.5_2keV')
    
    """

    em_field = add_xray_emissivity_field(e_min, e_max, filename=filename,
                                         with_metals=with_metals,
                                         constant_metallicity=constant_metallicity)

    def _luminosity_field(field, data):
        return data[em_field] * data["CellVolume"]
    field_name = "Xray_Luminosity_%s_%skeV" % (e_min, e_max)
    add_field(field_name, function=_luminosity_field,
              display_name=r"\rm{L}_{X}\/(%s-%s\/keV)" % (e_min, e_max),
              units=r"\rm{erg}\/\rm{s}^{-1}")
    return field_name

def add_xray_photon_emissivity_field(e_min, e_max, filename=None,
                                     with_metals=True,
                                     constant_metallicity=None):
    r"""Create an X-ray photon emissivity field for a given energy range.

    Parameters
    ----------
    e_min: float
        the minimum energy in keV for the energy band.
    e_min: float
        the maximum energy in keV for the energy band.

    Keyword Parameters
    ------------------
    filename: string
        Path to data file containing emissivity values.  If None,
        a file called xray_emissivity.h5 is used.  This file contains 
        emissivity tables for primordial elements and for metals at 
        solar metallicity for the energy range 0.1 to 100 keV.
        Default: None.
    with_metals: bool
        If True, use the metallicity field to add the contribution from 
        metals.  If False, only the emission from H/He is considered.
        Default: True.
    constant_metallicity: float
        If specified, assume a constant metallicity for the emission 
        from metals.  The *with_metals* keyword must be set to False 
        to use this.
        Default: None.

    This will create a field named "Xray_Photon_Emissivity_{e_min}_{e_max}keV".
    The units of the field are photons s^-1 cm^-3.

    Examples
    --------

    >>> from yt.mods import *
    >>> from yt.analysis_modules.spectral_integrator.api import *
    >>> add_xray_emissivity_field(0.5, 2)
    >>> pf = load(dataset)
    >>> p = ProjectionPlot(pf, 'x', "Xray_Emissivity_0.5_2keV")
    >>> p.save()

    """

    my_si = EmissivityIntegrator(filename=filename)
    energy_erg = np.power(10, my_si.log_E) * 1.60217646e-9

    em_0 = my_si._get_interpolator((my_si.emissivity_primordial[..., :] / energy_erg),
                                   e_min, e_max)
    em_Z = None
    if with_metals or constant_metallicity is not None:
        em_Z = my_si._get_interpolator((my_si.emissivity_metals[..., :] / energy_erg),
                                       e_min, e_max)

    def _emissivity_field(field, data):
        dd = {"log_nH" : np.log10(data["H_NumberDensity"]),
              "log_T"   : np.log10(data["Temperature"])}

        my_emissivity = np.power(10, em_0(dd))
        if em_Z is not None:
            if with_metals:
                my_Z = data["Metallicity"]
            elif constant_metallicity is not None:
                my_Z = constant_metallicity
            my_emissivity += my_Z * np.power(10, em_Z(dd))

        return data["H_NumberDensity"]**2 * my_emissivity

    field_name = "Xray_Photon_Emissivity_%s_%skeV" % (e_min, e_max)
    add_field(field_name, function=_emissivity_field,
              projection_conversion="cm",
              display_name=r"\epsilon_{X}\/(%s-%s\/keV)" % (e_min, e_max),
              units=r"\rm{photons}\/\rm{cm}^{-3}\/\rm{s}^{-1}",
              projected_units=r"\rm{photons}\/\rm{cm}^{-2}\/\rm{s}^{-1}")
    return field_name
