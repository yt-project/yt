"""
Integrator classes to deal with interpolation and integration of input spectral
bins.  Currently only supports Cloudy-style data.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2007-2011 Matthew Turk.  All Rights Reserved.

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

import numpy as np

from yt.funcs import *

from yt.data_objects.field_info_container import add_field
from yt.utilities.linear_interpolators import \
    UnilinearFieldInterpolator, \
    BilinearFieldInterpolator, \
    TrilinearFieldInterpolator

class SpectralFrequencyIntegrator(object):
    def __init__(self, table, field_names,
                 bounds, ev_bounds):
        """
        From a table, interpolate over field_names to get resultant luminosity.
        Table must be of the style such that it is ordered by
        ``[field_names[0], field_names[1], ev]``
        """
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
                        units=r"\rm{ergs}\/\rm{cm}^{-3}\/\rm{s}^{-1}",
                        projected_units=r"\rm{ergs}\/\rm{cm}^{-2}\/\rm{s}^{-1}")
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

