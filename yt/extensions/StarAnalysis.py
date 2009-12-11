"""
StarAnalysis - Functions to analyze stars.

Author: Stephen Skory <sskory@physics.ucsd.edu>
Affiliation: UC San Diego / CASS
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008-2009 Stephen Skory (and others).  All Rights Reserved.

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

import yt.lagos as lagos
from yt.logger import lagosLogger as mylog
import numpy as na
import math

YEAR = 3.155693e7

class StarFormationRate(object):
    def __init__(self, pf, data_source, bins=300):
        self._pf = pf
        self._data_source = data_source
        self.bin_count = bins
        # Set up for time conversion.
        self.cosm = lagos.EnzoCosmology(HubbleConstantNow = 
             (100.0 * self._pf['CosmologyHubbleConstantNow']),
             OmegaMatterNow = self._pf['CosmologyOmegaMatterNow'],
             OmegaLambdaNow = self._pf['CosmologyOmegaLambdaNow'],
             InitialRedshift = self._pf['CosmologyInitialRedshift'])
        # Find the time right now.
        self.time_now = self.cosm.ComputeTimeFromRedshift(
            self._pf["CosmologyCurrentRedshift"]) # seconds
        # Build the distribution.
        self.build_dist()

    def build_dist(self):
        """
        Build the data for plotting.
        """
        # Pick out the stars.
        ct = self._data_source["creation_time"]
        ct_stars = ct[ct > 0]
        mass_stars = self._data_source["ParticleMassMsun"][ct > 0]
        # Find the oldest stars in units of code time.
        tmin= min(ct_stars)
        # Multiply the end to prevent numerical issues.
        self.time_bins = na.linspace(tmin*0.99, self._pf['InitialTime'],
            num = self.bin_count + 1)
        # Figure out which bins the stars go into.
        inds = na.digitize(ct_stars, self.time_bins) - 1
        # Sum up the stars created in each time bin.
        self.mass_bins = na.zeros(self.bin_count + 1, dtype='float64')
        for index in na.unique(inds):
            self.mass_bins[index] += sum(mass_stars[inds == index])
        # Calculate the cumulative mass sum over time by forward adding.
        self.cum_mass_bins = self.mass_bins.copy()
        for index in xrange(self.bin_count):
            self.cum_mass_bins[index+1] += self.cum_mass_bins[index]
        # We will want the time taken between bins.
        self.time_bins_dt = self.time_bins[1:] - self.time_bins[:-1]
    
    def write_out(self, name="StarAnalysis.out"):
        """
        Write out the star analysis to a text file *name*. The columns are in
        order:
        1) Time (yrs)
        2) Look-back time (yrs)
        3) Redshift
        4) Star formation rate in this bin per year (Msol/yr)
        5) Star formation rate in this bin per year per Mpc**3 (Msol/yr/Mpc**3)
        6) Stars formed in this time bin (Msol)
        7) Cumulative stars formed up to this time bin (Msol)
        """
        fp = open(name, "w")
        vol = self._data_source.volume('mpc')
        tc = self._pf["Time"]
        # Use the center of the time_bin, not the left edge.
        for i, time in enumerate((self.time_bins[1:] + self.time_bins[:-1])/2.):
            line = "%1.5e %1.5e %1.5e %1.5e %1.5e %1.5e %1.5e\n" % \
            (time * tc / YEAR, # Time
            (self.time_now - time * tc)/YEAR, # Lookback time
            self.cosm.ComputeRedshiftFromTime(time * tc), # Redshift
            self.mass_bins[i] / (self.time_bins_dt[i] * tc / YEAR), # Msol/yr
            self.mass_bins[i] / (self.time_bins_dt[i] * tc / YEAR) / vol, # Msol/yr/vol
            self.mass_bins[i], # Msol in bin
            self.cum_mass_bins[i]) # cumulative
            fp.write(line)
        fp.close()
        
            