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
import h5py

import math, itertools

YEAR = 3.155693e7

class StarFormationRate(object):
    def __init__(self, pf, data_source=None, star_mass=None,
            star_creation_time=None, volume=None, bins=300):
        self._pf = pf
        self._data_source = data_source
        self.star_mass = star_mass
        self.star_creation_time = star_creation_time
        self.volume = volume
        self.bin_count = bins
        # Check to make sure we have the right set of informations.
        if data_source is None:
            if self.star_mass is None or self.star_creation_time is None or \
            self.volume is None:
                mylog.error(
                """
                If data_source is not provided, all of these paramters need to be set:
                star_mass (array, Msun),
                star_creation_time (array, code units),
                volume (float, Mpc**3).
                """)
                return None
            self.mode = 'provided'
        else:
            self.mode = 'data_source'
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
        if self.mode == 'data_source':
            ct = self._data_source["creation_time"]
            ct_stars = ct[ct > 0]
            mass_stars = self._data_source["ParticleMassMsun"][ct > 0]
        elif self.mode == 'provided':
            ct_stars = self.star_creation_time
            mass_stars = self.star_mass
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
        if self.mode == 'data_source':
            vol = self._data_source.volume('mpc')
        elif self.mode == 'provided':
            vol = self.volume
        tc = self._pf["Time"]
        # Use the center of the time_bin, not the left edge.
        for i, time in enumerate((self.time_bins[1:] + self.time_bins[:-1])/2.):
            line = "%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\n" % \
            (time * tc / YEAR, # Time
            (self.time_now - time * tc)/YEAR, # Lookback time
            self.cosm.ComputeRedshiftFromTime(time * tc), # Redshift
            self.mass_bins[i] / (self.time_bins_dt[i] * tc / YEAR), # Msol/yr
            self.mass_bins[i] / (self.time_bins_dt[i] * tc / YEAR) / vol, # Msol/yr/vol
            self.mass_bins[i], # Msol in bin
            self.cum_mass_bins[i]) # cumulative
            fp.write(line)
        fp.close()

CHABRIER = {
"Z0001" : "bc2003_hr_m22_chab_ssp.ised.h5", #/* 0.5% */
"Z0004" : "bc2003_hr_m32_chab_ssp.ised.h5", #/* 2% */
"Z004" : "bc2003_hr_m42_chab_ssp.ised.h5", #/* 20% */
"Z008" : "bc2003_hr_m52_chab_ssp.ised.h5", #/* 40% */
"Z02" : "bc2003_hr_m62_chab_ssp.ised.h5", #/* solar; 0.02 */
"Z05" : "bc2003_hr_m72_chab_ssp.ised.h5" #/* 250% */
}

SALPETER = {
"Z0001" : "bc2003_hr_m22_salp_ssp.ised.h5", #/* 0.5% */
"Z0004" : "bc2003_hr_m32_salp_ssp.ised.h5", #/* 2% */
"Z004" : "bc2003_hr_m42_salp_ssp.ised.h5", #/* 20% */
"Z008" : "bc2003_hr_m52_salp_ssp.ised.h5", #/* 40% */
"Z02" : "bc2003_hr_m62_salp_ssp.ised.h5", #/* solar; 0.02 */
"Z05" : "bc2003_hr_m72_salp_ssp.ised.h5" #/* 250% */
}

Zsun = 0.02

#/* dividing line of metallicity; linear in log(Z/Zsun) */
METAL1 = 0.01  # /* in units of Z/Zsun */
METAL2 = 0.0632
METAL3 = 0.2828
METAL4 = 0.6325
METAL5 = 1.5811
METALS = na.array([METAL1, METAL2, METAL3, METAL4, METAL5])

# Translate METALS array digitize to the table dicts
MtoD = na.array(["Z0001", "Z0004", "Z004", "Z008", "Z02",  "Z05"])

class BuildSED(object):
    def __init__(self, pf, data_source=None, star_mass=None,
            star_creation_time=None, star_metalicity_fraction=None, bins=300,
            bcdir="", model="chabrier"):
        self._pf = pf
        self._data_source = data_source
        self.star_mass = star_mass
        self.star_creation_time = star_creation_time
        self.star_metal = star_metalicity_fraction
        self.bin_count = bins
        self.bcdir = bcdir
        
        if model == "chabrier":
            self.model = CHABRIER
        elif model == "salpeter":
            self.model = SALPETER
        
        # Find the time right now.
        self.time_now = self.cosm.ComputeTimeFromRedshift(
            self._pf["CosmologyCurrentRedshift"]) # seconds

    def read_bclib(self):
        """
        Read in the age and wavelength bins, and the flux bins for each
        metalicity.
        """
        self.flux = {}
        for file in self.model:
            fname = self.bcdir + self.model[file]
            fp = h5py.File(fname, 'r')
            self.age = fp["agebins"][:] # 1D floats
            self.wavelength = fp["wavebins"][:] # 1D floats
            self.flux[file] = fp["flam"][:,:] # 2D floats
            fp.close()
    
    def calculate_spectrum(self):
        """
        For the set of stars, calculate the collective SED.
        """
        # Initialize values
        self.final_spec = na.zeros(self.wavelength.size, dtype='float64')
        self.outflux = na.zeros(self.wavelength.size, dtype='float64')
    
        count = 0
        totmass = 0.
        tot_tform = 0.
        tot_metal = 0.
        zform_ave = 0.
        
        # Age of star in years.
        dt = (self.time_now - self.star_creation_time * self._pf['Time']) / YEAR
        # Figure out which METALS bin the star goes into.
        Mindex = na.digitize(dt, METALS)
        # Replace the indices with strings.
        Mname = na.choose(Mindex, MtoD)
        # Figure out which age bin this star goes into.
        Aindex = na.digitize(dt, self.age)
        # Ratios used for the interpolation.
        ratio1 = (dt - self.age[Aindex-1]) / (self.age[Aindex] - self.age[Aindex-1])
        ratio2 = (self.age[Aindex] - dt) / (self.age[Aindex] - self.age[Aindex-1])
        for star in iterools.izip(dt, Mname, Aindex, ratio1, ratio2):
            flux = self.flux[star[1]][:,star[2]]
    
    def interpolate(self, dt, flux):
        #/* first, find out the neighboring age steps */
        while(age[i]<dt) i++;
        
        if(i==0)
        {  
          for(j=0;j<WAVEBIN;j++) outflux[j]=flux[j][0]*dt/age[0];
          printf("negative dt?\n");
        }
        else
        {
          for(j=0;j<WAVEBIN;j++)
        {
          ratio1=(dt-age[i-1])/(age[i]-age[i-1]);
          ratio2=(age[i]-dt)/(age[i]-age[i-1]);
         
          /* interpolate in log(flux), linear in time */
          outlog = ratio2*log10(flux[j][i-1]) + ratio1*log10(flux[j][i]); 
          
          outflux[j]=pow(10,outlog);
