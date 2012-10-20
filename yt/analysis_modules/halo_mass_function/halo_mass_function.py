"""
halo_mass_function - Halo Mass Function and supporting functions.

Author: Stephen Skory <s@skory.us>
Affiliation: UC San Diego / CASS
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008-2011 Stephen Skory (and others).  All Rights Reserved.

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
import math, time

from yt.funcs import *
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelDummy, \
    ParallelAnalysisInterface, \
    parallel_blocking_call

class HaloMassFcn(ParallelAnalysisInterface):
    def __init__(self, pf, halo_file=None, omega_matter0=None, omega_lambda0=None,
    omega_baryon0=0.05, hubble0=None, sigma8input=0.86, primordial_index=1.0,
    this_redshift=None, log_mass_min=None, log_mass_max=None, num_sigma_bins=360,
    fitting_function=4, mass_column=5):
        """
        Initalize a HaloMassFcn object to analyze the distribution of haloes
        as a function of mass.
        :param halo_file (str): The filename of the output of the Halo Profiler.
        Default=None.
        :param omega_matter0 (float): The fraction of the universe made up of
        matter (dark and baryonic). Default=None.
        :param omega_lambda0 (float): The fraction of the universe made up of
        dark energy. Default=None.
        :param omega_baryon0 (float): The fraction of the universe made up of
        ordinary baryonic matter. This should match the value
        used to create the initial conditions, using 'inits'. This is 
        *not* stored in the enzo datset so it must be checked by hand.
        Default=0.05.
        :param hubble0 (float): The expansion rate of the universe in units of
        100 km/s/Mpc. Default=None.
        :param sigma8input (float): The amplitude of the linear power
        spectrum at z=0 as specified by the rms amplitude of mass-fluctuations
        in a top-hat sphere of radius 8 Mpc/h. This should match the value
        used to create the initial conditions, using 'inits'. This is 
        *not* stored in the enzo datset so it must be checked by hand.
        Default=0.86.
        :param primoridal_index (float): This is the index of the mass power
        spectrum before modification by the transfer function. A value of 1
        corresponds to the scale-free primordial spectrum. This should match
        the value used to make the initial conditions using 'inits'. This is 
        *not* stored in the enzo datset so it must be checked by hand.
        Default=1.0.
        :param this_redshift (float): The current redshift. Default=None.
        :param log_mass_min (float): The log10 of the mass of the minimum of the
        halo mass range. Default=None.
        :param log_mass_max (float): The log10 of the mass of the maximum of the
        halo mass range. Default=None.
        :param num_sigma_bins (float): The number of bins (points) to use for
        the calculations and generated fit. Default=360.
        :param fitting_function (int): Which fitting function to use.
        1 = Press-schechter, 2 = Jenkins, 3 = Sheth-Tormen, 4 = Warren fit
        5 = Tinker
        Default=4.
        :param mass_column (int): The column of halo_file that contains the
        masses of the haloes. Default=4.
        """
        ParallelAnalysisInterface.__init__(self)
        self.pf = pf
        self.halo_file = halo_file
        self.omega_matter0 = omega_matter0
        self.omega_lambda0 = omega_lambda0
        self.omega_baryon0 = omega_baryon0
        self.hubble0 = hubble0
        self.sigma8input = sigma8input
        self.primordial_index = primordial_index
        self.this_redshift = this_redshift
        self.log_mass_min = log_mass_min
        self.log_mass_max = log_mass_max
        self.num_sigma_bins = num_sigma_bins
        self.fitting_function = fitting_function
        self.mass_column = mass_column
        
        # Determine the run mode.
        if halo_file is None:
            # We are hand-picking our various cosmological parameters
            self.mode = 'single'
        else:
            # Make the fit using the same cosmological parameters as the dataset.
            self.mode = 'haloes'
            self.omega_matter0 = self.pf.omega_matter
            self.omega_lambda0 = self.pf.omega_lambda
            self.hubble0 = self.pf.hubble_constant
            self.this_redshift = self.pf.current_redshift
            self.read_haloes()
            if self.log_mass_min == None:
                self.log_mass_min = math.log10(min(self.haloes))
            if self.log_mass_max == None:
                self.log_mass_max = math.log10(max(self.haloes))

        # Input error check.
        if self.mode == 'single':
            if omega_matter0 == None or omega_lambda0 == None or \
            hubble0 == None or this_redshift == None or log_mass_min == None or\
            log_mass_max == None:
                mylog.error("All of these parameters need to be set:")
                mylog.error("[omega_matter0, omega_lambda0, \
                hubble0, this_redshift, log_mass_min, log_mass_max]")
                mylog.error("[%s,%s,%s,%s,%s,%s]" % (omega_matter0,\
                omega_lambda0, hubble0, this_redshift,\
                log_mass_min, log_mass_max))
                return None
        
        # Poke the user to make sure they're doing it right.
        mylog.info(
        """
        Please make sure these are the correct values! They are
        not stored in enzo datasets, so must be entered by hand.
        sigma8input=%f primordial_index=%f omega_baryon0=%f
        """ % (self.sigma8input, self.primordial_index, self.omega_baryon0))
        time.sleep(1)
        
        # Do the calculations.
        self.sigmaM()
        self.dndm()
        
        if self.mode == 'haloes':
            self.bin_haloes()

    def write_out(self, prefix='HMF', fit=True, haloes=True):
        """
        Writes out the halo mass functions to file(s) with prefix *prefix*.
        """
        # First the fit file.
        if fit:
            fitname = prefix + '-fit.dat'
            fp = self.comm.write_on_root(fitname)
            line = \
            """#Columns:
#1. log10 of mass (Msolar, NOT Msolar/h)
#2. mass (Msolar/h)
#3. (dn/dM)*dM (differential number density of haloes, per Mpc^3 (NOT h^3/Mpc^3)
#4. cumulative number density of halos (per Mpc^3, NOT h^3/Mpc^3)
"""
            fp.write(line)
            for i in xrange(self.logmassarray.size - 1):
                line = "%e\t%e\t%e\t%e\n" % (self.logmassarray[i], self.massarray[i],
                self.dn_M_z[i], self.nofmz_cum[i])
                fp.write(line)
            fp.close()
        if self.mode == 'haloes' and haloes:
            haloname = prefix + '-haloes.dat'
            fp = self.comm.write_on_root(haloname)
            line = \
            """#Columns:
#1. log10 of mass (Msolar, NOT Msolar/h)
#2. mass (Msolar/h)
#3. cumulative number density of haloes (per Mpc^3, NOT h^3/Mpc^3)
"""
            fp.write(line)
            for i in xrange(self.logmassarray.size - 1):
                line = "%e\t%e\t%e\n" % (self.logmassarray[i], self.massarray[i],
                self.dis[i])
                fp.write(line)
            fp.close()
        
    def read_haloes(self):
        """
        Read in the virial masses of the haloes.
        """
        mylog.info("Reading halo masses from %s" % self.halo_file)
        f = open(self.halo_file,'r')
        line = f.readline()
        if line == "":
            self.haloes = np.array([])
            return
        while line[0] == '#':
            line = f.readline()
        self.haloes = []
        while line:
            line = line.split()
            mass = float(line[self.mass_column])
            if mass > 0:
                self.haloes.append(float(line[self.mass_column]))
            line = f.readline()
        f.close()
        self.haloes = np.array(self.haloes)

    def bin_haloes(self):
        """
        With the list of virial masses, find the halo mass function.
        """
        bins = np.logspace(self.log_mass_min,
            self.log_mass_max,self.num_sigma_bins)
        avgs = (bins[1:]+bins[:-1])/2.
        dis, bins = np.histogram(self.haloes,bins)
        # add right to left
        for i,b in enumerate(dis):
            dis[self.num_sigma_bins-i-3] += dis[self.num_sigma_bins-i-2]
            if i == (self.num_sigma_bins - 3): break

        self.dis = dis  / self.pf['CosmologyComovingBoxSize']**3.0 * self.hubble0**3.0

    def sigmaM(self):
        """
         Written by BWO, 2006 (updated 25 January 2007).
         Converted to Python by Stephen Skory December 2009.

         This routine takes in cosmological parameters and creates a file (array) with
         sigma(M) in it, which is necessary for various press-schechter type
         stuff.  In principle one can calculate it ahead of time, but it's far,
         far faster in the long run to calculate your sigma(M) ahead of time.
        
         Inputs: cosmology, user must set parameters
        
         Outputs: four columns of data containing the following information:

         1) log mass (Msolar)
         2) mass (Msolar/h)
         3) Radius (comoving Mpc/h)
         4) sigma (normalized) using Msun/h as the input
         
         The arrays output are used later.
        """
        
        # Set up the transfer function object.
        self.TF = TransferFunction(self.omega_matter0, self.omega_baryon0, 0.0, 0,
            self.omega_lambda0, self.hubble0, self.this_redshift);

        if self.TF.qwarn:
            mylog.error("You should probably fix your cosmology parameters!")

        # output arrays
        # 1) log10 of mass (Msolar, NOT Msolar/h)
        self.Rarray = np.empty(self.num_sigma_bins,dtype='float64')
        # 2) mass (Msolar/h)
        self.logmassarray = np.empty(self.num_sigma_bins, dtype='float64')
        # 3) spatial scale corresponding to that radius (Mpc/h)
        self.massarray = np.empty(self.num_sigma_bins, dtype='float64')
        # 4) sigma(M, z=0, where mass is in Msun/h)
        self.sigmaarray = np.empty(self.num_sigma_bins, dtype='float64')

        # get sigma_8 normalization
        R = 8.0;  # in units of Mpc/h (comoving)

        sigma8_unnorm = math.sqrt(self.sigma_squared_of_R(R));
        sigma_normalization = self.sigma8input / sigma8_unnorm;

        rho0 = self.omega_matter0 * 2.78e+11; # in units of h^2 Msolar/Mpc^3

        # spacing in mass of our sigma calculation
        dm = (float(self.log_mass_max) - self.log_mass_min)/self.num_sigma_bins;

        """
         loop over the total number of sigma_bins the user has requested. 
         For each bin, calculate mass and equivalent radius, and call
         sigma_squared_of_R to get the sigma(R) (equivalent to sigma(M)),
         normalize by user-specified sigma_8, and then write out.
        """
        for i in xrange(self.num_sigma_bins):
    
            # thislogmass is in units of Msolar, NOT Msolar/h
            thislogmass = self.log_mass_min +  i*dm
    
            # mass in units of h^-1 Msolar
            thismass = math.pow(10.0, thislogmass) * self.hubble0; 
    
            # radius is in units of h^-1 Mpc (comoving)
            thisradius = math.pow( 3.0*thismass / 4.0 / math.pi / rho0, 1.0/3.0 );
    
            R = thisradius; # h^-1 Mpc (comoving)
    
            self.Rarray[i] = thisradius;  # h^-1 Mpc (comoving)
            self.logmassarray[i] = thislogmass;  # Msun (NOT Msun/h)
            self.massarray[i] = thismass;  # Msun/h
    
            # get normalized sigma(R)
            self.sigmaarray[i] = math.sqrt(self.sigma_squared_of_R(R)) * sigma_normalization;
            # All done!

    def dndm(self):
        
        # constants - set these before calling any functions!
        rho0 = self.omega_matter0 * 2.78e+11; # in units of h^2 Msolar/Mpc^3
        self.delta_c0 = 1.69;  # critical density for turnaround (Press-Schechter)
        
        nofmz_cum = 0.0;  # keep track of cumulative number density
        
        # Loop over masses, going BACKWARD, and calculate dn/dm as well as the 
        # cumulative mass function.
        
        # output arrays
        # 5) (dn/dM)*dM (differential number density of halos, per Mpc^3 (NOT h^3/Mpc^3)
        self.dn_M_z = np.empty(self.num_sigma_bins, dtype='float64')
        # 6) cumulative number density of halos (per Mpc^3, NOT h^3/Mpc^3)
        self.nofmz_cum = np.zeros(self.num_sigma_bins, dtype='float64')
        
        for j in xrange(self.num_sigma_bins - 1):
            i = (self.num_sigma_bins - 2) - j
        
            thissigma = self.sigmaof_M_z(i, self.this_redshift);
            nextsigma = self.sigmaof_M_z(i+1, self.this_redshift);
            
            # calc dsigmadm - has units of h (since massarray has units of h^-1)
            dsigmadm = (nextsigma-thissigma) / (self.massarray[i+1] - self.massarray[i]);

            # calculate dn(M,z) (dn/dM * dM)
            # this has units of h^3 since rho0 has units of h^2, dsigmadm
            # has units of h, and massarray has units of h^-1
            dn_M_z = -1.0 / thissigma * dsigmadm * rho0 / self.massarray[i] * \
            self.multiplicityfunction(thissigma)*(self.massarray[i+1] - self.massarray[i]);

            # scale by h^3 to get rid of all factors of h
            dn_M_z *= math.pow(self.hubble0, 3.0);
            
            # keep track of cumulative number density
            if dn_M_z > 1.0e-20:
                nofmz_cum += dn_M_z;
            
            # Store this.
            self.nofmz_cum[i] = nofmz_cum
            self.dn_M_z[i] = dn_M_z
        

    def sigma_squared_of_R(self, R):
        """
        /* calculates sigma^2(R).  This is the routine where the magic happens (or
           whatever it is that we do here).  Integrates the sigma_squared_integrand
           parameter from R to infinity.  Calls GSL (gnu scientific library) to do
           the actual integration.  
        
           Note that R is in h^-1 Mpc (comoving)
        */
        """
        self.R = R
        result = integrate_inf(self.sigma_squared_integrand)

        sigmasquaredofR = result / 2.0 / math.pi / math.pi

        return sigmasquaredofR;

    def sigma_squared_integrand(self, k):
        """
        /* integrand for integral to get sigma^2(R). */
        """

        Rcom = self.R;  # this is R in comoving Mpc/h

        f = k*k*self.PofK(k)*np.power( abs(self.WofK(Rcom,k)), 2.0);

        return f

    def PofK(self, k):
        """
        /* returns power spectrum as a function of wavenumber k */
        """

        thisPofK = np.power(k, self.primordial_index) * np.power( self.TofK(k), 2.0);

        return thisPofK;

    def TofK(self, k):
        """
        /* returns transfer function as a function of wavenumber k. */
        """
        
        thisTofK = self.TF.TFmdm_onek_hmpc(k);

        return thisTofK;

    def WofK(self, R, k):
        """
        returns W(k), which is the fourier transform of the top-hat function.
        """

        x = R*k;

        thisWofK = 3.0 * ( np.sin(x) - x*np.cos(x) ) / (x*x*x);

        return thisWofK;

    def multiplicityfunction(self, sigma):
        """
        /* Multiplicity function - this is where the various fitting functions/analytic 
        theories are different.  The various places where I found these fitting functions
        are listed below.  */
        """
        
        nu = self.delta_c0 / sigma;
        
        if self.fitting_function==1:
            # Press-Schechter (This form from Jenkins et al. 2001, MNRAS 321, 372-384, eqtn. 5)
            thismult = math.sqrt(2.0/math.pi) * nu * math.exp(-0.5*nu*nu);
        
        elif self.fitting_function==2:
            # Jenkins et al. 2001, MNRAS 321, 372-384, eqtn. 9
            thismult = 0.315 * math.exp( -1.0 * math.pow( abs( math.log(1.0/sigma) + 0.61), 3.8 ) );
        
        elif self.fitting_function==3:
            # Sheth-Tormen 1999, eqtn 10, using expression from Jenkins et al. 2001, eqtn. 7
            A=0.3222;
            a=0.707;
            p=0.3;
            thismult = A*math.sqrt(2.0*a/math.pi)*(1.0+ math.pow( 1.0/(nu*nu*a), p) )*\
            nu * math.exp(-0.5*a*nu*nu);
        
        elif self.fitting_function==4:
            # LANL fitting function - Warren et al. 2005, astro-ph/0506395, eqtn. 5 
            A=0.7234; 
            a=1.625; 
            b=0.2538; 
            c=1.1982;
            thismult = A*( math.pow(sigma, -1.0*a) + b)*math.exp(-1.0*c / sigma / sigma );

        elif self.fitting_function==5:
            # Tinker et al. 2008, eqn 3, \Delta=300 # \Delta=200
            A = 0.2 #0.186
            a = 1.52 #1.47
            b = 2.25 #2.57
            c = 1.27 #1.19
            thismult = A * ( math.pow((sigma / b), -a) + 1) * \
                math.exp(-1 * c / sigma / sigma)
        
        else:
            mylog.error("Don't understand this.  Fitting function requested is %d\n",
            self.fitting_function)
            return None
        
        return thismult

    def sigmaof_M_z(self, sigmabin, redshift):
        """
        /* sigma(M, z) */
        """
        
        thissigma = self.Dofz(redshift) * self.sigmaarray[sigmabin];
        
        return thissigma;

    def Dofz(self, redshift):
        """
        /* Growth function */
        """

        thisDofz = self.gofz(redshift) / self.gofz(0.0) / (1.0+redshift);

        return thisDofz;


    def gofz(self, redshift):
        """
        /* g(z) - I don't think this has any other name*/
        """

        thisgofz = 2.5 * self.omega_matter_of_z(redshift) / \
        ( math.pow( self.omega_matter_of_z(redshift), 4.0/7.0 ) - \
          self.omega_lambda_of_z(redshift) + \
          ( (1.0 + self.omega_matter_of_z(redshift) / 2.0) * \
          (1.0 + self.omega_lambda_of_z(redshift) / 70.0) ))

        return thisgofz;


    def omega_matter_of_z(self,redshift):
        """
        /* Omega matter as a function of redshift */
        """
        
        thisomofz = self.omega_matter0 * math.pow( 1.0+redshift, 3.0) / \
            math.pow( self.Eofz(redshift), 2.0 );
        
        return thisomofz;

    def omega_lambda_of_z(self,redshift):
        """
        /* Omega lambda as a function of redshift */
        """

        thisolofz = self.omega_lambda0 / math.pow( self.Eofz(redshift), 2.0 )

        return thisolofz;

    def Eofz(self, redshift):
        """
        /* E(z) - I don't think this has any other name */
        """
        thiseofz = math.sqrt( self.omega_lambda0 \
            + (1.0 - self.omega_lambda0 - self.omega_matter0)*math.pow( 1.0+redshift, 2.0) \
            + self.omega_matter0 * math.pow( 1.0+redshift, 3.0)  );

        return thiseofz;


""" 
/* Fitting Formulae for CDM + Baryon + Massive Neutrino (MDM) cosmologies. */
/* Daniel J. Eisenstein & Wayne Hu, Institute for Advanced Study */

/* There are two primary routines here, one to set the cosmology, the
other to construct the transfer function for a single wavenumber k. 
You should call the former once (per cosmology) and the latter as 
many times as you want. */

/* TFmdm_set_cosm() -- User passes all the cosmological parameters as
	arguments; the routine sets up all of the scalar quantites needed 
	computation of the fitting formula.  The input parameters are: 
	1) omega_matter -- Density of CDM, baryons, and massive neutrinos,
				in units of the critical density. 
	2) omega_baryon -- Density of baryons, in units of critical. 
	3) omega_hdm    -- Density of massive neutrinos, in units of critical 
	4) degen_hdm    -- (Int) Number of degenerate massive neutrino species 
	5) omega_lambda -- Cosmological constant 
	6) hubble       -- Hubble constant, in units of 100 km/s/Mpc 
	7) redshift     -- The redshift at which to evaluate */

/* TFmdm_onek_mpc() -- User passes a single wavenumber, in units of Mpc^-1.
	Routine returns the transfer function from the Eisenstein & Hu
	fitting formula, based on the cosmology currently held in the 
	internal variables.  The routine returns T_cb (the CDM+Baryon
	density-weighted transfer function), although T_cbn (the CDM+
	Baryon+Neutrino density-weighted transfer function) is stored
	in the global variable tf_cbnu. */

/* We also supply TFmdm_onek_hmpc(), which is identical to the previous
	routine, but takes the wavenumber in units of h Mpc^-1. */

/* We hold the internal scalar quantities in global variables, so that
the user may access them in an external program, via "extern" declarations. */

/* Please note that all internal length scales are in Mpc, not h^-1 Mpc! */
"""

class TransferFunction(object):
    def __init__(self, omega_matter, omega_baryon, omega_hdm,
	    degen_hdm, omega_lambda, hubble, redshift):
        """
        /* This routine takes cosmological parameters and a redshift and sets up
        all the internal scalar quantities needed to compute the transfer function. */
        /* INPUT: omega_matter -- Density of CDM, baryons, and massive neutrinos,
                        in units of the critical density. */
        /* 	  omega_baryon -- Density of baryons, in units of critical. */
        /* 	  omega_hdm    -- Density of massive neutrinos, in units of critical */
        /* 	  degen_hdm    -- (Int) Number of degenerate massive neutrino species */
        /*        omega_lambda -- Cosmological constant */
        /* 	  hubble       -- Hubble constant, in units of 100 km/s/Mpc */
        /*        redshift     -- The redshift at which to evaluate */
        /* OUTPUT: Returns 0 if all is well, 1 if a warning was issued.  Otherwise,
            sets many global variables for use in TFmdm_onek_mpc() */
        """
        self.qwarn = 0;
        self.theta_cmb = 2.728/2.7 # Assuming T_cmb = 2.728 K
    
        # Look for strange input
        if (omega_baryon<0.0):
            mylog.error("TFmdm_set_cosm(): Negative omega_baryon set to trace amount.\n")
            self.qwarn = 1
        if (omega_hdm<0.0):
            mylog.error("TFmdm_set_cosm(): Negative omega_hdm set to trace amount.\n")
            self.qwarn = 1;
        if (hubble<=0.0):
            mylog.error("TFmdm_set_cosm(): Negative Hubble constant illegal.\n")
            return None
        elif (hubble>2.0):
            mylog.error("TFmdm_set_cosm(): Hubble constant should be in units of 100 km/s/Mpc.\n");
            self.qwarn = 1;
        if (redshift<=-1.0):
            mylog.error("TFmdm_set_cosm(): Redshift < -1 is illegal.\n");
            return None
        elif (redshift>99.0):
            mylog.error("TFmdm_set_cosm(): Large redshift entered.  TF may be inaccurate.\n");
            self.qwarn = 1;

        if (degen_hdm<1): degen_hdm=1;
        self.num_degen_hdm = degen_hdm;	
        # Have to save this for TFmdm_onek_mpc()
        # This routine would crash if baryons or neutrinos were zero,
        # so don't allow that.
        if (omega_baryon<=0): omega_baryon=1e-5;
        if (omega_hdm<=0): omega_hdm=1e-5;
    
        self.omega_curv = 1.0-omega_matter-omega_lambda;
        self.omhh = omega_matter*SQR(hubble);
        self.obhh = omega_baryon*SQR(hubble);
        self.onhh = omega_hdm*SQR(hubble);
        self.f_baryon = omega_baryon/omega_matter;
        self.f_hdm = omega_hdm/omega_matter;
        self.f_cdm = 1.0-self.f_baryon-self.f_hdm;
        self.f_cb = self.f_cdm+self.f_baryon;
        self.f_bnu = self.f_baryon+self.f_hdm;
    
        # Compute the equality scale.
        self.z_equality = 25000.0*self.omhh/SQR(SQR(self.theta_cmb)) # Actually 1+z_eq
        self.k_equality = 0.0746*self.omhh/SQR(self.theta_cmb);
    
        # Compute the drag epoch and sound horizon
        z_drag_b1 = 0.313*math.pow(self.omhh,-0.419)*(1+0.607*math.pow(self.omhh,0.674));
        z_drag_b2 = 0.238*math.pow(self.omhh,0.223);
        self.z_drag = 1291*math.pow(self.omhh,0.251)/(1.0+0.659*math.pow(self.omhh,0.828))* \
            (1.0+z_drag_b1*math.pow(self.obhh,z_drag_b2));
        self.y_drag = self.z_equality/(1.0+self.z_drag);
    
        self.sound_horizon_fit = 44.5*math.log(9.83/self.omhh)/math.sqrt(1.0+10.0*math.pow(self.obhh,0.75));
    
        # Set up for the free-streaming & infall growth function 
        self.p_c = 0.25*(5.0-math.sqrt(1+24.0*self.f_cdm));
        self.p_cb = 0.25*(5.0-math.sqrt(1+24.0*self.f_cb));
    
        omega_denom = omega_lambda+SQR(1.0+redshift)*(self.omega_curv+\
                omega_matter*(1.0+redshift));
        self.omega_lambda_z = omega_lambda/omega_denom;
        self.omega_matter_z = omega_matter*SQR(1.0+redshift)*(1.0+redshift)/omega_denom;
        self.growth_k0 = self.z_equality/(1.0+redshift)*2.5*self.omega_matter_z/ \
            (math.pow(self.omega_matter_z,4.0/7.0)-self.omega_lambda_z+ \
            (1.0+self.omega_matter_z/2.0)*(1.0+self.omega_lambda_z/70.0));
        self.growth_to_z0 = self.z_equality*2.5*omega_matter/(math.pow(omega_matter,4.0/7.0) \
            -omega_lambda + (1.0+omega_matter/2.0)*(1.0+omega_lambda/70.0));
        self.growth_to_z0 = self.growth_k0/self.growth_to_z0;	
        
        # Compute small-scale suppression
        self.alpha_nu = self.f_cdm/self.f_cb*(5.0-2.*(self.p_c+self.p_cb))/(5.-4.*self.p_cb)* \
        math.pow(1+self.y_drag,self.p_cb-self.p_c)* \
        (1+self.f_bnu*(-0.553+0.126*self.f_bnu*self.f_bnu))/ \
        (1-0.193*math.sqrt(self.f_hdm*self.num_degen_hdm)+0.169*self.f_hdm*math.pow(self.num_degen_hdm,0.2))* \
        (1+(self.p_c-self.p_cb)/2*(1+1/(3.-4.*self.p_c)/(7.-4.*self.p_cb))/(1+self.y_drag));
        self.alpha_gamma = math.sqrt(self.alpha_nu);
        self.beta_c = 1/(1-0.949*self.f_bnu);
        # Done setting scalar variables
        self.hhubble = hubble # Need to pass Hubble constant to TFmdm_onek_hmpc()
        

    def TFmdm_onek_mpc(self,  kk):
        """
        /* Given a wavenumber in Mpc^-1, return the transfer function for the
        cosmology held in the global variables. */
        /* Input: kk -- Wavenumber in Mpc^-1 */
        /* Output: The following are set as global variables:
            growth_cb -- the transfer function for density-weighted
                    CDM + Baryon perturbations. 
            growth_cbnu -- the transfer function for density-weighted
                    CDM + Baryon + Massive Neutrino perturbations. */
        /* The function returns growth_cb */
        """
    
        self.qq = kk/self.omhh*SQR(self.theta_cmb);
    
        # Compute the scale-dependent growth functions
        self.y_freestream = 17.2*self.f_hdm*(1+0.488*math.pow(self.f_hdm,-7.0/6.0))* \
            SQR(self.num_degen_hdm*self.qq/self.f_hdm);
        temp1 = math.pow(self.growth_k0, 1.0-self.p_cb);
        temp2 = np.power(self.growth_k0/(1+self.y_freestream),0.7);
        self.growth_cb = np.power(1.0+temp2, self.p_cb/0.7)*temp1;
        self.growth_cbnu = np.power(np.power(self.f_cb,0.7/self.p_cb)+temp2, self.p_cb/0.7)*temp1;
    
        # Compute the master function
        self.gamma_eff = self.omhh*(self.alpha_gamma+(1-self.alpha_gamma)/ \
            (1+SQR(SQR(kk*self.sound_horizon_fit*0.43))));
        self.qq_eff = self.qq*self.omhh/self.gamma_eff;
    
        tf_sup_L = np.log(2.71828+1.84*self.beta_c*self.alpha_gamma*self.qq_eff);
        tf_sup_C = 14.4+325/(1+60.5*np.power(self.qq_eff,1.11));
        self.tf_sup = tf_sup_L/(tf_sup_L+tf_sup_C*SQR(self.qq_eff));
    
        self.qq_nu = 3.92*self.qq*math.sqrt(self.num_degen_hdm/self.f_hdm);
        self.max_fs_correction = 1+1.2*math.pow(self.f_hdm,0.64)*math.pow(self.num_degen_hdm,0.3+0.6*self.f_hdm)/ \
            (np.power(self.qq_nu,-1.6)+np.power(self.qq_nu,0.8));
        self.tf_master = self.tf_sup*self.max_fs_correction;
    
        # Now compute the CDM+HDM+baryon transfer functions
        tf_cb = self.tf_master*self.growth_cb/self.growth_k0;
        tf_cbnu = self.tf_master*self.growth_cbnu/self.growth_k0;
        return tf_cb


    def TFmdm_onek_hmpc(self, kk):
        """
        /* Given a wavenumber in h Mpc^-1, return the transfer function for the
        cosmology held in the global variables. */
        /* Input: kk -- Wavenumber in h Mpc^-1 */
        /* Output: The following are set as global variables:
            growth_cb -- the transfer function for density-weighted
                    CDM + Baryon perturbations. 
            growth_cbnu -- the transfer function for density-weighted
                    CDM + Baryon + Massive Neutrino perturbations. */
        /* The function returns growth_cb */
        """
        return self.TFmdm_onek_mpc(kk*self.hhubble);

def SQR(a):
    return a*a

def integrate_inf(fcn, error=1e-3, initial_guess=10):
    """
    Integrate a function *fcn* from zero to infinity, stopping when the answer
    changes by less than *error*. Hopefully someday we can do something
    better than this!
    """
    xvals = np.logspace(0,np.log10(initial_guess), initial_guess+1)-.9
    yvals = fcn(xvals)
    xdiffs = xvals[1:] - xvals[:-1]
    # Trapezoid rule, but with different dxes between values, so np.trapz
    # will not work.
    areas = (yvals[1:] + yvals[:-1]) * xdiffs / 2.0
    area0 = np.sum(areas)
    # Next guess.
    next_guess = 10 * initial_guess
    xvals = np.logspace(0,np.log10(next_guess), 2*initial_guess**2+1)-.99
    yvals = fcn(xvals)
    xdiffs = xvals[1:] - xvals[:-1]
    # Trapezoid rule.
    areas = (yvals[1:] + yvals[:-1]) * xdiffs / 2.0
    area1 = np.sum(areas)
    # Now we refine until the error is smaller than *error*.
    diff = area1 - area0
    area_final = area1
    area_last = area1
    one_pow = 3
    while diff > error:
        next_guess *= 10
        xvals = np.logspace(0,np.log10(next_guess), one_pow*initial_guess**one_pow+1) - (1 - 0.1**one_pow)
        yvals = fcn(xvals)
        xdiffs = xvals[1:] - xvals[:-1]
        # Trapezoid rule.
        areas = (yvals[1:] + yvals[:-1]) * xdiffs / 2.0
        area_next = np.sum(areas)
        diff = area_next - area_last
        area_last = area_next
        one_pow+=1
    return area_last
