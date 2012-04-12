"""
Definitions for ART files


Author: Christopher Erick Moody <juxtaposicion@gmail.com>
Affiliation: UC Santa Cruz
Homepage: http://yt-project.org/
License:
  Copyright (C) 2012 Matthew Turk.  All Rights Reserved.

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

# These functions are ART-specific

import numpy as na

field_list = [ 'Density','TotalEnergy',
             'XMomentumDensity','YMomentumDensity','ZMomentumDensity',
             'Pressure','Gamma','GasEnergy',
             'MetalDensitySNII', 'MetalDensitySNIa',
             'PotentialNew','PotentialOld']

amr_header_struct = [
    ('>i','pad byte'),
    ('>256s','jname'),
    ('>i','pad byte'),
    
    ('>i','pad byte'),
    ('>i','istep'),
    ('>d','t'),
    ('>d','dt'),
    ('>f','aexpn'),
    ('>f','ainit'),
    ('>i','pad byte'),
    
    ('>i','pad byte'),
    ('>f','boxh'),
    ('>f','Om0'),
    ('>f','Oml0'),
    ('>f','Omb0'),
    ('>f','hubble'),
    ('>i','pad byte'),
    
    ('>i','pad byte'),
    ('>i','nextras'),
    ('>i','pad byte'),

    ('>i','pad byte'),
    ('>f','extra1'),
    ('>f','extra2'),
    ('>i','pad byte'),

    ('>i','pad byte'),
    ('>256s','lextra'),
    ('>256s','lextra'),
    ('>i','pad byte'),
    
    ('>i', 'pad byte'),
    ('>i', 'min_level'),
    ('>i', 'max_level'),
    ('>i', 'pad byte'),
    ]


particle_header_struct = [
    ('>i','pad'),
    ('45s','header'), 
    ('>f','aexpn'),
    ('>f','aexp0'),
    ('>f','amplt'),
    ('>f','astep'),

    ('>i','istep'),
    ('>f','partw'),
    ('>f','tintg'),

    ('>f','Ekin'),
    ('>f','Ekin1'),
    ('>f','Ekin2'),
    ('>f','au0'),
    ('>f','aeu0'),


    ('>i','Nrow'),
    ('>i','Ngridc'),
    ('>i','Nspecies'),
    ('>i','Nseed'),

    ('>f','Om0'),
    ('>f','Oml0'),
    ('>f','hubble'),
    ('>f','Wp5'),
    ('>f','Ocurv'),
    ('>f','Omb0'),
    ('>%ds'%(396),'extras'),
    ('>f','unknown'),

    ('>i','pad')]

#Taken from Kravtsov's thesis
paramteres = {}
parameters["Y_p"] = 0.245
parameters["wmu"] = 4.0/(8.0-5.0*parameters["Y_p"])
parameters["gamma"] = 5./3.
parameters["T_CMB0"] = 2.726  
parameters["T_min"] = 300.0 #T floor in K
parameters['ng'] = 128 # of 0 level cells in 1d 

def update_parameters(params):
    """Using code constants and supplied header values evaluate new dependent parameters"""
    
    #these are all defined in ART_analysis.F or 
    #the Kravtsov thesis
    Om0 = params['Om0']
    hubble = params['hubble']
    dummy = 100.0 * hubble * na.sqrt(Om0)
    ng = params['ng']
    wmu = params["wmu"]
    boxh = params['boxh'] 
    
    #distance unit #boxh is units of h^-1 Mpc
    params["r0"] = params["boxh"] / params['ng']
    r0 = params["r0"]
    #time, yrs
    params["t0"] = 2.0 / dummy * 3.0856e19 / 3.15e7
    #velocity velocity units in km/s
    params["v0"] = 50.0*params["r0"]*\
            na.sqrt(params["Om0"])
    #density = 3H0^2 * Om0 / (8*pi*G) - unit of density in Msun/Mpc^3
    params["rho0"] = 2.776e11 * hubble**2.0 * Om0
    rho0 = params["rho0"]
    #Pressure = rho0 * v0**2 - unit of pressure in g/cm/s^2
    params["P0"] = 4.697e-16 * Om0**2.0 * r0**2.0 * hubble**2.0
    #T_0 = unit of temperature in K and in keV)
    #T_0 = 2.61155 * r0**2 * wmu * Om0 ! [keV]
    params["T_0"] = 3.03e5 * r0**2.0 * wmu * Om0 # [K]
    #S_0 = unit of entropy in keV * cm^2
    params["S_0"] = 52.077 * wmu**(5.0/3.0) * hubble**(-4.0/3.0)*Om0**(1.0/3.0)*r0**2.0
    
    #mass conversion (Mbox = rho0 * Lbox^3, Mbox_code = Ng^3
    #     for non-cosmological run aM0 must be defined during initialization
    #     [aM0] = [Msun]
    params["aM0"] = rho0 * (boxh/hubble)**3.0 / ng**3.0
    return params

