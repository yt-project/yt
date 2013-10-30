"""
Unit test the photon_simulator analysis module.
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.frontends.stream.api import load_uniform_grid
from yt.testing import *
from yt.utilities.physical_constants import cm_per_kpc, \
     K_per_keV, cm_per_mpc, mp
from yt.utilities.cosmology import Cosmology
from yt.analysis_modules.api import PhotonList, EventList, \
     XSpecThermalModel, XSpecAbsorbModel, ThermalPhotonModel
import os
import xspec

def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"
                
@requires_module("xspec")
def test_beta_model():

    # Set up the beta model and stream dataset
    R = 1000.
    r_c = 100.
    rho_c = 1.673e-26
    beta = 1.
    T = 4.
    nx = 256
    nH = 0.1
    Zmet = 0.3
    X_H = 0.75
    nenp0 = (rho_c/mp)**2*0.5*(1.+X_H)*X_H

    ddims = (nx,nx,nx)
    
    x, y, z = np.mgrid[-R:R:nx*1j,
                       -R:R:nx*1j,
                       -R:R:nx*1j]
    
    r = np.sqrt(x**2+y**2+z**2)

    dens = np.zeros(ddims)
    dens[r <= R] = rho_c*(1.+(r[r <= R]/r_c)**2)**(-1.5*beta)
    dens[r > R] = 0.0
    temp = T*K_per_keV*np.ones(ddims)

    data = {}
    data["Density"] = dens
    data["Temperature"] = temp
    data["x-velocity"] = np.zeros(ddims)
    data["y-velocity"] = np.zeros(ddims)
    data["z-velocity"] = np.zeros(ddims)

    bbox = np.array([[-0.5,0.5],[-0.5,0.5],[-0.5,0.5]])
    
    pf = load_uniform_grid(data, ddims, 2*R*cm_per_kpc, bbox=bbox)

    # Grab a sphere data object
    
    sphere = pf.h.sphere(pf.domain_center, 1.0/pf["mpc"])

    # Create the photons

    ARF = os.environ["YT_DATA_DIR"]+"xray_data/chandra_ACIS-S3_onaxis_arf.fits"
    RMF = os.environ["YT_DATA_DIR"]+"xray_data/chandra_ACIS-S3_onaxis_rmf.fits"
                    
    A = 6000.
    exp_time = 1.0e5
    redshift = 0.05
    cosmo = Cosmology()
    DA = cosmo.AngularDiameterDistance(0.0,redshift)*cm_per_mpc
    EM = 4.*np.pi*nenp0*0.196022123*(r_c*cm_per_kpc)**3
    norm = 1.0e-14*EM/(4.*np.pi*DA**2*(1.+redshift)**2)

    apec_model = XSpecThermalModel("apec", 0.01, 20.0, 10000)
    abs_model  = XSpecAbsorbModel("TBabs", nH)
        
    thermal_model = ThermalPhotonModel(apec_model, Zmet=Zmet)
    photons = PhotonList.from_scratch(sphere, redshift, A, exp_time,
                                      thermal_model, cosmology=cosmo)
    
    events = photons.project_photons([0.0,0.0,1.0],
                                     responses=[ARF,RMF],
                                     absorb_model=abs_model)
    events.write_spectrum("spec_chandra.fits", clobber=True)

    # Now fit the resulting spectrum
    
    spec = xspec.Spectrum("spec_chandra.fits")
    xspec.Fit.statMethod = "cstat"
    
    spec.ignore("**-0.5")
    spec.ignore("7.0-**")
    
    m = xspec.Model("tbabs*apec")
    m.TBabs.nH = 0.09
    m.apec.kT = 5.
    m.apec.Abundanc = 0.2
    m.apec.Abundanc.frozen = False
    m.apec.Redshift = redshift
    
    xspec.Fit.renorm()
    xspec.Fit.perform()
    xspec.Fit.error("1-3,5")
    
    assert(T > m.apec.kT.error[0] and T < m.apec.kT.error[1])
    assert(Zmet > m.apec.Abundanc.error[0] and Zmet < m.apec.Abundanc.error[1])
    assert(norm > m.apec.norm.error[0] and norm < m.apec.norm.error[1])
