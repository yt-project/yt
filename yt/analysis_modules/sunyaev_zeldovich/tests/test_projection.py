"""
Unit test the sunyaev_zeldovich analysis module.
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.frontends.stream.api import load_uniform_grid
from yt.funcs import get_pbar
from yt.utilities.physical_constants import cm_per_kpc, K_per_keV, \
     mh, cm_per_km, kboltz, Tcmb, hcgs, clight, sigma_thompson
from yt.testing import *
from yt.analysis_modules.sunyaev_zeldovich.projection import SZProjection, I0
import numpy as np
try:
    import SZpack
except ImportError:
    pass

mue = 1./0.88
freqs = np.array([30., 90., 240.])
    
def setup():
    pass

def full_szpack3d(pf, xo):
    data = pf.h.grids[0]
    dz = pf.h.get_smallest_dx()*pf.units["cm"]
    nx,ny,nz = data["Density"].shape
    dn = np.zeros((nx,ny,nz))
    Dtau = sigma_thompson*data["Density"]/(mh*mue)*dz
    Te = data["Temperature"]/K_per_keV
    betac = data["z-velocity"]/clight
    pbar = get_pbar("Computing 3-D cell-by-cell S-Z signal for comparison.", nx) 
    for i in xrange(nx):
        pbar.update(i)
        for j in xrange(ny):
            for k in xrange(nz):
                dn[i,j,k] = SZpack.compute_3d(xo, Dtau[i,j,k],
                                              Te[i,j,k], betac[i,j,k],
                                              1.0, 0.0, 0.0, 1.0e-5)
    pbar.finish()
    return I0*xo**3*np.sum(dn, axis=2)

def setup_cluster():

    R = 1000.
    r_c = 100.
    rho_c = 1.673e-26
    beta = 1.
    T0 = 4.
    nx,ny,nz = 16,16,16
    c = 0.17
    a_c = 30.
    a = 200.
    v0 = 300.*cm_per_km
    ddims = (nx,ny,nz)
    
    x, y, z = np.mgrid[-R:R:nx*1j,
                       -R:R:ny*1j,
                       -R:R:nz*1j]

    r = np.sqrt(x**2+y**2+z**2)

    dens = np.zeros(ddims)
    dens = rho_c*(1.+(r/r_c)**2)**(-1.5*beta)
    temp = T0*K_per_keV/(1.+r/a)*(c+r/a_c)/(1.+r/a_c)
    velz = v0*temp/(T0*K_per_keV)

    data = {}
    data["Density"] = dens
    data["Temperature"] = temp
    data["x-velocity"] = np.zeros(ddims)
    data["y-velocity"] = np.zeros(ddims)
    data["z-velocity"] = velz

    bbox = np.array([[-0.5,0.5],[-0.5,0.5],[-0.5,0.5]])
    
    L = 2*R*cm_per_kpc
    dl = L/nz

    pf = load_uniform_grid(data, ddims, L, bbox=bbox)

    return pf

pf = setup_cluster()

@requires_module("SZpack")
def test_projection():
    nx,ny,nz = pf.domain_dimensions
    xinit = 1.0e9*hcgs*freqs/(kboltz*Tcmb)
    szprj = SZProjection(pf, freqs, mue=mue, high_order=True)
    szprj.on_axis(2, nx=nx)
    deltaI = np.zeros((3,nx,ny))
    for i in xrange(3):
        deltaI[i,:,:] = full_szpack3d(pf, xinit[i])
        yield assert_almost_equal, deltaI[i,:,:], szprj["%d_GHz" % int(freqs[i])], 6
        
