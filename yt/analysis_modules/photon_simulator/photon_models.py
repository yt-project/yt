"""
Classes for specific photon models

The algorithms used here are based off of the method used by the
PHOX code (http://www.mpa-garching.mpg.de/~kdolag/Phox/),
developed by Veronica Biffi and Klaus Dolag. References for
PHOX may be found at:

Biffi, V., Dolag, K., Bohringer, H., & Lemson, G. 2012, MNRAS, 420, 3545
http://adsabs.harvard.edu/abs/2012MNRAS.420.3545B

Biffi, V., Dolag, K., Bohringer, H. 2013, MNRAS, 428, 1395
http://adsabs.harvard.edu/abs/2013MNRAS.428.1395B

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from yt.funcs import *
from yt.utilities.physical_constants import \
     mp, cm_per_km, K_per_keV, cm_per_mpc
from yt.utilities.parallel_tools.parallel_analysis_interface import \
     communication_system

N_TBIN = 10000
TMIN = 8.08e-2
TMAX = 50.

comm = communication_system.communicators[-1]

class PhotonModel(object):

    def __init__(self):
        pass

    def __call__(self, data_source, parameters):
        photons = {}
        return photons

class ThermalPhotonModel(PhotonModel):
    r"""
    Initialize a ThermalPhotonModel from a thermal spectrum. 
    
    Parameters
    ----------

    spectral_model : `SpectralModel`
        A thermal spectral model instance, either of `XSpecThermalModel`
        or `TableApecModel`. 
    X_H : float, optional
        The hydrogen mass fraction.
    Zmet : float or string, optional
        The metallicity. If a float, assumes a constant metallicity throughout.
        If a string, is taken to be the name of the metallicity field.
    """
    def __init__(self, spectral_model, X_H=0.75, Zmet=0.3):
        self.X_H = X_H
        self.Zmet = Zmet
        self.spectral_model = spectral_model

    def __call__(self, data_source, parameters):
        
        pf = data_source.pf

        exp_time = parameters["FiducialExposureTime"]
        area = parameters["FiducialArea"]
        redshift = parameters["FiducialRedshift"]
        D_A = parameters["FiducialAngularDiameterDistance"]*cm_per_mpc
        dist_fac = 1.0/(4.*np.pi*D_A*D_A*(1.+redshift)**3)
                
        vol_scale = pf.units["cm"]**(-3)/np.prod(pf.domain_width)
        
        num_cells = data_source["Temperature"].shape[0]
        start_c = comm.rank*num_cells/comm.size
        end_c = (comm.rank+1)*num_cells/comm.size
        
        kT = data_source["Temperature"][start_c:end_c].copy()/K_per_keV
        vol = data_source["CellVolume"][start_c:end_c].copy()
        dx = data_source["dx"][start_c:end_c].copy()
        EM = (data_source["Density"][start_c:end_c].copy()/mp)**2
        EM *= 0.5*(1.+self.X_H)*self.X_H*vol
    
        data_source.clear_data()
    
        x = data_source["x"][start_c:end_c].copy()
        y = data_source["y"][start_c:end_c].copy()
        z = data_source["z"][start_c:end_c].copy()
    
        data_source.clear_data()
        
        vx = data_source["x-velocity"][start_c:end_c].copy()
        vy = data_source["y-velocity"][start_c:end_c].copy()
        vz = data_source["z-velocity"][start_c:end_c].copy()
    
        if isinstance(self.Zmet, basestring):
            metalZ = data_source[self.Zmet][start_c:end_c].copy()
        else:
            metalZ = self.Zmet*np.ones(EM.shape)
        
        data_source.clear_data()

        idxs = np.argsort(kT)
        dshape = idxs.shape

        kT_bins = np.linspace(TMIN, max(kT[idxs][-1], TMAX), num=N_TBIN+1)
        dkT = kT_bins[1]-kT_bins[0]
        kT_idxs = np.digitize(kT[idxs], kT_bins)
        kT_idxs = np.minimum(np.maximum(1, kT_idxs), N_TBIN) - 1
        bcounts = np.bincount(kT_idxs).astype("int")
        bcounts = bcounts[bcounts > 0]
        n = int(0)
        bcell = []
        ecell = []
        for bcount in bcounts:
            bcell.append(n)
            ecell.append(n+bcount)
            n += bcount
        kT_idxs = np.unique(kT_idxs)
        
        self.spectral_model.prepare()
        energy = self.spectral_model.ebins
    
        cell_em = EM[idxs]*vol_scale
        cell_vol = vol[idxs]*vol_scale
    
        number_of_photons = np.zeros(dshape, dtype='uint64')
        energies = []
    
        u = np.random.random(cell_em.shape)
        
        pbar = get_pbar("Generating Photons", dshape[0])

        for i, ikT in enumerate(kT_idxs):

            ncells = int(bcounts[i])
            ibegin = bcell[i]
            iend = ecell[i]
            kT = kT_bins[ikT] + 0.5*dkT
        
            em_sum_c = cell_em[ibegin:iend].sum()
            em_sum_m = (metalZ[ibegin:iend]*cell_em[ibegin:iend]).sum()
            
            cspec, mspec = self.spectral_model.get_spectrum(kT)
            cspec *= dist_fac*em_sum_c/vol_scale
            mspec *= dist_fac*em_sum_m/vol_scale
        
            cumspec_c = np.cumsum(cspec)
            counts_c = cumspec_c[:]/cumspec_c[-1]
            counts_c = np.insert(counts_c, 0, 0.0)
            tot_ph_c = cumspec_c[-1]*area*exp_time

            cumspec_m = np.cumsum(mspec)
            counts_m = cumspec_m[:]/cumspec_m[-1]
            counts_m = np.insert(counts_m, 0, 0.0)
            tot_ph_m = cumspec_m[-1]*area*exp_time
        
            for icell in xrange(ibegin, iend):
            
                cell_norm_c = tot_ph_c*cell_em[icell]/em_sum_c
                cell_n_c = np.uint64(cell_norm_c) + np.uint64(np.modf(cell_norm_c)[0] >= u[icell])
            
                cell_norm_m = tot_ph_m*metalZ[icell]*cell_em[icell]/em_sum_m
                cell_n_m = np.uint64(cell_norm_m) + np.uint64(np.modf(cell_norm_m)[0] >= u[icell])
            
                cell_n = cell_n_c + cell_n_m

                if cell_n > 0:
                    number_of_photons[icell] = cell_n
                    randvec_c = np.random.uniform(size=cell_n_c)
                    randvec_c.sort()
                    randvec_m = np.random.uniform(size=cell_n_m)
                    randvec_m.sort()
                    cell_e_c = np.interp(randvec_c, counts_c, energy)
                    cell_e_m = np.interp(randvec_m, counts_m, energy)
                    energies.append(np.concatenate([cell_e_c,cell_e_m]))
                
                pbar.update(icell)

        pbar.finish()
            
        active_cells = number_of_photons > 0
        idxs = idxs[active_cells]
        
        photons = {}

        src_ctr = parameters["center"]
        
        photons["x"] = (x[idxs]-src_ctr[0])*pf.units["kpc"]
        photons["y"] = (y[idxs]-src_ctr[1])*pf.units["kpc"]
        photons["z"] = (z[idxs]-src_ctr[2])*pf.units["kpc"]
        photons["vx"] = vx[idxs]/cm_per_km
        photons["vy"] = vy[idxs]/cm_per_km
        photons["vz"] = vz[idxs]/cm_per_km
        photons["dx"] = dx[idxs]*pf.units["kpc"]
        photons["NumberOfPhotons"] = number_of_photons[active_cells]
        photons["Energy"] = np.concatenate(energies)
    
        return photons
