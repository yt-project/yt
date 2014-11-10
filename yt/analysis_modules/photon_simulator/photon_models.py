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
from yt.utilities.physical_constants import mp, kboltz
from yt.utilities.parallel_tools.parallel_analysis_interface import \
     communication_system, parallel_objects

n_kT = 10000
kT_min = 8.08e-2
kT_max = 50.

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
        
        ds = data_source.ds

        exp_time = parameters["FiducialExposureTime"]
        area = parameters["FiducialArea"]
        redshift = parameters["FiducialRedshift"]
        D_A = parameters["FiducialAngularDiameterDistance"].in_cgs()
        dist_fac = 1.0/(4.*np.pi*D_A.value*D_A.value*(1.+redshift)**3)
                
        vol_scale = 1.0/np.prod(ds.domain_width.in_cgs().to_ndarray())

        my_kT_min, my_kT_max = data_source.quantities.extrema("kT")

        self.spectral_model.prepare()
        energy = self.spectral_model.ebins

        citer = data_source.chunks(["kT","cell_volume","density",
                                    "x","y","z","dx","velocity_x",
                                    "velocity_y","velocity_z"], "io")

        photons = {}
        photons["x"] = []
        photons["y"] = []
        photons["z"] = []
        photons["vx"] = []
        photons["vy"] = []
        photons["vz"] = []
        photons["dx"] = []
        photons["Energy"] = []
        photons["NumberOfPhotons"] = []

        for chunk in parallel_objects(citer):

            kT = chunk["kT"].v
            if len(kT) == 0:
                continue
            vol = chunk["cell_volume"].in_cgs().v
            EM = (chunk["density"]/mp).v**2
            EM *= 0.5*(1.+self.X_H)*self.X_H*vol

            if isinstance(self.Zmet, basestring):
                metalZ = chunk[self.Zmet].v
            else:
                metalZ = self.Zmet*chunk["ones"]

            idxs = np.argsort(kT)
            dshape = idxs.shape

            kT_bins = np.linspace(kT_min, max(my_kT_max, kT_max), num=n_kT+1)
            dkT = kT_bins[1]-kT_bins[0]
            kT_idxs = np.digitize(kT[idxs], kT_bins)
            kT_idxs = np.minimum(np.maximum(1, kT_idxs), n_kT) - 1
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

            cell_em = EM[idxs]*vol_scale

            u = np.random.random(cell_em.shape)

            pbar = get_pbar("Generating Photons", n_kT)

            for i, ikT in enumerate(kT_idxs):

                ibegin = bcell[i]
                iend = ecell[i]
                kT = kT_bins[ikT] + 0.5*dkT
        
                em_sum_c = cell_em[ibegin:iend].sum()
                em_sum_m = (metalZ[ibegin:iend]*cell_em[ibegin:iend]).sum()

                cspec, mspec = self.spectral_model.get_spectrum(kT)
                cspec *= dist_fac*em_sum_c/vol_scale
                mspec *= dist_fac*em_sum_m/vol_scale

                cumspec_c = np.cumsum(cspec.ndarray_view())
                counts_c = cumspec_c[:]/cumspec_c[-1]
                counts_c = np.insert(counts_c, 0, 0.0)
                tot_ph_c = cumspec_c[-1]*area.value*exp_time.value

                cumspec_m = np.cumsum(mspec.ndarray_view())
                counts_m = cumspec_m[:]/cumspec_m[-1]
                counts_m = np.insert(counts_m, 0, 0.0)
                tot_ph_m = cumspec_m[-1]*area.value*exp_time.value

                for icell in xrange(ibegin, iend):
            
                    cell_norm_c = tot_ph_c*cell_em[icell]/em_sum_c
                    cell_n_c = np.uint64(cell_norm_c) + np.uint64(np.modf(cell_norm_c)[0] >= u[icell])
            
                    cell_norm_m = tot_ph_m*metalZ[icell]*cell_em[icell]/em_sum_m
                    cell_n_m = np.uint64(cell_norm_m) + np.uint64(np.modf(cell_norm_m)[0] >= u[icell])
            
                    cell_n = cell_n_c + cell_n_m

                    if cell_n > 0:
                        randvec_c = np.random.uniform(size=cell_n_c)
                        randvec_c.sort()
                        randvec_m = np.random.uniform(size=cell_n_m)
                        randvec_m.sort()
                        cell_e_c = np.interp(randvec_c, counts_c, energy)
                        cell_e_m = np.interp(randvec_m, counts_m, energy)
                        photons["x"].append(chunk["x"][icell])
                        photons["y"].append(chunk["y"][icell])
                        photons["z"].append(chunk["z"][icell])
                        photons["dx"].append(chunk["dx"][icell])
                        photons["vx"].append(chunk["velocity_x"][icell])
                        photons["vy"].append(chunk["velocity_y"][icell])
                        photons["vz"].append(chunk["velocity_z"][icell])
                        photons["NumberOfPhotons"].append(cell_n)
                        photons["Energy"].append(np.concatenate([cell_e_c,cell_e_m]))

            
                pbar.update(i)

            pbar.finish()

        src_ctr = parameters["center"]

        photons["x"] = (ds.arr(photons["x"])-src_ctr[0]).in_units("kpc")
        photons["y"] = (ds.arr(photons["y"])-src_ctr[1]).in_units("kpc")
        photons["z"] = (ds.arr(photons["z"])-src_ctr[2]).in_units("kpc")
        photons["vx"] = ds.arr(photons["vx"]).in_units("km/s")
        photons["vy"] = ds.arr(photons["vy"]).in_units("km/s")
        photons["vz"] = ds.arr(photons["vz"]).in_units("km/s")
        photons["dx"] = ds.arr(photons["dx"]).in_units("kpc")
        photons["Energy"] = ds.arr(np.concatenate(photons["Energy"]), "keV")
    
        return photons
