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
from yt.units.yt_array import uconcatenate

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
    photons_per_chunk : integer
        The maximum number of photons that are allocated per chunk. Increase or decrease
        as needed.
    """
    def __init__(self, spectral_model, X_H=0.75, Zmet=0.3, photons_per_chunk=10000000):
        self.X_H = X_H
        self.Zmet = Zmet
        self.spectral_model = spectral_model
        self.photons_per_chunk = photons_per_chunk

    def __call__(self, data_source, parameters):
        
        ds = data_source.ds

        exp_time = parameters["FiducialExposureTime"]
        area = parameters["FiducialArea"]
        redshift = parameters["FiducialRedshift"]
        D_A = parameters["FiducialAngularDiameterDistance"].in_cgs()
        dist_fac = 1.0/(4.*np.pi*D_A.value*D_A.value*(1.+redshift)**3)
        src_ctr = parameters["center"]

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

        spectral_norm = area.v*exp_time.v*dist_fac/vol_scale

        for chunk in parallel_objects(citer):

            kT = chunk["kT"].v
            num_cells = len(kT)
            if num_cells == 0:
                continue
            vol = chunk["cell_volume"].in_cgs().v
            EM = (chunk["density"]/mp).v**2
            EM *= 0.5*(1.+self.X_H)*self.X_H*vol

            if isinstance(self.Zmet, basestring):
                metalZ = chunk[self.Zmet].v
            else:
                metalZ = self.Zmet

            idxs = np.argsort(kT)

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

            number_of_photons = np.zeros(num_cells, dtype="uint64")
            energies = np.zeros(self.photons_per_chunk)

            start_e = 0
            end_e = 0

            pbar = get_pbar("Generating photons for chunk ", num_cells)

            for ibegin, iend, ikT in zip(bcell, ecell, kT_idxs):

                kT = kT_bins[ikT] + 0.5*dkT

                n_current = iend-ibegin

                cem = cell_em[ibegin:iend]

                em_sum_c = cem.sum()
                if isinstance(self.Zmet, basestring):
                    em_sum_m = (metalZ[ibegin:iend]*cem).sum()
                else:
                    em_sum_m = metalZ*em_sum_c

                cspec, mspec = self.spectral_model.get_spectrum(kT)

                cumspec_c = np.cumsum(cspec.d)
                tot_ph_c = cumspec_c[-1]*spectral_norm*em_sum_c
                cumspec_c /= cumspec_c[-1]
                cumspec_c = np.insert(cumspec_c, 0, 0.0)

                cumspec_m = np.cumsum(mspec.d)
                tot_ph_m = cumspec_m[-1]*spectral_norm*em_sum_m
                cumspec_m /= cumspec_m[-1]
                cumspec_m = np.insert(cumspec_m, 0, 0.0)

                u = np.random.random(size=n_current)

                cell_norm_c = tot_ph_c*cem/em_sum_c
                cell_n_c = np.uint64(cell_norm_c) + np.uint64(np.modf(cell_norm_c)[0] >= u)
            
                if isinstance(self.Zmet, basestring):
                    cell_norm_m = tot_ph_m*metalZ[ibegin:iend]*cem/em_sum_m
                else:
                    cell_norm_m = tot_ph_m*metalZ*cem/em_sum_m
                cell_n_m = np.uint64(cell_norm_m) + np.uint64(np.modf(cell_norm_m)[0] >= u)
            
                number_of_photons[ibegin:iend] = cell_n_c + cell_n_m

                end_e += int((cell_n_c+cell_n_m).sum())

                if end_e > self.photons_per_chunk:
                    raise RuntimeError("Number of photons generated for this chunk "+
                                       "exceeds photons_per_chunk (%d)! " % self.photons_per_chunk +
                                       "Increase photons_per_chunk!")

                energies[start_e:end_e] = _generate_energies(cell_n_c, cell_n_m, cumspec_c, cumspec_m, energy)
            
                start_e = end_e

                pbar.update(iend)

            pbar.finish()

            active_cells = number_of_photons > 0
            idxs = idxs[active_cells]

            mylog.info("Number of photons generated for this chunk: %d" % int(number_of_photons.sum()))
            mylog.info("Number of cells with photons: %d" % int(active_cells.sum()))

            photons["NumberOfPhotons"].append(number_of_photons[active_cells])
            photons["Energy"].append(ds.arr(energies[:end_e].copy(), "keV"))
            photons["x"].append((chunk["x"][idxs]-src_ctr[0]).in_units("kpc"))
            photons["y"].append((chunk["y"][idxs]-src_ctr[1]).in_units("kpc"))
            photons["z"].append((chunk["z"][idxs]-src_ctr[2]).in_units("kpc"))
            photons["vx"].append(chunk["velocity_x"][idxs].in_units("km/s"))
            photons["vy"].append(chunk["velocity_y"][idxs].in_units("km/s"))
            photons["vz"].append(chunk["velocity_z"][idxs].in_units("km/s"))
            photons["dx"].append(chunk["dx"][idxs].in_units("kpc"))

        for key in photons:
            if len(photons[key]) > 0:
                photons[key] = uconcatenate(photons[key])

        return photons

def _generate_energies(cell_n_c, cell_n_m, counts_c, counts_m, energy):
    energies = np.array([])
    for cn_c, cn_m in zip(cell_n_c, cell_n_m):
        if cn_c > 0:
            randvec_c = np.random.uniform(size=cn_c)
            randvec_c.sort()
            cell_e_c = np.interp(randvec_c, counts_c, energy)
            energies = np.append(energies, cell_e_c)
        if cn_m > 0: 
            randvec_m = np.random.uniform(size=cn_m)
            randvec_m.sort()
            cell_e_m = np.interp(randvec_m, counts_m, energy)
            energies = np.append(energies, cell_e_m)
    return energies
