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

from yt.extern.six import string_types
import numpy as np
from yt.funcs import mylog, get_pbar
from yt.units.yt_array import YTArray
from yt.utilities.physical_constants import mp
from yt.utilities.parallel_tools.parallel_analysis_interface import \
     parallel_objects
from yt.units.yt_array import uconcatenate

n_kT = 10000
kT_min = 8.08e-2
kT_max = 50.

photon_units = {"Energy":"keV",
                "dx":"kpc"}
for ax in "xyz":
    photon_units[ax] = "kpc"
    photon_units["v"+ax] = "km/s"

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
    method : string, optional
        The method used to generate the photon energies from the spectrum:
        "invert_cdf": Invert the cumulative distribution function of the spectrum.
        "accept_reject": Acceptance-rejection method using the spectrum. 
        The first method should be sufficient for most cases. 
    prng : NumPy `RandomState` object or numpy.random                                                                            
        A pseudo-random number generator. Typically will only be specified                                                            
        if you have a reason to generate the same set of random numbers, such as for a 
        test. Default is the numpy.random module.                                                                      
    """
    def __init__(self, spectral_model, X_H=0.75, Zmet=0.3, 
                 photons_per_chunk=10000000, method="invert_cdf",
                 prng=np.random):
        self.X_H = X_H
        self.Zmet = Zmet
        self.spectral_model = spectral_model
        self.photons_per_chunk = photons_per_chunk
        self.method = method
        self.prng = prng

    def __call__(self, data_source, parameters):

        ds = data_source.ds

        exp_time = parameters["FiducialExposureTime"]
        area = parameters["FiducialArea"]
        redshift = parameters["FiducialRedshift"]
        D_A = parameters["FiducialAngularDiameterDistance"].in_cgs()
        dist_fac = 1.0/(4.*np.pi*D_A.value*D_A.value*(1.+redshift)**2)
        src_ctr = parameters["center"]

        my_kT_min, my_kT_max = data_source.quantities.extrema("kT")

        self.spectral_model.prepare_spectrum(redshift)
        emid = self.spectral_model.emid
        ebins = self.spectral_model.ebins
        nchan = len(emid)

        citer = data_source.chunks([], "io")

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

        spectral_norm = area.v*exp_time.v*dist_fac

        tot_num_cells = data_source.ires.shape[0]

        pbar = get_pbar("Generating photons ", tot_num_cells)

        cell_counter = 0

        for chunk in parallel_objects(citer):

            kT = chunk["kT"].v
            num_cells = len(kT)
            if num_cells == 0:
                continue
            vol = chunk["cell_volume"].in_cgs().v
            EM = (chunk["density"]/mp).in_cgs().v**2
            EM *= 0.5*(1.+self.X_H)*self.X_H*vol

            if isinstance(self.Zmet, string_types):
                metalZ = chunk[self.Zmet].v
            else:
                metalZ = self.Zmet*np.ones(num_cells)

            idxs = np.argsort(kT)

            kT_bins = np.linspace(kT_min, max(my_kT_max.v, kT_max), num=n_kT+1)
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

            cell_em = EM[idxs]*spectral_norm

            number_of_photons = np.zeros(num_cells, dtype="uint64")
            energies = np.zeros(self.photons_per_chunk)

            start_e = 0
            end_e = 0

            for ibegin, iend, ikT in zip(bcell, ecell, kT_idxs):

                kT = kT_bins[ikT] + 0.5*dkT

                n_current = iend-ibegin

                cem = cell_em[ibegin:iend]

                cspec, mspec = self.spectral_model.get_spectrum(kT)

                tot_ph_c = cspec.d.sum()
                tot_ph_m = mspec.d.sum()

                u = self.prng.uniform(size=n_current)

                cell_norm_c = tot_ph_c*cem
                cell_norm_m = tot_ph_m*metalZ[ibegin:iend]*cem
                cell_norm = np.modf(cell_norm_c + cell_norm_m)
                cell_n = np.uint64(cell_norm[1]) + np.uint64(cell_norm[0] >= u)

                number_of_photons[ibegin:iend] = cell_n

                end_e += int(cell_n.sum())

                if end_e > self.photons_per_chunk:
                    raise RuntimeError("Number of photons generated for this chunk "+
                                       "exceeds photons_per_chunk (%d)! " % self.photons_per_chunk +
                                       "Increase photons_per_chunk!")

                if self.method == "invert_cdf":
                    cumspec_c = np.cumsum(cspec.d)
                    cumspec_m = np.cumsum(mspec.d)
                    cumspec_c = np.insert(cumspec_c, 0, 0.0)
                    cumspec_m = np.insert(cumspec_m, 0, 0.0)

                ei = start_e
                for cn, Z in zip(number_of_photons[ibegin:iend], metalZ[ibegin:iend]):
                    if cn == 0: continue
                    # The rather verbose form of the few next statements is a
                    # result of code optimization and shouldn't be changed
                    # without checking for perfomance degradation. See
                    # https://bitbucket.org/yt_analysis/yt/pull-requests/1766
                    # for details.
                    if self.method == "invert_cdf":
                        cumspec = cumspec_c
                        cumspec += Z * cumspec_m
                        norm_factor = 1.0 / cumspec[-1]
                        cumspec *= norm_factor
                        randvec = self.prng.uniform(size=cn)
                        randvec.sort()
                        cell_e = np.interp(randvec, cumspec, ebins)
                    elif self.method == "accept_reject":
                        tot_spec = cspec.d
                        tot_spec += Z * mspec.d
                        norm_factor = 1.0 / tot_spec.sum()
                        tot_spec *= norm_factor
                        eidxs = self.prng.choice(nchan, size=cn, p=tot_spec)
                        cell_e = emid[eidxs]
                    energies[ei:ei+cn] = cell_e
                    cell_counter += 1
                    pbar.update(cell_counter)
                    ei += cn

                start_e = end_e

            active_cells = number_of_photons > 0
            idxs = idxs[active_cells]

            photons["NumberOfPhotons"].append(number_of_photons[active_cells])
            photons["Energy"].append(ds.arr(energies[:end_e].copy(), "keV"))
            photons["x"].append((chunk["x"][idxs]-src_ctr[0]).in_units("kpc"))
            photons["y"].append((chunk["y"][idxs]-src_ctr[1]).in_units("kpc"))
            photons["z"].append((chunk["z"][idxs]-src_ctr[2]).in_units("kpc"))
            photons["vx"].append(chunk["velocity_x"][idxs].in_units("km/s"))
            photons["vy"].append(chunk["velocity_y"][idxs].in_units("km/s"))
            photons["vz"].append(chunk["velocity_z"][idxs].in_units("km/s"))
            photons["dx"].append(chunk["dx"][idxs].in_units("kpc"))

        pbar.finish()

        for key in photons:
            if len(photons[key]) > 0:
                photons[key] = uconcatenate(photons[key])
            elif key == "NumberOfPhotons":
                photons[key] = np.array([])
            else:
                photons[key] = YTArray([], photon_units[key])

        mylog.info("Number of photons generated: %d" % int(np.sum(photons["NumberOfPhotons"])))
        mylog.info("Number of cells with photons: %d" % len(photons["x"]))

        self.spectral_model.cleanup_spectrum()

        return photons
