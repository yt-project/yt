"""
StarAnalysis - Functions to analyze stars.



"""
from __future__ import print_function

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os
import numpy as np
import h5py
import math
import itertools

from yt.config import ytcfg
from yt.funcs import \
    iterable, get_pbar
from yt.utilities.cosmology import \
    Cosmology
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.physical_constants import \
    sec_per_year, \
    speed_of_light_cgs


YEAR = sec_per_year  # sec / year
LIGHT = speed_of_light_cgs  # cm / s


class StarFormationRate(object):

    r"""Calculates the star formation rate for a given population of
    star particles.

    Parameters
    ----------
    ds : EnzoDataset object
    data_source : AMRRegion object, optional
        The region from which stars are extracted for analysis. If this
        is not supplied, the next three must be, otherwise the next
        three do not need to be specified.
    star_mass : Ordered array or list of floats
        The mass of the stars to be analyzed in units of Msun.
    star_creation_time : Ordered array or list of floats
        The creation time for the stars in code units.
    bins : Integer
        The number of time bins used for binning the stars. Default = 300.
    star_filter : A user-defined filtering rule for stars.
        See: http://yt-project.org/docs/dev/analyzing/filtering.html
        Default: ct>0

    Examples
    --------

    >>> import yt
    >>> from yt.analysis_modules.star_analysis.api import StarFormationRate
    >>> ds = yt.load("Enzo_64/RD0006/RedshiftOutput0006")
    >>> sp = ds.sphere([0.5, 0.5, 0.5], 0.1)
    >>> sfr = StarFormationRate(ds, sp)
    """

    def __init__(self, ds, data_source=None, star_mass=None,
                 star_creation_time=None, bins=300,
                 star_filter=None):
        self._ds = ds
        self._data_source = data_source
        self._filter = star_filter
        self.filter_provided = self._filter is not None
        self.star_mass = np.array(star_mass)
        self.star_creation_time = np.array(star_creation_time)
        self.bin_count = bins
        # Check to make sure we have the right set of informations.
        if data_source is None:
            if self.star_mass is None or self.star_creation_time is None:
                mylog.error(
                    """
                If data_source is not provided, all of these parameters
                need to be set:
                  star_mass (array, Msun),
                  star_creation_time (array, code units),
                """)
                return None
            self.mode = 'provided'
        else:
            self.mode = 'data_source'
        # Set up for time conversion.
        self.cosm = Cosmology(
            hubble_constant=self._ds.hubble_constant,
            omega_matter=self._ds.omega_matter,
            omega_lambda=self._ds.omega_lambda)
        # Find the time right now.
        self.time_now = self._ds.current_time
        # Build the distribution.
        self.build_dist()
        # Attach some convenience arrays.
        self.attach_arrays()

    def build_dist(self):
        """
        Build the data for plotting.
        """
        # Pick out the stars.
        if self.filter_provided:
            ct = self._filter['creation_time']
            mass_stars = self._data_source[self._filter, "particle_mass"]
        else:
            if self.mode == 'data_source':
                ct = self._data_source['creation_time']
                if ct is None:
                    errmsg = 'data source must have particle_age!'
                    mylog.error(errmsg)
                    raise RuntimeError(errmsg)
                mask = ct > 0
                if not any(mask):
                    errmsg = 'all particles have age < 0'
                    mylog.error(errmsg)
                    raise RuntimeError(errmsg)
                # type = self._data_source['particle_type']
                ct_stars = ct[mask]
                mass_stars = self._data_source[
                    'particle_mass'][mask].in_units('Msun')
                del mask
            elif self.mode == 'provided':
                ct_stars = self.star_creation_time
                mass_stars = self.star_mass
        # Find the oldest stars in units of code time.
        tmin = min(ct_stars.in_units("s"))
        # Multiply the end to prevent numerical issues.
        self.time_bins = np.linspace(
            tmin * 1.01, self._ds.current_time.in_units("s"),
            num=self.bin_count + 1)
        # Figure out which bins the stars go into.
        inds = np.digitize(ct_stars.in_units("s"), self.time_bins) - 1
        # Sum up the stars created in each time bin.
        self.mass_bins = np.zeros(self.bin_count + 1, dtype='float64')
        for index in np.unique(inds):
            self.mass_bins[index] += sum(mass_stars[inds == index].tolist())
        # Calculate the cumulative mass sum over time by forward adding.
        self.cum_mass_bins = self.mass_bins.copy()
        for index in xrange(self.bin_count):
            self.cum_mass_bins[index + 1] += self.cum_mass_bins[index]
        # We will want the time taken between bins.
        self.time_bins_dt = self.time_bins[:-1] - self.time_bins[1:]

    def attach_arrays(self):
        """
        Attach convenience arrays to the class for easy access.
        """
        if self.mode == 'data_source':
            try:
                vol = self._data_source[
                    'cell_volume'].in_units('Mpccm**3').sum()
            except AttributeError:
                # If we're here, this is probably a HOPHalo object, and we
                # can get the volume this way.
                ds = self._data_source.get_sphere()
                vol = ds['cell_volume'].in_units('Mpccm**3').sum()
        elif self.mode == 'provided':
            vol = self['cell_volume'].in_units('Mpccm**3').sum()
        # tc = self._ds["Time"] # code time to seconds conversion factor
        self.time = []
        self.lookback_time = []
        self.redshift = []
        self.Msol_yr = []
        self.Msol_yr_vol = []
        self.Msol = []
        self.Msol_cumulative = []
        co = Cosmology()
        # Use the center of the time_bin, not the left edge.
        for i, time in enumerate((self.time_bins[1:] +
                                  self.time_bins[:-1]) / 2.):
            self.time.append(time.in_units('yr'))
            self.lookback_time.append(
                self.time_now.in_units('yr') - time.in_units('yr'))
            self.redshift.append(co.z_from_t(time.in_units('s')))
            self.Msol_yr.append(self.mass_bins[i] /
                                (self.time_bins_dt[i].in_units('yr')))
            # changed vol from mpc to mpccm used in literature
            self.Msol_yr_vol.append(
                self.mass_bins[i] / (self.time_bins_dt[i].in_units('yr')) / vol
            )
            self.Msol.append(self.mass_bins[i])
            self.Msol_cumulative.append(self.cum_mass_bins[i])
        self.time = np.array(self.time)
        self.lookback_time = np.array(self.lookback_time)
        self.redshift = np.array(self.redshift)
        self.Msol_yr = np.array(self.Msol_yr)
        self.Msol_yr_vol = np.array(self.Msol_yr_vol)
        self.Msol = np.array(self.Msol)
        self.Msol_cumulative = np.array(self.Msol_cumulative)

    def write_out(self, name="StarFormationRate.out"):
        r"""Write out the star analysis to a text file *name*. The columns are in
        order.

        The columns in the output file are:
           1. Time (yrs)
           2. Look-back time (yrs)
           #3. Redshift
           4. Star formation rate in this bin per year (Msol/yr)
           5. Star formation rate in this bin per year
              per Mpc**3 (Msol/yr/Mpc**3)
           6. Stars formed in this time bin (Msol)
           7. Cumulative stars formed up to this time bin (Msol)

        Parameters
        ----------
        name : String
            The name of the file to write to. Default = StarFormationRate.out.

        Examples
        --------
        >>> import yt
        >>> from yt.analysis_modules.star_analysis.api import StarFormationRate
        >>> ds = yt.load("Enzo_64/RD0006/RedshiftOutput0006")
        >>> sp = ds.sphere([0.5, 0.5, 0.5], 0.1)
        >>> sfr = StarFormationRate(ds, sp)
        >>> sfr.write_out("stars-SFR.out")
        """
        fp = open(name, "w")
        fp.write(
            "#time\tlookback\tredshift\tMsol/yr\tMsol/yr/Mpc3\tMsol\tcumMsol\t\n")
        for i, time in enumerate(self.time):
            line = "%1.5e %1.5e %1.5e %1.5e %1.5e %1.5e %1.5e\n" % \
                (time,  # Time
                 self.lookback_time[i],  # Lookback time
                 self.redshift[i],  # Redshift
                 self.Msol_yr[i],  # Msol/yr
                 self.Msol_yr_vol[i],  # Msol/yr/vol
                 self.Msol[i],  # Msol in bin
                 self.Msol_cumulative[i])  # cumulative
            fp.write(line)
        fp.close()

# Begin Synthetic Spectrum Stuff.

CHABRIER = {
    "Z0001": "bc2003_hr_m22_chab_ssp.ised.h5",  # /* 0.5% */
    "Z0004": "bc2003_hr_m32_chab_ssp.ised.h5",  # /* 2% */
    "Z004": "bc2003_hr_m42_chab_ssp.ised.h5",  # /* 20% */
    "Z008": "bc2003_hr_m52_chab_ssp.ised.h5",  # /* 40% */
    "Z02": "bc2003_hr_m62_chab_ssp.ised.h5",  # /* solar; 0.02 */
    "Z05": "bc2003_hr_m72_chab_ssp.ised.h5"  # /* 250% */
}

SALPETER = {
    "Z0001": "bc2003_hr_m22_salp_ssp.ised.h5",  # /* 0.5% */
    "Z0004": "bc2003_hr_m32_salp_ssp.ised.h5",  # /* 2% */
    "Z004": "bc2003_hr_m42_salp_ssp.ised.h5",  # /* 20% */
    "Z008": "bc2003_hr_m52_salp_ssp.ised.h5",  # /* 40% */
    "Z02": "bc2003_hr_m62_salp_ssp.ised.h5",  # /* solar; 0.02 */
    "Z05": "bc2003_hr_m72_salp_ssp.ised.h5"  # /* 250% */
}

Zsun = 0.02

# /* dividing line of metallicity; linear in log(Z/Zsun) */
METAL1 = 0.01  # /* in units of Z/Zsun */
METAL2 = 0.0632
METAL3 = 0.2828
METAL4 = 0.6325
METAL5 = 1.5811
METALS = np.array([METAL1, METAL2, METAL3, METAL4, METAL5])

# Translate METALS array digitize to the table dicts
MtoD = np.array(["Z0001", "Z0004", "Z004", "Z008", "Z02",  "Z05"])

"""
This spectrum code is based on code from Ken Nagamine, converted from C to
Python. I've also reversed the order of elements in the flux arrays to be in
C-ordering, for faster memory access."""


class SpectrumBuilder(object):

    r"""Initialize the data to build a summed flux spectrum for a
    collection of stars using the models of Bruzual & Charlot (2003).
    This function loads the necessary data tables into memory and
    must be called before analyzing any star particles.

    Parameters
    ----------
    ds : EnzoDataset object
    bcdir : String
        Path to directory containing Bruzual & Charlot h5 fit files.
    model : String
        Choice of Initial Metalicity Function model, 'chabrier' or
        'salpeter'. Default = 'chabrier'.

    Examples
    --------
    >>> import yt
    >>> from yt.analysis_modules.star_analysis.api import SpectrumBuilder
    >>> ds = yt.load("Enzo_64/RD0006/RedshiftOutput0006")
    >>> spec = SpectrumBuilder(ds, "bc", model="salpeter")
    """

    def __init__(self, ds, bcdir="", model="chabrier", time_now=None,
                 star_filter=None):
        self._ds = ds
        if not os.path.isdir(bcdir):
            bcdir = os.path.join(ytcfg.get("yt", "test_data_dir"), bcdir)
            if not os.path.isdir(bcdir):
                raise RuntimeError("Failed to locate %s" % bcdir)
        self.bcdir = bcdir
        self._filter = star_filter
        self.filter_provided = self._filter is not None
        if model == "chabrier":
            self.model = CHABRIER
        elif model == "salpeter":
            self.model = SALPETER
        # Set up for time conversion.
        self.cosm = Cosmology(
            hubble_constant=self._ds.hubble_constant,
            omega_matter=self._ds.omega_matter,
            omega_lambda=self._ds.omega_lambda)
        # Find the time right now.

        if time_now is None:
            self.time_now = self._ds.current_time
        else:
            self.time_now = time_now

        # Read the tables.
        self.read_bclib()

    def read_bclib(self):
        """
        Read in the age and wavelength bins, and the flux bins for each
        metallicity.
        """
        self.flux = {}
        for file in self.model:
            fname = self.bcdir + "/" + self.model[file]
            fp = h5py.File(fname, 'r')
            self.age = fp["agebins"][:]  # 1D floats
            self.wavelength = fp["wavebins"][:]  # 1D floats
            self.flux[file] = fp["flam"][:, :]  # 2D floats, [agebin, wavebin]
            fp.close()

    def calculate_spectrum(self, data_source=None, star_mass=None,
                           star_creation_time=None,
                           star_metallicity_fraction=None,
                           star_metallicity_constant=None, min_age=0.):
        r"""For the set of stars, calculate the collective spectrum.
        Attached to the output are several useful objects:
        final_spec: The collective spectrum in units of flux binned in
                    wavelength.
        wavelength: The wavelength for the spectrum bins, in Angstroms.
        total_mass: Total mass of all the stars.
        avg_mass: Average mass of all the stars.
        avg_metal: Average metallicity of all the stars.

        Parameters
        ----------
        data_source : AMRRegion object, optional
            The region from which stars are extracted for analysis. If this is
            not specified, the next three parameters must be supplied.
        star_mass : Array or list of floats
            An array of star masses in Msun units.
        star_creation_time : Array or list of floats
            An array of star creation times in code units.
        star_metallicity_fraction : Array or list of floats
            An array of star metallicity fractions, in code
                units (which is not Z/Zsun, rather just Z).
        star_metallicity_constant : Float
            If desired, override the star
            metallicity fraction of all the stars to the given value.
        min_age : Float
            Removes young stars younger than this number (in years)
            from the spectrum. Default: 0 (all stars).

        Examples
        --------
        >>> import yt
        >>> from yt.analysis_modules.star_analysis.api import SpectrumBuilder
        >>> ds = yt.load("Enzo_64/RD0006/RedshiftOutput0006")
        >>> spec = SpectrumBuilder(ds, "bc", model="salpeter")
        >>> sp = ds.sphere([0.5, 0.5, 0.5], 0.1)
        >>> spec.calculate_spectrum(data_source=sp, min_age = 1.e6)
        """
        # Initialize values
        self.final_spec = np.zeros(self.wavelength.size, dtype='float64')
        self._data_source = data_source
        if iterable(star_mass):
            self.star_mass = star_mass
        else:
            self.star_mass = [star_mass]
        if iterable(star_creation_time):
            self.star_creation_time = star_creation_time
        else:
            self.star_creation_time = [star_creation_time]
        if iterable(star_metallicity_fraction):
            self.star_metal = star_metallicity_fraction
        else:
            self.star_metal = [star_metallicity_fraction]
        self.min_age = min_age

        # Check to make sure we have the right set of data.
        if data_source is None:
            if self.star_mass is None or self.star_creation_time is None or \
                    (star_metallicity_fraction is None and
                     star_metallicity_constant is None):
                mylog.error(
                    """
                If data_source is not provided, all of these paramters
                need to be set:
                   star_mass (array, Msun),
                   star_creation_time (array, code units),
                And one of:
                   star_metallicity_fraction (array, code units).
                --OR--
                   star_metallicity_constant (float, code units).
                """)
                return None
            if star_metallicity_constant is not None:
                self.star_metal = np.ones(self.star_mass.size, dtype='float64') * \
                    star_metallicity_constant
            if star_metallicity_fraction is not None:
                self.star_metal = star_metallicity_fraction
        else:
            # Get the data we need.
            if self.filter_provided:
                ct = self._filter['creation_time']
                # mass_stars = self._data_source[self._filter, "particle_mass"]
                if star_metallicity_constant is not None:
                    self.star_metal = np.ones(self.star_mass.size, dtype='float64') * \
                        star_metallicity_constant
                else:
                    self.star_metal = self._data_source[
                        self._filter, "metallicity_fraction"].in_units('Zsun')
            else:
                ct = self._data_source["creation_time"]
                if ct is None:
                    errmsg = 'data source must have particle_age!'
                    mylog.error(errmsg)
                    raise RuntimeError(errmsg)
                mask = ct > 0
                if not any(mask):
                    errmsg = 'all particles have age < 0'
                    mylog.error(errmsg)
                    raise RuntimeError(errmsg)
                # type = self._data_source['particle_type']
                self.star_creation_time = ct[mask]
                self.star_mass = self._data_source[
                    'particle_mass'][mask].in_units('Msun')
                if star_metallicity_constant is not None:
                    self.star_metal = np.ones(self.star_mass.size, dtype='float64') * \
                        star_metallicity_constant
                else:
                    self.star_metal = self._data_source[
                        "metallicity_fraction"][mask].in_units('Zsun')
        # Age of star in years.
        dt = (self.time_now - self.star_creation_time).in_units('yr').tolist()
        dt = np.maximum(dt, 0.0)
        # Remove young stars
        sub = dt >= self.min_age
        if len(sub) == 0:
            return
        self.star_metal = self.star_metal[sub]
        dt = dt[sub]
        self.star_creation_time = self.star_creation_time[sub]
        # Figure out which METALS bin the star goes into.
        Mindex = np.digitize(self.star_metal, METALS)
        # Replace the indices with strings.
        Mname = MtoD[Mindex]
        # Figure out which age bin this star goes into.
        Aindex = np.digitize(dt, self.age)
        # Ratios used for the interpolation.
        ratio1 = (dt - self.age[Aindex - 1]) / \
            (self.age[Aindex] - self.age[Aindex - 1])
        ratio2 = (self.age[Aindex] - dt) / \
            (self.age[Aindex] - self.age[Aindex - 1])
        # Sort the stars by metallicity and then by age, which should reduce
        # memory access time by a little bit in the loop.
        indexes = np.arange(self.star_metal.size)
        sort = np.asarray([indexes[i]
                           for i in np.lexsort([indexes, Aindex, Mname])])
        Mname = Mname[sort]
        Aindex = Aindex[sort]
        ratio1 = ratio1[sort]
        ratio2 = ratio2[sort]
        self.star_mass = self.star_mass[sort]
        self.star_creation_time = self.star_creation_time[sort]
        self.star_metal = self.star_metal[sort]

        # Interpolate the flux for each star, adding to the total by weight.
        pbar = get_pbar("Calculating fluxes", len(self.star_mass))
        for i, star in enumerate(itertools.izip(Mname, Aindex, ratio1, ratio2,
                                                self.star_mass)):
            # Pick the right age bin for the right flux array.
            flux = self.flux[star[0]][star[1], :]
            # Get the one just before the one above.
            flux_1 = self.flux[star[0]][star[1] - 1, :]
            # interpolate in log(flux), linear in time.
            int_flux = star[3] * np.log10(flux_1) + star[2] * np.log10(flux)
            # Add this flux to the total, weighted by mass.
            self.final_spec += np.power(10., int_flux) * star[4]
            pbar.update(i)
        pbar.finish()

        # Normalize.
        self.total_mass = self.star_mass.sum()
        self.avg_mass = self.star_mass.mean()
        tot_metal = (self.star_metal * self.star_mass).sum()
        if tot_metal > 0:
            self.avg_metal = math.log10(
                (tot_metal / self.total_mass).in_units('Zsun'))
        else:
            self.avg_metal = -99

        # Below is an attempt to do the loop using vectors and matrices,
        # however it doesn't appear to be much faster, probably due to all
        # the gymnastics that have to be done to do element-by-element
        # multiplication for matricies.
        # I'm keeping it in here in case I come up with a more
        # elegant way that actually is faster.
#         for metal_name in MtoD:
#             # Pick out our stars in this metallicity bin.
#             select = (Mname == metal_name)
#             A = Aindex[select]
#             if A.size == 0: continue
#             r1 = ratio1[select]
#             r2 = ratio2[select]
#             sm = self.star_mass[select]
#             # From the flux array for this metal, and our selection, build
#             # a new flux array just for the ages of these stars, in the
#             # same order as the selection of stars.
#             this_flux = np.matrix(self.flux[metal_name][A])
#             # Make one for the last time step for each star in the same fashion
#             # as above.
#             this_flux_1 = np.matrix(self.flux[metal_name][A-1])
#             # This is kind of messy, but we're going to multiply this_fluxes
#             # by the appropriate ratios and add it together to do the
#             # interpolation in log(flux) and linear in time.
#             print r1.size
#             r1 = np.matrix(r1.tolist() * self.wavelength.size).reshape(self.wavelength.size,r1.size).T
#             r2 = np.matrix(r2.tolist() * self.wavelength.size).reshape(self.wavelength.size,r2.size).T
#             print this_flux_1.shape, r1.shape
#             int_flux = np.multiply(np.log10(this_flux_1),r1) \
#                 + np.multiply(np.log10(this_flux),r2)
#             # Weight the fluxes by mass.
#             sm = np.matrix(sm.tolist()*self.wavelength.size).reshape(self.wavelength.size,sm.size).T
#             int_flux = np.multiply(np.power(10., int_flux), sm)
#             # Sum along the columns, converting back to an array, adding
#             # to the full spectrum.
#             self.final_spec += np.array(int_flux.sum(axis=0))[0,:]

    def write_out(self, name="sum_flux.out"):
        r"""Write out the summed flux to a file.

        The output file from this function has two columns: Wavelength
        (Angstrom) and Flux (Luminosity per unit wavelength, L_sun Ang^-1,
        L_sun = 3.826 * 10^33 ergs s^-1.).

        Parameters
        ----------
        name : String
            Name of file to write to. Default = "sum_flux.out"

        Examples
        --------
        >>> import yt
        >>> from yt.analysis_modules.star_analysis.api import SpectrumBuilder
        >>> ds = yt.load("Enzo_64/RD0006/RedshiftOutput0006")
        >>> sp = ds.sphere([0.5, 0.5, 0.5], 0.1)
        >>> spec = SpectrumBuilder(ds, "bc", model="salpeter")
        >>> spec.calculate_spectrum(data_source=sp, min_age = 1.e6)
        >>> spec.write_out("spec.out")
        """
        fp = open(name, 'w')
        for i, wave in enumerate(self.wavelength):
            fp.write("%1.5e\t%1.5e\n" % (wave, self.final_spec[i]))
        fp.close()

    def write_out_SED(self, name="sum_SED.out", flux_norm=5200.):
        r"""Write out the summed SED to a file. The file has two columns:
        1) Wavelength (Angstrom)
        2) Relative flux normalized to the flux at *flux_norm*.
        It also will attach to the SpectrumBuilder object
        an array *f_nu* which is the normalized flux,
        identical to the disk output.

        Parameters
        ----------
        name : String
            Name of file to write to. Default = "sum_SED.out"
        flux_norm : Float
            Wavelength of the flux to normalize the distribution against.
            Default = 5200 Ang.

        Examples
        --------
        >>> import yt
        >>> from yt.analysis_modules.star_analysis.api import SpectrumBuilder
        >>> ds = yt.load("Enzo_64/RD0006/RedshiftOutput0006")
        >>> spec = SpectrumBuilder(ds, "bc", model="salpeter")
        >>> sp = ds.sphere([0.5, 0.5, 0.5], 0.1)
        >>> spec.calculate_spectrum(data_source=sp, min_age = 1.e6)
        >>> spec.write_out_SED(name = "SED.out", flux_norm = 6000.)
        """
        # find the f_nu closest to flux_norm
        fn_wavelength = np.argmin(abs(self.wavelength - flux_norm))
        f_nu = self.final_spec * np.power(self.wavelength, 2.) / LIGHT
        # Normalize f_nu
        self.f_nu = f_nu / f_nu[fn_wavelength]
        # Write out.
        fp = open(name, 'w')
        for i, wave in enumerate(self.wavelength):
            fp.write("%1.5e\t%1.5e\n" % (wave, self.f_nu[i]))
        fp.close()
