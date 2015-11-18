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
from yt.utilities.on_demand_imports import _h5py as h5py
import math

from yt.config import ytcfg
from yt.extern.six.moves import zip as izip
from yt.funcs import \
    get_pbar
from yt.units import \
    g, s, Zsun
from yt.units.yt_array import YTArray, YTQuantity
from yt.utilities.cosmology import \
    Cosmology
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.physical_constants import \
    speed_of_light_cgs


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
    volume : Float
        The comoving volume of the region for the specified list of stars.
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
                 star_creation_time=None, bins=300, volume=None,
                 star_filter=None):
        self._ds = ds
        self._data_source = data_source
        self._filter = star_filter
        self.ds_provided = self._data_source is not None
        self.filter_provided = self._filter is not None
        self.bin_count = bins

        # Set up for time conversion.
        self.cosm = Cosmology(
            hubble_constant=self._ds.hubble_constant,
            omega_matter=self._ds.omega_matter,
            omega_lambda=self._ds.omega_lambda)
        # Find the time right now.
        self.time_now = self._ds.current_time

        if not self.ds_provided:
            # Check to make sure we have the right set of informations.
            if star_mass is None or star_creation_time is None \
                    or volume is None:
                mylog.error("""
                    If data_source is not provided, all of these parameters
                    need to be set:
                        star_mass (array, Msun),
                        star_creation_time (array, code units),
                        volume (float, cMpc**3).""")
                return None

            if isinstance(star_mass, YTArray):
                assert star_mass.units.same_dimensions_as(g.units)
            elif star_mass is not None:
                star_mass = YTArray(star_mass, 'Msun')
            self.star_mass = star_mass

            if isinstance(star_creation_time, YTArray):
                assert star_creation_time.units.same_dimensions_as(s.units)
            elif star_creation_time is not None:
                star_creation_time = self._ds.arr(star_creation_time,
                                                  'code_time')
            self.star_creation_time = star_creation_time

            if isinstance(volume, YTQuantity):
                assert volume.units.same_dimensions_as(
                    self._ds.quan(1.0, 'Mpccm**3').units
                )
            elif volume is not None:
                volume = self._ds.quan(volume, 'Mpccm**3')
            self.volume = volume

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
            if self.ds_provided:
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
            else:
                ct_stars = self.star_creation_time
                mass_stars = self.star_mass
        # Find the oldest stars in units of code time.
        tmin = ct_stars.min().in_units("s")
        # Multiply the end to prevent numerical issues.
        self.time_bins = np.linspace(
            tmin * 1.01, self._ds.current_time.in_units("s"),
            num=self.bin_count + 1)
        # Figure out which bins the stars go into.
        inds = np.digitize(ct_stars.in_units("s"), self.time_bins) - 1
        # Sum up the stars created in each time bin.
        self.mass_bins = YTArray(
            np.zeros(self.bin_count + 1, dtype='float64'), "Msun"
        )
        for index in np.unique(inds):
            self.mass_bins[index] += (mass_stars[inds == index]).sum()
        # We will want the time taken between bins.
        self.time_bins_dt = self.time_bins[1:] - self.time_bins[:-1]

    def attach_arrays(self):
        """
        Attach convenience arrays to the class for easy access.
        """
        if self.ds_provided:
            try:
                vol = self._data_source[
                    'cell_volume'].in_units('Mpccm ** 3').sum()
            except AttributeError:
                # If we're here, this is probably a HOPHalo object, and we
                # can get the volume this way.
                ds = self._data_source.get_sphere()
                vol = ds['cell_volume'].in_units('Mpccm ** 3').sum()
        else:
            vol = self.volume.in_units('Mpccm ** 3')

        # Use the center of the time_bin, not the left edge.
        self.time = 0.5 * \
            (self.time_bins[1:] + self.time_bins[:-1]).in_units('yr')
        self.lookback_time = self.time_now - self.time  # now in code_time...
        self.redshift = self.cosm.z_from_t(self.time)

        self.Msol_yr = (
            self.mass_bins[:-1] / self.time_bins_dt[:]).in_units('Msun/yr')
        # changed vol from mpc to mpccm used in literature
        self.Msol_yr_vol = self.Msol_yr / vol

        self.Msol = self.mass_bins[:-1].in_units("Msun")
        self.Msol_cumulative = self.Msol.cumsum()

    def write_out(self, name="StarFormationRate.out"):
        r"""Write out the star analysis to a text file *name*. The columns are in
        order.

        The columns in the output file are:
           1. Time (yrs)
           2. Look-back time (yrs)
           3. Redshift
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
                (time.in_units("yr"),  # Time
                 self.lookback_time[i].in_units('yr'),  # Lookback time
                 self.redshift[i],  # Redshift
                 self.Msol_yr[i].in_units("Msun/yr"),
                 self.Msol_yr_vol[i],
                 self.Msol[i].in_units("Msun"),  # Msol in bin
                 self.Msol_cumulative[i].in_units("Msun"))  # cumulative
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
            self.age = YTArray(fp["agebins"][:], 'yr')  # 1D floats
            self.wavelength = fp["wavebins"][:]  # 1D floats
            self.flux[file] = fp["flam"][:, :]  # 2D floats, [agebin, wavebin]
            fp.close()

    def calculate_spectrum(self, data_source=None, star_mass=None,
                           star_creation_time=None,
                           star_metallicity_fraction=None,
                           star_metallicity_constant=None,
                           min_age=YTQuantity(0.0, 'yr')):

        r"""For the set of stars, calculate the collective spectrum.
        Attached to the output are several useful objects:

        Attributes
        ----------
        final_spec: array
            The collective spectrum in units of flux binned in wavelength.
        wavelength: array
            The wavelength for the spectrum bins, in Angstroms.
        total_mass: float
            Total mass of all the stars.
        avg_mass: float
            Average mass of all the stars.
        avg_metal: float
            Average metallicity of all the stars.

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
        >>> spec.calculate_spectrum(data_source=sp, min_age=1.e6)
        """

        # Initialize values
        self.final_spec = np.zeros(self.wavelength.size, dtype='float64')
        self._data_source = data_source

        if isinstance(star_mass, YTArray):
            assert star_mass.units.same_dimensions_as(g.units)
        elif star_mass is not None:
            star_mass = YTArray(star_mass, 'Msun')
        self.star_mass = star_mass

        if isinstance(star_creation_time, YTArray):
            assert star_creation_time.units.same_dimensions_as(s.units)
        elif star_creation_time is not None:
            star_creation_time = self._ds.arr(star_creation_time,
                                              'code_time')
        self.star_creation_time = star_creation_time

        if isinstance(star_metallicity_fraction, YTArray):
            assert \
                star_metallicity_fraction.units.same_dimensions_as(Zsun.units)
        elif star_metallicity_fraction is not None:
            star_metallicity_fraction = self._ds.arr(
                star_metallicity_fraction, 'code_metallicity'
            )
        self.star_metallicity_fraction = star_metallicity_fraction

        if isinstance(min_age, YTQuantity):
            assert min_age.units.same_dimensions_as(s.units)
        elif min_age is not None:
            min_age = YTQuantity(min_age, 'yr')
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

            if star_metallicity_fraction is not None:
                self.star_metal = star_metallicity_fraction
            else:
                self.star_metal = \
                    self._ds.arr(np.ones_like(self.star_mass) *
                                 star_metallicity_constant, 'Zsun')
        else:
            # Get the data we need.
            if self.filter_provided:
                ct = self._filter['creation_time']
                # mass_stars = self._data_source[self._filter, "particle_mass"]
                if star_metallicity_constant is None:
                    self.star_metal = self._data_source[
                        self._filter, "metallicity_fraction"].in_units('Zsun')
                else:
                    self.star_metal = \
                        self._ds.arr(np.ones_like(
                            self._data_source[self._filter,
                                              "metallicity_fraction"]) *
                        star_metallicity_constant, "Zsun")
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
                    self.star_metal = self._ds.arr(
                        np.ones_like(self.star_mass) *
                        star_metallicity_constant, 'Zsun')
                else:
                    self.star_metal = self._data_source[
                        "metallicity_fraction"][mask].in_units('Zsun')
        # Age of star in years.
        dt = (self.time_now - self.star_creation_time).in_units('yr')
        dt[dt < 0.0] = 0.0
        # Remove young stars
        sub = dt >= self.min_age
        if len(sub) == 0:
            return
        self.star_metal = self.star_metal[sub]
        dt = dt[sub]
        self.star_creation_time = self.star_creation_time[sub]
        # Figure out which METALS bin the star goes into.
        Mindex = np.digitize(self.star_metal.in_units('Zsun'), METALS)
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
        for i, star in enumerate(izip(Mname, Aindex, ratio1, ratio2,
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
        f_nu = self.final_spec * np.power(self.wavelength, 2.) \
            / speed_of_light_cgs
        # Normalize f_nu
        self.f_nu = f_nu / f_nu[fn_wavelength]
        # Write out.
        fp = open(name, 'w')
        for i, wave in enumerate(self.wavelength):
            fp.write("%1.5e\t%1.5e\n" % (wave, self.f_nu[i]))
        fp.close()
