"""
Projection class for the Sunyaev-Zeldovich effect. Requires SZpack (at least
version 1.1.1) to be downloaded and installed:

http://www.chluba.de/SZpack/

For details on the computations involved please refer to the following references:

Chluba, Nagai, Sazonov, Nelson, MNRAS, 2012, arXiv:1205.5778
Chluba, Switzer, Nagai, Nelson, MNRAS, 2012, arXiv:1211.3206
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.physical_constants import sigma_thompson, clight, hcgs, kboltz, mh, Tcmb
from yt.utilities.fits_image import FITSImageBuffer
from yt.fields.local_fields import add_field
from yt.fields.field_plugin_registry import register_field_plugin
from yt.data_objects.image_array import ImageArray
from yt.funcs import fix_axis, mylog, iterable, get_pbar
from yt.utilities.definitions import inv_axis_names
from yt.visualization.volume_rendering.camera import off_axis_projection
from yt.utilities.parallel_tools.parallel_analysis_interface import \
     communication_system, parallel_root_only
from yt.visualization.plot_window import StandardCenter
from yt import units

import numpy as np

I0 = (2*(kboltz*Tcmb)**3/((hcgs*clight)**2)/units.sr).in_units("MJy/steradian")

try:
    import SZpack
except ImportError:
    pass

vlist = "xyz"

@register_field_plugin
def setup_SZ_fields(registry, ftype = "gas", slice_info = None):

    def _tempkev(field, data):
        return (kboltz*data["gas","temperature"]).in_units("keV")
    registry.add_field(("gas","TempkeV"), function=_tempkev, units="keV")

    def _t_squared(field, data):
        return data["gas","density"]*data["gas","TempkeV"]*data["gas","TempkeV"]
    registry.add_field(("gas","TSquared"), function=_t_squared, units="g*keV**2/cm**3")

    def _beta_perp_squared(field, data):
        return data["gas","density"]*data["gas","velocity_magnitude"]**2/clight/clight - data["gas","BetaParSquared"]
    registry.add_field(("gas","BetaPerpSquared"), function=_beta_perp_squared, units="g/cm**3")

    def _beta_par_squared(field, data):
        return data["gas","BetaPar"]**2/data["gas","density"]
    registry.add_field(("gas","BetaParSquared"), function=_beta_par_squared, units="cm**3/g")

    def _t_beta_par(field, data):
        return data["gas","TempkeV"]*data["gas","BetaPar"]
    registry.add_field(("gas","TBetaPar"), function=_t_beta_par, units="keV*g/cm**3")

    def _t_sz(field, data):
        return data["gas","density"]*data["gas","TempkeV"]
    registry.add_field(("gas","TeSZ"), function=_t_sz, units="keV*g/cm**3")

class SZProjection(object):
    r""" Initialize a SZProjection object.

    Parameters
    ----------
    pf : parameter_file
        The parameter file.
    freqs : array_like
        The frequencies (in GHz) at which to compute the SZ spectral distortion.
    mue : float, optional
        Mean molecular weight for determining the electron number density.
    high_order : boolean, optional
        Should we calculate high-order moments of velocity and temperature?

    Examples
    --------
    >>> freqs = [90., 180., 240.]
    >>> szprj = SZProjection(pf, freqs, high_order=True)
    """
    def __init__(self, pf, freqs, mue=1.143, high_order=False):

        self.pf = pf
        self.num_freqs = len(freqs)
        self.high_order = high_order
        self.freqs = np.array(freqs)
        self.mueinv = 1./mue
        self.xinit = hcgs*pf.arr(self.freqs, "GHz")/(kboltz*Tcmb)
        self.freq_fields = ["%d_GHz" % (int(freq)) for freq in freqs]
        self.data = {}

        self.display_names = {}
        self.display_names["TeSZ"] = r"$\mathrm{T_e}$"
        self.display_names["Tau"] = r"$\mathrm{\tau}$"

        for f, field in zip(self.freqs, self.freq_fields):
            self.display_names[field] = r"$\mathrm{\Delta{I}_{%d\ GHz}}$" % (int(f))

    def on_axis(self, axis, center="c", width=(1, "unitary"), nx=800, source=None):
        r""" Make an on-axis projection of the SZ signal.

        Parameters
        ----------
        axis : integer or string
            The axis of the simulation domain along which to make the SZprojection.
        center : array_like or string, optional
            The center of the projection.
        width : float or tuple
            The width of the projection.
        nx : integer, optional
            The dimensions on a side of the projection image.
        source : yt.data_objects.api.AMRData, optional
            If specified, this will be the data source used for selecting regions to project.

        Examples
        --------
        >>> szprj.on_axis("y", center="max", width=(1.0, "mpc"), source=my_sphere)
        """
        axis = fix_axis(axis)

        if center == "c":
            ctr = self.pf.domain_center
        elif center == "max":
            v, ctr = self.pf.h.find_max("density")
        else:
            ctr = center

        def _beta_par(field, data):
            axis = data.get_field_parameter("axis")
            # Load these, even though we will only use one
            for ax in 'xyz':
                data['velocity_%s' % ax]
            vpar = data["density"]*data["velocity_%s" % (vlist[axis])]
            return vpar/clight
        add_field(("gas","BetaPar"), function=_beta_par, units="g/cm**3")

        proj = self.pf.h.proj("density", axis, center=ctr, data_source=source)
        proj.data_source.set_field_parameter("axis", axis)
        frb = proj.to_frb(width, nx)
        dens = frb["density"]
        Te = frb["TeSZ"]/dens
        bpar = frb["BetaPar"]/dens
        omega1 = frb["TSquared"]/dens/(Te*Te) - 1.
        bperp2 = np.zeros((nx,nx))
        sigma1 = np.zeros((nx,nx))
        kappa1 = np.zeros((nx,nx))
        if self.high_order:
            bperp2 = frb["BetaPerpSquared"]/dens
            sigma1 = frb["TBetaPar"]/dens/Te - bpar
            kappa1 = frb["BetaParSquared"]/dens - bpar*bpar
        tau = sigma_thompson*dens*self.mueinv/mh

        nx,ny = frb.buff_size
        self.bounds = frb.bounds
        self.dx = (frb.bounds[1]-frb.bounds[0])/nx
        self.dy = (frb.bounds[3]-frb.bounds[2])/ny
        self.nx = nx

        self._compute_intensity(tau, Te, bpar, omega1, sigma1, kappa1, bperp2)

    def off_axis(self, L, center="c", width=(1, "unitary"), nx=800, source=None):
        r""" Make an off-axis projection of the SZ signal.

        Parameters
        ----------
        L : array_like
            The normal vector of the projection.
        center : array_like or string, optional
            The center of the projection.
        width : float or tuple
            The width of the projection.
        nx : integer, optional
            The dimensions on a side of the projection image.
        source : yt.data_objects.api.AMRData, optional
            If specified, this will be the data source used for selecting regions to project.
            Currently unsupported in yt 2.x.

        Examples
        --------
        >>> L = np.array([0.5, 1.0, 0.75])
        >>> szprj.off_axis(L, center="c", width=(2.0, "mpc"))
        """
        if iterable(width):
            w = width[0]/self.pf.units[width[1]]
        else:
            w = width
        if center == "c":
            ctr = self.pf.domain_center
        elif center == "max":
            v, ctr = self.pf.h.find_max("density")
        else:
            ctr = center

        if source is not None:
            mylog.error("Source argument is not currently supported for off-axis S-Z projections.")
            raise NotImplementedError

        def _beta_par(field, data):
            vpar = data["density"]*(data["velocity_x"]*L[0]+
                                    data["velocity_y"]*L[1]+
                                    data["velocity_z"]*L[2])
            return vpar/clight
        add_field("BetaPar", function=_beta_par,

        dens    = off_axis_projection(self.pf, ctr, L, w, nx, "density")
        Te      = off_axis_projection(self.pf, ctr, L, w, nx, "TeSZ")/dens
        bpar    = off_axis_projection(self.pf, ctr, L, w, nx, "BetaPar")/dens
        omega1  = off_axis_projection(self.pf, ctr, L, w, nx, "TSquared")/dens
        omega1  = omega1/(Te*Te) - 1.
        if self.high_order:
            bperp2  = off_axis_projection(self.pf, ctr, L, w, nx, "BetaPerpSquared")/dens
            sigma1  = off_axis_projection(self.pf, ctr, L, w, nx, "TBetaPar")/dens
            sigma1  = sigma1/Te - bpar
            kappa1  = off_axis_projection(self.pf, ctr, L, w, nx, "BetaParSquared")/dens
            kappa1 -= bpar
        else:
            bperp2 = np.zeros((nx,nx))
            sigma1 = np.zeros((nx,nx))
            kappa1 = np.zeros((nx,nx))
        tau = sigma_thompson*dens*self.mueinv/mh

        self.bounds = np.array([-0.5*w, 0.5*w, -0.5*w, 0.5*w])
        self.dx = w/nx
        self.dy = w/nx
        self.nx = nx

        self._compute_intensity(tau, Te, bpar, omega1, sigma1, kappa1, bperp2)

    def _compute_intensity(self, tau, Te, bpar, omega1, sigma1, kappa1, bperp2):

        # Bad hack, but we get NaNs if we don't do something like this
        small_beta = np.abs(bpar) < 1.0e-20
        bpar[small_beta] = 1.0e-20

        comm = communication_system.communicators[-1]

        nx, ny = self.nx,self.nx
        signal = np.zeros((self.num_freqs,nx,ny))
        xo = np.zeros((self.num_freqs))

        k = int(0)

        start_i = comm.rank*nx/comm.size
        end_i = (comm.rank+1)*nx/comm.size

        pbar = get_pbar("Computing SZ signal.", nx*nx)

        for i in xrange(start_i, end_i):
            for j in xrange(ny):
                xo[:] = self.xinit[:]
                SZpack.compute_combo_means(xo, tau[i,j], Te[i,j],
                                           bpar[i,j], omega1[i,j],
                                           sigma1[i,j], kappa1[i,j], bperp2[i,j])
                signal[:,i,j] = xo[:]
                pbar.update(k)
                k += 1

        signal = comm.mpi_allreduce(signal)

        pbar.finish()

        for i, field in enumerate(self.freq_fields):
            self.data[field] = ImageArray(I0*self.xinit[i]**3*signal[i,:,:])
        self.data["Tau"] = ImageArray(tau)
        self.data["TeSZ"] = ImageArray(Te)

    @parallel_root_only
    def write_fits(self, filename, sky_center=None, sky_scale=None, clobber=True):
        r""" Export images to a FITS file. Writes the SZ distortion in all
        specified frequencies as well as the mass-weighted temperature and the
        optical depth. Distance units are in kpc, unless *sky_center*
        and *scale* are specified. 

        Parameters
        ----------
        filename : string
            The name of the FITS file to be written. 
        sky_center : tuple of floats, optional
            The center of the observation in (RA, Dec) in degrees. Only used if
            converting to sky coordinates.          
        sky_scale : float, optional
            Scale between degrees and kpc. Only used if
            converting to sky coordinates.
        clobber : boolean, optional
            If the file already exists, do we overwrite?

        Examples
        --------
        >>> # This example just writes out a FITS file with kpc coords
        >>> szprj.write_fits("SZbullet.fits", clobber=False)
        >>> # This example uses sky coords
        >>> sky_scale = 1./3600. # One arcsec per kpc
        >>> sky_center = (30., 45.) # In degrees
        >>> szprj.write_fits("SZbullet.fits", sky_center=sky_center, sky_scale=sky_scale)
        """

        deltas = np.array([self.dx*self.pf.units["kpc"],
                           self.dy*self.pf.units["kpc"]])

        if sky_center is None:
            center = [0.0]*2
            units = "kpc"
        else:
            center = sky_center
            units = "deg"
            deltas *= sky_scale
            
        fib = FITSImageBuffer(self.data, fields=self.data.keys(),
                              center=center, units=units,
                              scale=deltas)
        fib.update_all_headers("Time", self.pf.current_time)
        fib.writeto(filename, clobber=clobber)
        
    @parallel_root_only
    def write_png(self, filename_prefix, cmap_name="algae",
                  log_fields=None):
        r""" Export images to PNG files. Writes the SZ distortion in all
        specified frequencies as well as the mass-weighted temperature and the
        optical depth. Distance units are in kpc.

        Parameters
        ----------
        filename_prefix : string
            The prefix of the image filenames.

        Examples
        --------
        >>> szprj.write_png("SZsloshing")
        """
        import matplotlib
        import matplotlib.pyplot as plt
        if log_fields is None: log_fields = {}
        ticks_font = matplotlib.font_manager.FontProperties(family='serif',size=16)
        extent = tuple([bound*self.pf.units["kpc"] for bound in self.bounds])
        for field, image in self.items():
            data = image.copy()
            vmin, vmax = image.min(), image.max()
            negative = False
            crossover = False
            if vmin < 0 and vmax < 0:
                data *= -1
                negative = True                                        
            if log_fields.has_key(field):
                log_field = log_fields[field]
            else:
                log_field = True
            if log_field:
                formatter = matplotlib.ticker.LogFormatterMathtext()        
                norm = matplotlib.colors.LogNorm()
                if vmin < 0 and vmax > 0:
                    crossover = True
                    linthresh = min(vmax, -vmin)/100.
                    norm=matplotlib.colors.SymLogNorm(linthresh,
                                                      vmin=vmin, vmax=vmax)
            else:
                norm = None
                formatter = None
            filename = filename_prefix+"_"+field+".png"
            cbar_label = self.display_names[field]
            if self.units[field] is not None:
                cbar_label += " ("+self.units[field]+")"
            fig = plt.figure(figsize=(10.0,8.0))
            ax = fig.add_subplot(111)
            cax = ax.imshow(data, norm=norm, extent=extent, cmap=cmap_name, origin="lower")
            for label in ax.get_xticklabels():
                label.set_fontproperties(ticks_font)
            for label in ax.get_yticklabels():
                label.set_fontproperties(ticks_font)                      
            ax.set_xlabel(r"$\mathrm{x\ (kpc)}$", fontsize=16)
            ax.set_ylabel(r"$\mathrm{y\ (kpc)}$", fontsize=16)
            cbar = fig.colorbar(cax, format=formatter)
            cbar.ax.set_ylabel(cbar_label, fontsize=16)
            if negative:
                cbar.ax.set_yticklabels(["-"+label.get_text()
                                         for label in cbar.ax.get_yticklabels()])
            if crossover:
                yticks = list(-10**np.arange(np.floor(np.log10(-vmin)),
                                             np.rint(np.log10(linthresh))-1, -1)) + [0] + \
                         list(10**np.arange(np.rint(np.log10(linthresh)),
                                            np.ceil(np.log10(vmax))+1))
                cbar.set_ticks(yticks)
            for label in cbar.ax.get_yticklabels():
                label.set_fontproperties(ticks_font)                 
            fig.tight_layout()
            plt.savefig(filename)
            
    @parallel_root_only
    def write_hdf5(self, filename):
        r"""Export the set of S-Z fields to a set of HDF5 datasets.

        Parameters
        ----------
        filename : string
            This file will be opened in "write" mode.

        Examples
        --------
        >>> szprj.write_hdf5("SZsloshing.h5")
        """
        import h5py
        f = h5py.File(filename, "w")
        for field, data in self.items():
            f.create_dataset(field,data=data.ndarray_view())
        f.close()

    def keys(self):
        return self.data.keys()

    def items(self):
        return self.data.items()

    def values(self):
        return self.data.values()

    def has_key(self, key):
        return key in self.data.keys()

    def __getitem__(self, key):
        return self.data[key]

    @property
    def shape(self):
        return (self.nx,self.nx)
