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

from yt.config import \
    ytcfg
from yt.utilities.physical_constants import sigma_thompson, clight, hcgs, kboltz, mh, Tcmb
from yt.funcs import fix_axis, get_pbar
from yt.visualization.volume_rendering.off_axis_projection import \
    off_axis_projection
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    communication_system, parallel_root_only
from yt import units
from yt.utilities.on_demand_imports import _astropy

import numpy as np

I0 = (2*(kboltz*Tcmb)**3/((hcgs*clight)**2)/units.sr).in_units("MJy/steradian")

try:
    import SZpack
except ImportError:
    pass

vlist = "xyz"
def setup_sunyaev_zeldovich_fields(ds):
    def _t_squared(field, data):
        return data["gas","density"]*data["gas","kT"]*data["gas","kT"]
    ds.add_field(("gas", "t_squared"), function = _t_squared,
                 units="g*keV**2/cm**3")

    def _beta_par_squared(field, data):
        return data["gas","beta_par"]**2/data["gas","density"]
    ds.add_field(("gas","beta_par_squared"), function = _beta_par_squared,
                 units="g/cm**3")

    def _beta_perp_squared(field, data):
        return data["gas","density"]*data["gas","velocity_magnitude"]**2/clight/clight - data["gas","beta_par_squared"]
    ds.add_field(("gas","beta_perp_squared"), function = _beta_perp_squared,
                 units="g/cm**3")

    def _t_beta_par(field, data):
        return data["gas","kT"]*data["gas","beta_par"]
    ds.add_field(("gas","t_beta_par"), function = _t_beta_par,
                 units="keV*g/cm**3")

    def _t_sz(field, data):
        return data["gas","density"]*data["gas","kT"]
    ds.add_field(("gas","t_sz"), function = _t_sz,
                 units="keV*g/cm**3")

def generate_beta_par(L):
    def _beta_par(field, data):
        vpar = data["density"]*(data["velocity_x"]*L[0]+
                                data["velocity_y"]*L[1]+
                                data["velocity_z"]*L[2])
        return vpar/clight
    return _beta_par

class SZProjection(object):
    r""" Initialize a SZProjection object.

    Parameters
    ----------
    ds : Dataset
        The dataset
    freqs : array_like
        The frequencies (in GHz) at which to compute the SZ spectral distortion.
    mue : float, optional
        Mean molecular weight for determining the electron number density.
    high_order : boolean, optional
        Should we calculate high-order moments of velocity and temperature?

    Examples
    --------
    >>> freqs = [90., 180., 240.]
    >>> szprj = SZProjection(ds, freqs, high_order=True)
    """
    def __init__(self, ds, freqs, mue=1.143, high_order=False):

        self.ds = ds
        self.num_freqs = len(freqs)
        self.high_order = high_order
        self.freqs = ds.arr(freqs, "GHz")
        self.mueinv = 1./mue
        self.xinit = hcgs*self.freqs.in_units("Hz")/(kboltz*Tcmb)
        self.freq_fields = ["%d_GHz" % (int(freq)) for freq in freqs]
        self.data = {}

        self.display_names = {}
        self.display_names["TeSZ"] = r"$\mathrm{T_e}$"
        self.display_names["Tau"] = r"$\mathrm{\tau}$"

        for f, field in zip(self.freqs, self.freq_fields):
            self.display_names[field] = r"$\mathrm{\Delta{I}_{%d\ GHz}}$" % int(f)

    def on_axis(self, axis, center="c", width=(1, "unitary"), nx=800, source=None):
        r""" Make an on-axis projection of the SZ signal.

        Parameters
        ----------
        axis : integer or string
            The axis of the simulation domain along which to make the SZprojection.
        center : A sequence of floats, a string, or a tuple.
            The coordinate of the center of the image. If set to 'c', 'center' or
            left blank, the plot is centered on the middle of the domain. If set to
            'max' or 'm', the center will be located at the maximum of the
            ('gas', 'density') field. Centering on the max or min of a specific
            field is supported by providing a tuple such as ("min","temperature") or
            ("max","dark_matter_density"). Units can be specified by passing in *center*
            as a tuple containing a coordinate and string unit name or by passing
            in a YTArray. If a list or unitless array is supplied, code units are
            assumed.
        width : tuple or a float.
            Width can have four different formats to support windows with variable
            x and y widths.  They are:

            ==================================     =======================
            format                                 example
            ==================================     =======================
            (float, string)                        (10,'kpc')
            ((float, string), (float, string))     ((10,'kpc'),(15,'kpc'))
            float                                  0.2
            (float, float)                         (0.2, 0.3)
            ==================================     =======================

            For example, (10, 'kpc') requests a plot window that is 10 kiloparsecs
            wide in the x and y directions, ((10,'kpc'),(15,'kpc')) requests a
            window that is 10 kiloparsecs wide along the x axis and 15
            kiloparsecs wide along the y axis.  In the other two examples, code
            units are assumed, for example (0.2, 0.3) requests a plot that has an
            x width of 0.2 and a y width of 0.3 in code units.  If units are
            provided the resulting plot axis labels will use the supplied units.
        nx : integer, optional
            The dimensions on a side of the projection image.
        source : yt.data_objects.data_containers.YTSelectionContainer, optional
            If specified, this will be the data source used for selecting regions to project.

        Examples
        --------
        >>> szprj.on_axis("y", center="max", width=(1.0, "Mpc"), source=my_sphere)
        """
        axis = fix_axis(axis, self.ds)
        ctr, dctr = self.ds.coordinates.sanitize_center(center, axis)
        width = self.ds.coordinates.sanitize_width(axis, width, None)

        L = np.zeros(3)
        L[axis] = 1.0

        beta_par = generate_beta_par(L)
        self.ds.add_field(("gas","beta_par"), function=beta_par, units="g/cm**3")
        setup_sunyaev_zeldovich_fields(self.ds)
        proj = self.ds.proj("density", axis, center=ctr, data_source=source)
        frb = proj.to_frb(width[0], nx, height=width[1])
        dens = frb["density"]
        Te = frb["t_sz"]/dens
        bpar = frb["beta_par"]/dens
        omega1 = frb["t_squared"]/dens/(Te*Te) - 1.
        bperp2 = np.zeros((nx,nx))
        sigma1 = np.zeros((nx,nx))
        kappa1 = np.zeros((nx,nx))
        if self.high_order:
            bperp2 = frb["beta_perp_squared"]/dens
            sigma1 = frb["t_beta_par"]/dens/Te - bpar
            kappa1 = frb["beta_par_squared"]/dens - bpar*bpar
        tau = sigma_thompson*dens*self.mueinv/mh

        nx,ny = frb.buff_size
        self.bounds = frb.bounds
        self.dx = (frb.bounds[1]-frb.bounds[0])/nx
        self.dy = (frb.bounds[3]-frb.bounds[2])/ny
        self.nx = nx

        self._compute_intensity(np.array(tau), np.array(Te), np.array(bpar),
                                np.array(omega1), np.array(sigma1),
                                np.array(kappa1), np.array(bperp2))

        self.ds.field_info.pop(("gas","beta_par"))

    def off_axis(self, L, center="c", width=(1.0, "unitary"), depth=(1.0,"unitary"),
                 nx=800, nz=800, north_vector=None, no_ghost=False, source=None):
        r""" Make an off-axis projection of the SZ signal.

        Parameters
        ----------
        L : array_like
            The normal vector of the projection.
        center : A sequence of floats, a string, or a tuple.
            The coordinate of the center of the image. If set to 'c', 'center' or
            left blank, the plot is centered on the middle of the domain. If set to
            'max' or 'm', the center will be located at the maximum of the
            ('gas', 'density') field. Centering on the max or min of a specific
            field is supported by providing a tuple such as ("min","temperature") or
            ("max","dark_matter_density"). Units can be specified by passing in *center*
            as a tuple containing a coordinate and string unit name or by passing
            in a YTArray. If a list or unitless array is supplied, code units are
            assumed.
        width : tuple or a float.
            Width can have four different formats to support windows with variable
            x and y widths.  They are:

            ==================================     =======================
            format                                 example
            ==================================     =======================
            (float, string)                        (10,'kpc')
            ((float, string), (float, string))     ((10,'kpc'),(15,'kpc'))
            float                                  0.2
            (float, float)                         (0.2, 0.3)
            ==================================     =======================

            For example, (10, 'kpc') requests a plot window that is 10 kiloparsecs
            wide in the x and y directions, ((10,'kpc'),(15,'kpc')) requests a
            window that is 10 kiloparsecs wide along the x axis and 15
            kiloparsecs wide along the y axis.  In the other two examples, code
            units are assumed, for example (0.2, 0.3) requests a plot that has an
            x width of 0.2 and a y width of 0.3 in code units.  If units are
            provided the resulting plot axis labels will use the supplied units.
        depth : A tuple or a float
            A tuple containing the depth to project through and the string
            key of the unit: (width, 'unit').  If set to a float, code units
            are assumed
        nx : integer, optional
            The dimensions on a side of the projection image.
        nz : integer, optional
            Deprecated, this is still in the function signature for API
            compatibility
        north_vector : a sequence of floats
            A vector defining the 'up' direction in the plot.  This
            option sets the orientation of the slicing plane.  If not
            set, an arbitrary grid-aligned north-vector is chosen.
        no_ghost: bool, optional
            Optimization option for off-axis cases. If True, homogenized bricks will
            extrapolate out from grid instead of interpolating from
            ghost zones that have to first be calculated.  This can
            lead to large speed improvements, but at a loss of
            accuracy/smoothness in resulting image.  The effects are
            less notable when the transfer function is smooth and
            broad. Default: True
        source : yt.data_objects.data_containers.YTSelectionContainer, optional
            If specified, this will be the data source used for selecting regions 
            to project.

        Examples
        --------
        >>> L = np.array([0.5, 1.0, 0.75])
        >>> szprj.off_axis(L, center="c", width=(2.0, "Mpc"))
        """
        wd = self.ds.coordinates.sanitize_width(L, width, depth)
        w = tuple(el.in_units('code_length').v for el in wd)
        ctr, dctr = self.ds.coordinates.sanitize_center(center, L)
        res = (nx, nx)

        if source is None:
            source = self.ds

        beta_par = generate_beta_par(L)
        self.ds.add_field(("gas","beta_par"), function=beta_par, units="g/cm**3")
        setup_sunyaev_zeldovich_fields(self.ds)

        dens = off_axis_projection(source, ctr, L, w, res, "density",
                                   north_vector=north_vector, no_ghost=no_ghost)
        Te = off_axis_projection(source, ctr, L, w, res, "t_sz",
                                 north_vector=north_vector, no_ghost=no_ghost)/dens
        bpar = off_axis_projection(source, ctr, L, w, res, "beta_par",
                                   north_vector=north_vector, no_ghost=no_ghost)/dens
        omega1 = off_axis_projection(source, ctr, L, w, res, "t_squared",
                                     north_vector=north_vector, no_ghost=no_ghost)/dens
        omega1 = omega1/(Te*Te) - 1.
        if self.high_order:
            bperp2 = off_axis_projection(source, ctr, L, w, res, "beta_perp_squared", 
                                         north_vector=north_vector, no_ghost=no_ghost)/dens
            sigma1 = off_axis_projection(source, ctr, L, w, res, "t_beta_par", 
                                         north_vector=north_vector, no_ghost=no_ghost)/dens
            sigma1 = sigma1/Te - bpar
            kappa1 = off_axis_projection(source, ctr, L, w, res, "beta_par_squared", 
                                         north_vector=north_vector, no_ghost=no_ghost)/dens
            kappa1 -= bpar
        else:
            bperp2 = np.zeros((nx,nx))
            sigma1 = np.zeros((nx,nx))
            kappa1 = np.zeros((nx,nx))
        tau = sigma_thompson*dens*self.mueinv/mh

        self.bounds = (-0.5*wd[0], 0.5*wd[0], -0.5*wd[1], 0.5*wd[1])
        self.dx = wd[0]/nx
        self.dy = wd[1]/nx
        self.nx = nx

        self._compute_intensity(np.array(tau), np.array(Te), np.array(bpar),
                                np.array(omega1), np.array(sigma1),
                                np.array(kappa1), np.array(bperp2))

        self.ds.field_info.pop(("gas","beta_par"))

    def _compute_intensity(self, tau, Te, bpar, omega1, sigma1, kappa1, bperp2):

        # Bad hack, but we get NaNs if we don't do something like this
        small_beta = np.abs(bpar) < 1.0e-20
        bpar[small_beta] = 1.0e-20

        comm = communication_system.communicators[-1]

        nx, ny = self.nx,self.nx
        signal = np.zeros((self.num_freqs,nx,ny))
        xo = np.zeros(self.num_freqs)

        k = int(0)

        start_i = comm.rank*nx//comm.size
        end_i = (comm.rank+1)*nx//comm.size

        pbar = get_pbar("Computing SZ signal.", nx*nx)

        for i in range(start_i, end_i):
            for j in range(ny):
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
            self.data[field] = I0*self.xinit[i]**3*signal[i,:,:]
        self.data["Tau"] = self.ds.arr(tau, "dimensionless")
        self.data["TeSZ"] = self.ds.arr(Te, "keV")

    @parallel_root_only
    def write_fits(self, filename, sky_scale=None, sky_center=None, clobber=True):
        r""" Export images to a FITS file. Writes the SZ distortion in all
        specified frequencies as well as the mass-weighted temperature and the
        optical depth. Distance units are in kpc, unless *sky_center*
        and *scale* are specified. 

        Parameters
        ----------
        filename : string
            The name of the FITS file to be written. 
        sky_scale : tuple
            Conversion between an angle unit and a length unit, if sky
            coordinates are desired, e.g. (1.0, "arcsec/kpc")
        sky_center : tuple, optional
            The (RA, Dec) coordinate in degrees of the central pixel. Must
            be specified with *sky_scale*.
        clobber : boolean, optional
            If the file already exists, do we overwrite?

        Examples
        --------
        >>> # This example just writes out a FITS file with kpc coords
        >>> szprj.write_fits("SZbullet.fits", clobber=False)
        >>> # This example uses sky coords
        >>> sky_scale = (1., "arcsec/kpc") # One arcsec per kpc
        >>> sky_center = (30., 45., "deg")
        >>> szprj.write_fits("SZbullet.fits", sky_center=sky_center, sky_scale=sky_scale)
        """
        from yt.utilities.fits_image import FITSImageData

        dx = self.dx.in_units("kpc")
        dy = dx

        w = _astropy.pywcs.WCS(naxis=2)
        w.wcs.crpix = [0.5*(self.nx+1)]*2
        w.wcs.cdelt = [dx.v,dy.v]
        w.wcs.crval = [0.0,0.0]
        w.wcs.cunit = ["kpc"]*2
        w.wcs.ctype = ["LINEAR"]*2

        fib = FITSImageData(self.data, fields=self.data.keys(), wcs=w)
        if sky_scale is not None and sky_center is not None:
            fib.create_sky_wcs(sky_center, sky_scale)
        fib.writeto(filename, clobber=clobber)

    @parallel_root_only
    def write_png(self, filename_prefix, cmap_name=None,
                  axes_units="kpc", log_fields=None):
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
        if cmap_name is None:
            cmap_name = ytcfg.get("yt", "default_colormap")
        
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        if log_fields is None: log_fields = {}
        ticks_font = matplotlib.font_manager.FontProperties(family='serif',size=16)
        extent = tuple([bound.in_units(axes_units).value for bound in self.bounds])
        for field, image in self.items():
            data = image.copy()
            vmin, vmax = image.min(), image.max()
            negative = False
            crossover = False
            if vmin < 0 and vmax < 0:
                data *= -1
                negative = True
            if field in log_fields:
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
            units = self.data[field].units.latex_representation()
            if units is not None and units != "":
                cbar_label += r'$\ \ ('+units+r')$'
            fig = plt.figure(figsize=(10.0,8.0))
            ax = fig.add_subplot(111)
            cax = ax.imshow(data.d, norm=norm, extent=extent, cmap=cmap_name, origin="lower")
            for label in ax.get_xticklabels():
                label.set_fontproperties(ticks_font)
            for label in ax.get_yticklabels():
                label.set_fontproperties(ticks_font)
            ax.set_xlabel(r"$\mathrm{x\ (%s)}$" % axes_units, fontsize=16)
            ax.set_ylabel(r"$\mathrm{y\ (%s)}$" % axes_units, fontsize=16)
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
        for field, data in self.items():
            data.write_hdf5(filename, dataset_name=field)

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
