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

from yt.utilities.physical_constants import sigma_thompson, clight, hcgs, kboltz, mp
from yt.data_objects.image_array import ImageArray
from yt.data_objects.field_info_container import add_field
from yt.funcs import fix_axis, mylog, iterable, get_pbar
from yt.utilities.definitions import inv_axis_names
from yt.visualization.image_writer import write_fits, write_projection
from yt.visualization.volume_rendering.camera import off_axis_projection
from yt.utilities.parallel_tools.parallel_analysis_interface import \
     communication_system, parallel_root_only
import numpy as np

Tcmb = 2.726
I0 = 2*(kboltz*Tcmb)**3/((hcgs*clight)**2)*1.0e17
        
try:
    import SZpack
except:
    raise ImportError("SZpack not installed. It can be obtained from from http://www.chluba.de/SZpack/.")

vlist = "xyz"

def _t_squared(field, data):
    return data["Density"]*data["TempkeV"]*data["TempkeV"]
add_field("TSquared", function=_t_squared)

def _beta_perp_squared(field, data):
    return data["Density"]*(data["VelocityMagnitude"]**2/clight/clight - data["BetaParSquared"])
add_field("BetaPerpSquared", function=_beta_perp_squared)

def _beta_par_squared(field, data):
    return data["Density"]*data["BetaPar"]**2
add_field("BetaParSquared", function=_beta_par_squared)

def _t_beta_par(field, data):
    return data["Density"]*data["TempkeV"]*data["BetaPar"]
add_field("TBetaPar", function=_t_beta_par)

def _t_sz(field, data):
    return data["Density"]*data["TempkeV"]
add_field("TeSZ", function=_t_sz)

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
        self.xinit = hcgs*freqs*1.0e9/(kboltz*Tcmb)
        self.freq_fields = ["%d_GHz" % (int(freq)) for freq in freqs]
        self.field_dict = {}

        self.units = {}
        self.units["TeSZ"] = r"$\mathrm{keV}$"
        self.units["Tau"] = None

        self.display_names = {}
        self.display_names["TeSZ"] = r"$\mathrm{T_e}$"
        self.display_names["Tau"] = r"$\mathrm{\tau}$"

        for f, field in zip(self.freqs, self.freq_fields):
            self.units[field] = r"$\mathrm{MJy\ sr^{-1}}$"
            self.display_names[field] = r"$\mathrm{\Delta{I}_{%d\ GHz}}$" % (int(freq))
            
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

        def _beta_par(field, data):
            axis = data.get_field_parameter("axis")
            vpar = data["Density"]*data["%s-velocity" % (vlist[axis])]
            return vpar/clight
        add_field("BetaPar", function=_beta_par)    

        proj = self.pf.h.proj(axis, "Density", source=source)
        proj.set_field_parameter("axis", axis)
        frb = proj.to_frb(width, nx)
        dens = frb["Density"]
        Te = frb["TeSZ"]/dens
        bpar = frb["BetaPar"]/dens
        omega1 = frb["TSquared"]/dens/(Te*Te) - 1.
        if self.high_order:
            bperp2 = frb["BetaPerpSquared"]/dens
            sigma1 = frb["TBetaPar"]/dens/Te - bpar
            kappa1 = frb["BetaParSquared"]/dens - bpar
        else:
            bperp2 = np.zeros((nx,nx))
            sigma1 = np.zeros((nx,nx))
            kappa1 = np.zeros((nx,nx))
        tau = sigma_thompson*dens*self.mueinv/mp

        nx,ny = frb.buff_size
        self.bounds = frb.bounds
        self.dx = (frb.bounds[1]-frb.bounds[0])/nx
        self.dy = (frb.bounds[3]-frb.bounds[2])/ny
        
        self._compute_intensity(tau, Te, bpar, omega1, sigma1, kappa1, bperp2)
                                                                                                                
    def off_axis(self, L, center="c", width=(1, "unitary"), nx=800):
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
            ctr = self.pf.h.find_max("Density")
        else:
            ctr = center
            
        def _beta_par(field, data):
            vpar = data["Density"]*(data["x-velocity"]*L[0]+
                                    data["y-velocity"]*L[1]+
                                    data["z-velocity"]*L[2])
            return vpar/clight
        add_field("BetaPar", function=_beta_par)

        dens    = off_axis_projection(self.pf, ctr, L, w, nx, "Density")
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
        tau = sigma_thompson*dens*self.mueinv/mp

        self.bounds = np.array([-0.5*w, 0.5*w, -0.5*w, 0.5*w])
        self.dx = w/nx
        self.dy = w/nx
        
        self._compute_intensity(tau, Te, bpar, omega1, sigma1, kappa1, bperp2)

    def _compute_intensity(self, tau, Te, bpar, omega1, sigma1, kappa1, bperp2):

        comm = communication_system.communicators[-1]
        
        nx, ny = tau.shape
        signal = np.zeros((self.num_freqs,nx,nx))
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
                signal[:,i,j] = -xo[:]
                pbar.update(k)
                k += 1

        signal = comm.mpi_allreduce(signal)
        
        pbar.finish()
                
        for i, field in enumerate(self.freq_fields):
            self.field_dict[field] = ImageArray(I0*self.xinit[i]**3*signal[i,:,:])
        self.field_dict["Tau"] = ImageArray(tau)
        self.field_dict["TeSZ"] = ImageArray(Te)

    @parallel_root_only
    def write_fits(self, filename_prefix, clobber=True):
        r""" Export images to a FITS file. Writes the SZ distortion in all
        specified frequencies as well as the mass-weighted temperature and the
        optical depth. Distance units are in kpc.  
        
        Parameters
        ----------
        filename_prefix : string
            The prefix of the FITS filename.
        clobber : boolean, optional
            If the file already exists, do we overwrite?
                    
        Examples
        --------
        >>> szprj.write_fits("SZbullet", clobber=False)
        """
        coords = {}
        coords["dx"] = self.dx*self.pf.units["kpc"]
        coords["dy"] = self.dy*self.pf.units["kpc"]
        coords["xctr"] = 0.0
        coords["yctr"] = 0.0
        coords["units"] = "kpc"
        other_keys = {"Time" : self.pf.current_time}
        write_fits(self.field_dict, filename_prefix, clobber=clobber, coords=coords,
                   other_keys=other_keys)

    @parallel_root_only
    def write_png(self, filename_prefix):
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
        extent = tuple([bound*self.pf.units["kpc"] for bound in self.bounds])
        for field, image in self.field_dict.items():
            filename=filename_prefix+"_"+field+".png"
            label = self.display_names[field]
            if self.units[field] is not None:
                label += " ("+self.units[field]+")"
            write_projection(image, filename, colorbar_label=label, take_log=False,
                             extent=extent, xlabel=r"$\mathrm{x\ (kpc)}$",
                             ylabel=r"$\mathrm{y\ (kpc)}$")

    def keys(self):
        return self.field_dict.keys()

    def has_key(self, key):
        return key in self.field_dict.keys()

    def __getitem__(self, key):
        return self.field_dict[key]
