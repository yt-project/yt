from yt.utilities.physical_constants import sigma_thompson, clight, hcgs, kboltz, mp
from yt.data_objects.image_array import ImageArray
from yt.data_objects.field_info_container import add_field
from yt.funcs import fix_axis, mylog, iterable, get_pbar
from yt.definitions import inv_axis_names
from yt.visualization.volume_rendering.camera import off_axis_projection
import numpy as np

Tcmb = 2.726

try:
    import SZpack
except:
    raise ImportError("SZpack not installed.")

vlist = "xyz"
nvec = int(0)

def _t_squared(field, data):
    return data["TempkeV"]*data["TempkeV"]
add_field("TSquared", function=_t_squared)

def _beta_perp_squared(field, data):
    return data["VelocityMagnitude"]**2/clight/clight - data["BetaParSquared"]
add_field("BetaPerpSquared", function=_beta_perp_squared)

def _beta_par(field, data):
    axis = data.get_field_parameter("axis")
    if iterable(nvec):
        vpar = (data["x-velocity"]*nvec[0]+
                data["y-velocity"]*nvec[1]+
                data["z-velocity"]*nvec[2])
    else:
        vpar = data["%s-velocity" % (vlist[nvec])]
    return vpar/clight
add_field("BetaPar", function=_beta_par)

def _beta_par_squared(field, data):
    return data["BetaPar"]**2
add_field("BetaParSquared", function=_beta_par_squared)

def _t_beta_par(field, data):
    return data["TempkeV"]*data["BetaPar"]
add_field("TBetaPar", function=_t_beta_par)

def SZProjection(pf, axis, freqs, center="c", width=(1, "unitary"), nx=800, mue=None):

    global nvec # Ugly!
    
    if mue is None:
        mueinv = 0.875
    else:
        mueinv = 1./mue
        
    num_freqs = len(freqs)
    freq_fields = ["%d_GHz" % (int(freq)) for freq in freqs]
    xinit = hcgs*freqs*1.0e9/(kboltz*Tcmb)

    if isinstance(axis, basestring) or not iterable(axis):
        axis = fix_axis(axis)
        nvec = axis
        proj1 = pf.h.proj(axis, "TempkeV", weight_field="Density")
        proj2 = pf.h.proj(axis, "Density")
        frb1 = proj1.to_frb(width, nx)
        frb2 = proj2.to_frb(width, nx)
        Te = frb1["TempkeV"]
        bpar = frb1["BetaPar"]
        bperp2 = frb1["BetaPerpSquared"]
        omega1 = frb1["TSquared"]/(Te*Te) - 1.
        sigma1 = frb1["TBetaPar"]/Te - bpar
        kappa1 = frb1["BetaParSquared"] - bpar
        tau = sigma_thompson*frb2["Density"]*mueinv/mp
    else:
        nvec = axis
        if iterable(width):
            w = width[0]/pf.units[width[1]]
        else:
            w = width
        Te      = off_axis_projection(pf, center, axis, w, nx, "TempkeV", weight="Density")
        bpar    = off_axis_projection(pf, center, axis, w, nx, "BetaPar", weight="Density")
        bperp2  = off_axis_projection(pf, center, axis, w, nx, "BetaPerpSquared", weight="Density")
        omega1  = off_axis_projection(pf, center, axis, w, nx, "TSquared", weight="Density")
        omega1  = omega1/(Te*Te) - 1.
        sigma1  = off_axis_projection(pf, center, axis, w, nx, "TBetaPar", weight="Density")
        sigma1  = sigma1/Te - bpar
        kappa1  = off_axis_projection(pf, center, axis, w, nx, "BetaParSquared", weight="Density")
        kappa1 -= bpar
        tau     = off_axis_projection(pf, center, axis, w, nx, "Density")
        tau    *= sigma_thompson*mueinv/mp
        
    SZsignal = np.zeros((num_freqs,nx,nx))
    xo = np.zeros((num_freqs))
    
    k = int(0)

    pbar = get_pbar("Computing SZ signal.", nx*nx)
    
    for i in xrange(nx):
	for j in xrange(nx):
            xo[:] = xinit[:]
	    SZpack.compute_combo_means(xo, tau[i,j], Te[i,j],
                                       bpar[i,j], omega[i,j],
                                       sigma[i,j], kappa[i,j], bperp2[i,j])
            SZsignal[:,i,j] = -xo[:]
            pbar.update(k)
            k += 1

    pbar.finish()
    
    outimg = {}
    for i in xrange(num_freqs) :
	outimg[freq_fields[i]] = ImageArray(SZsignal[i,:,:])
	
    return outimg
