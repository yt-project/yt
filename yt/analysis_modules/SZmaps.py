from yt.utilities.physical_constants import sigma_thompson, clight, hcgs, kboltz, mp
from yt.data_objects.image_array import ImageArray
import numpy as np

Tcmb = 2.726
mueinv = 0.875

try:
    import SZpack
except:
    raise ImportError

def _t_squared(field, data):
    return data["TempkeV"]*data["TempkeV"]
add_field("TSquared", function=_t_squared)

def _beta_perp_squared(field, data):
    axis = data.get_field_parameter("axis")
    if axis == "x":
	vv = np.sqrt(data["y-velocity"]**2+data["z-velocity"]**2)
    elif axis == "y":
	vv = np.sqrt(data["x-velocity"]**2+data["z-velocity"]**2)
    elif axis == "z":
	vv = np.sqrt(data["x-velocity"]**2+data["y-velocity"]**2)
    return vv/clight/clight
add_field("BetaPerpSquared", function=_beta_perp_squared)

def _beta_par(field, data):
    axis = data.get_field_parameter("axis")
    return data["%s-velocity" % (axis)]/clight
add_field("BetaPar", function=_beta_par)

def _beta_par_squared(field, data):
    return data["BetaPar"]**2
add_field("BetaParSquared", function=_beta_par_squared)

def _t_beta_par(field, data):
    return data["TempkeV"]*data["BetaPar"]
add_field("TBetaPar", function=_t_beta_par)

vlist = 'xyz'

def SZProjection(pf, axis, freqs, width=(1, "unitary"), nx=800, ny=800)

    num_freqs = len(freqs)
    freq_fields = ["%d_GHz" % (int(freq)) for freq in freqs]
    xo = hcgs*freqs*1.0e9/(kboltz*Tcmb)

    proj1 = pf.h.proj("TempkeV", weight_field="Density")
    proj2 = pf.h.proj("Density")

    if axis in vlist:
	vfield = "velocity_%s" % (axis)
	proj1.set_field_parameter("axis", axis)
    elif axis in xrange(0,3) :
	vfield = "velocity_%s" % (vlist[axis])
	proj1.set_field_parameter("axis", vlist[axis])
    
    frb1 = proj1.to_frb(width, n)
    frb2 = proj2.to_frb(width, n)
    
    TeSZ = frb1["TempkeV"]
    omega1 = frb1["Tsquared"]/(TeSZ*TeSZ) - 1.
    sigma1 = frb1["TBetaPar"]/TeSZ - betac_par
    kappa1 = frb1["BetaParSquared"] - betac_par
    
    frb1["tau"] = sigma_thompson*frb2["Density"]*mueinv/mp
    frb1["omega1"] = ImageArray(omega1)
    frb1["kappa1"] = ImageArray(kappa1)
    frb1["sigma1"] = ImageArray(sigma1)

    SZsignal = np.zeros((num_freqs,nx,ny))
    omega = np.zeros((3))
    sigma = np.zeros((3))
    
    for i in xrange(nx):

	for j in xrange(ny):
		
	    tau = frb1["tau"][i,j]
	    Te = frb1["TempkeV"][i,j]
	    bpar = frb1["BetaPar"][i,j]
	    bperp2 = frb["BetaPerpSquared"][i,j]
	    omega[0] = frb1["omega1"][i,j]
	    sigma[0] = frb1["sigma1"][i,j]
	    kappa = frb1["kappa1"][i,j]
	
	    SZsignal[:,i,j] = SZpack.compute_combo_means_ex(xo, tau, Te, bpar, omega,
							    sigma, kappa, bperp2)

    for i in xrange(num_freqs) :
	frb1[freq_fields[i]] = ImageArray(SZsignal[i,:,:])
	
    return frb1
