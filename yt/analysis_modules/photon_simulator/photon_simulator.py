"""
Classes for generating lists of photons and detected events
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
from yt.utilities.physical_constants import clight
from yt.utilities.cosmology import Cosmology
from yt.utilities.orientation import Orientation
from yt.utilities.parallel_tools.parallel_analysis_interface import \
     communication_system, parallel_root_only, get_mpi_type, \
     parallel_capable
from yt.units.yt_array import YTQuantity, YTArray, uconcatenate
import h5py
from yt.utilities.on_demand_imports import _astropy

comm = communication_system.communicators[-1]

def parse_value(value, default_units):
    if isinstance(value, YTQuantity):
        return value.in_units(default_units)
    elif iterable(value):
        return YTQuantity(value[0], value[1]).in_units(default_units)
    else:
        return YTQuantity(value, default_units)

class PhotonList(object):

    def __init__(self, photons, parameters, cosmo, p_bins):
        self.photons = photons
        self.parameters = parameters
        self.cosmo = cosmo
        self.p_bins = p_bins
        self.num_cells = len(photons["x"])
        
    def keys(self):
        return self.photons.keys()
    
    def items(self):
        ret = []
        for k, v in self.photons.items():
            if k == "Energy":
                ret.append((k, self[k]))
            else:
                ret.append((k,v))
        return ret
    
    def values(self):
        ret = []
        for k, v in self.photons.items():
            if k == "Energy":
                ret.append(self[k])
            else:
                ret.append(v)
        return ret
                                
    def __getitem__(self, key):
        if key == "Energy":
            return [self.photons["Energy"][self.p_bins[i]:self.p_bins[i+1]]
                    for i in xrange(self.num_cells)]
        else:
            return self.photons[key]

    def __repr__(self):
        return self.photons.__repr__()

    @classmethod
    def from_file(cls, filename):
        r"""
        Initialize a PhotonList from the HDF5 file *filename*.
        """

        photons = {}
        parameters = {}
        
        f = h5py.File(filename, "r")

        parameters["FiducialExposureTime"] = YTQuantity(f["/fid_exp_time"].value, "s")
        parameters["FiducialArea"] = YTQuantity(f["/fid_area"].value, "cm**2")
        parameters["FiducialRedshift"] = f["/fid_redshift"].value
        parameters["FiducialAngularDiameterDistance"] = YTQuantity(f["/fid_d_a"].value, "Mpc")
        parameters["Dimension"] = f["/dimension"].value
        parameters["Width"] = YTQuantity(f["/width"].value, "kpc")
        parameters["HubbleConstant"] = f["/hubble"].value
        parameters["OmegaMatter"] = f["/omega_matter"].value
        parameters["OmegaLambda"] = f["/omega_lambda"].value

        num_cells = f["/x"][:].shape[0]
        start_c = comm.rank*num_cells/comm.size
        end_c = (comm.rank+1)*num_cells/comm.size
        
        photons["x"] = YTArray(f["/x"][start_c:end_c], "kpc")
        photons["y"] = YTArray(f["/y"][start_c:end_c], "kpc")
        photons["z"] = YTArray(f["/z"][start_c:end_c], "kpc")
        photons["dx"] = YTArray(f["/dx"][start_c:end_c], "kpc")
        photons["vx"] = YTArray(f["/vx"][start_c:end_c], "km/s")
        photons["vy"] = YTArray(f["/vy"][start_c:end_c], "km/s")
        photons["vz"] = YTArray(f["/vz"][start_c:end_c], "km/s")

        n_ph = f["/num_photons"][:]
        
        if comm.rank == 0:
            start_e = np.uint64(0)
        else:
            start_e = n_ph[:start_c].sum()
        end_e = start_e + np.uint64(n_ph[start_c:end_c].sum())

        photons["NumberOfPhotons"] = n_ph[start_c:end_c]

        p_bins = np.cumsum(photons["NumberOfPhotons"])
        p_bins = np.insert(p_bins, 0, [np.uint64(0)])
        
        photons["Energy"] = YTArray(f["/energy"][start_e:end_e], "keV")
        
        f.close()

        cosmo = Cosmology(hubble_constant=parameters["HubbleConstant"],
                          omega_matter=parameters["OmegaMatter"],
                          omega_lambda=parameters["OmegaLambda"])

        return cls(photons, parameters, cosmo, p_bins)
    
    @classmethod
    def from_scratch(cls, data_source, redshift, area,
                     exp_time, photon_model, parameters=None,
                     center=None, dist=None, cosmology=None):
        r"""
        Initialize a PhotonList from a photon model. The redshift, collecting area,
        exposure time, and cosmology are stored in the *parameters* dictionary which
        is passed to the *photon_model* function. 

        Parameters
        ----------

        data_source : `yt.data_objects.data_containers.YTSelectionContainer`
            The data source from which the photons will be generated.
        redshift : float
            The cosmological redshift for the photons.
        area : float
            The collecting area to determine the number of photons in cm^2.
        exp_time : float
            The exposure time to determine the number of photons in seconds.
        photon_model : function
            A function that takes the *data_source* and the *parameters*
            dictionary and returns a *photons* dictionary. Must be of the
            form: photon_model(data_source, parameters)
        parameters : dict, optional
            A dictionary of parameters to be passed to the user function. 
        center : string or array_like, optional
            The origin of the photons. Accepts "c", "max", or a coordinate.                
        dist : tuple, optional
            The angular diameter distance in the form (value, unit), used
            mainly for nearby sources. This may be optionally supplied
            instead of it being determined from the *redshift* and given *cosmology*.
        cosmology : `yt.utilities.cosmology.Cosmology`, optional
            Cosmological information. If not supplied, we try to get
            the cosmology from the dataset. Otherwise, \LambdaCDM with
            the default yt parameters is assumed.

        Examples
        --------

        This is the simplest possible example, where we call the built-in thermal model:

        >>> thermal_model = ThermalPhotonModel(apec_model, Zmet=0.3)
        >>> redshift = 0.05
        >>> area = 6000.0
        >>> time = 2.0e5
        >>> sp = ds.sphere("c", (500., "kpc"))
        >>> my_photons = PhotonList.from_user_model(sp, redshift, area,
        ...                                         time, thermal_model)

        If you wish to make your own photon model function, it must take as its
        arguments the *data_source* and the *parameters* dictionary. However you
        determine them, the *photons* dict needs to have the following items, corresponding
        to cells which have photons:

        "x" : the x-position of the cell relative to the source center in kpc, YTArray
        "y" : the y-position of the cell relative to the source center in kpc, YTArray
        "z" : the z-position of the cell relative to the source center in kpc, YTArray
        "vx" : the x-velocity of the cell in km/s, YTArray
        "vy" : the y-velocity of the cell in km/s, YTArray
        "vz" : the z-velocity of the cell in km/s, YTArray
        "dx" : the width of the cell in kpc, YTArray
        "NumberOfPhotons" : the number of photons in the cell, NumPy array of unsigned 64-bit integers
        "Energy" : the source rest-frame energies of the photons, YTArray

        The last array is not the same size as the others because it contains the energies in all of
        the cells in a single 1-D array. The first photons["NumberOfPhotons"][0] elements are
        for the first cell, the next photons["NumberOfPhotons"][1] are for the second cell, and so on.

        The following is a simple example where a point source with a single line emission
        spectrum of photons is created. More complicated examples which actually
        create photons based on the fields in the dataset could be created. 

        >>> import numpy as np
        >>> import yt
        >>> from yt.analysis_modules.photon_simulator import *
        >>> def line_func(source, parameters):
        ...
        ...     ds = source.ds
        ... 
        ...     num_photons = parameters["num_photons"]
        ...     E0  = parameters["line_energy"] # Energies are in keV
        ...     sigE = parameters["line_sigma"] 
        ...     src_ctr = parameters["center"]
        ...
        ...     energies = norm.rvs(loc=E0, scale=sigE, size=num_photons)
        ...
        ...     # Place everything in the center cell
        ...     for i, ax in enumerate("xyz"):
        ...         photons[ax] = (ds.domain_center[0]-src_ctr[0]).in_units("kpc")
        ...     photons["vx"] = ds.arr([0], "km/s")
        ...     photons["vy"] = ds.arr([0], "km/s")
        ...     photons["vz"] = ds.arr([100.0], "km/s")
        ...     photons["dx"] = ds.find_field_values_at_point("dx", ds.domain_center).in_units("kpc")
        ...     photons["NumberOfPhotons"] = np.array(num_photons*np.ones(1), dtype="uint64")
        ...     photons["Energy"] = ds.arr(energies, "keV")
        >>>
        >>> redshift = 0.05
        >>> area = 6000.0
        >>> time = 2.0e5
        >>> parameters = {"num_photons" : 10000, "line_energy" : 5.0,
        ...               "line_sigma" : 0.1}
        >>> ddims = (128,128,128)
        >>> random_data = {"density":(np.random.random(ddims),"g/cm**3")}
        >>> ds = yt.load_uniform_grid(random_data, ddims)
        >>> dd = ds.all_data
        >>> my_photons = PhotonList.from_user_model(dd, redshift, area,
        ...                                         time, line_func,
        ...                                         parameters=parameters)

        """

        ds = data_source.ds

        if parameters is None:
             parameters = {}
        if cosmology is None:
            hubble = getattr(ds, "hubble_constant", None)
            omega_m = getattr(ds, "omega_matter", None)
            omega_l = getattr(ds, "omega_lambda", None)
            if hubble == 0: hubble = None
            if hubble is not None and \
               omega_m is not None and \
               omega_l is not None:
                cosmo = Cosmology(hubble_constant=hubble,
                                  omega_matter=omega_m,
                                  omega_lambda=omega_l)
            else:
                cosmo = Cosmology()
        else:
            cosmo = cosmology
        mylog.info("Cosmology: h = %g, omega_matter = %g, omega_lambda = %g" %
                   (cosmo.hubble_constant, cosmo.omega_matter, cosmo.omega_lambda))
        if dist is None:
            D_A = cosmo.angular_diameter_distance(0.0,redshift).in_units("Mpc")
        else:
            D_A = parse_value(dist, "Mpc")
            redshift = 0.0

        if center == "c":
            parameters["center"] = ds.domain_center
        elif center == "max":
            parameters["center"] = ds.find_max("density")[-1]
        elif iterable(center):
            if isinstance(center, YTArray):
                parameters["center"] = center.in_units("code_length")
            else:
                parameters["center"] = ds.arr(center, "code_length")
        elif center is None:
            parameters["center"] = data_source.get_field_parameter("center")
            
        parameters["FiducialExposureTime"] = parse_value(exp_time, "s")
        parameters["FiducialArea"] = parse_value(area, "cm**2")
        parameters["FiducialRedshift"] = redshift
        parameters["FiducialAngularDiameterDistance"] = D_A
        parameters["HubbleConstant"] = cosmo.hubble_constant
        parameters["OmegaMatter"] = cosmo.omega_matter
        parameters["OmegaLambda"] = cosmo.omega_lambda

        dimension = 0
        width = 0.0
        for i, ax in enumerate("xyz"):
            le, re = data_source.quantities.extrema(ax)
            delta_min, delta_max = data_source.quantities.extrema("d%s"%ax)
            le -= 0.5*delta_max
            re += 0.5*delta_max
            width = max(width, re-parameters["center"][i], parameters["center"][i]-le)
            dimension = max(dimension, int(width/delta_min))
        parameters["Dimension"] = 2*dimension
        parameters["Width"] = 2.*width.in_units("kpc")
                
        photons = photon_model(data_source, parameters)

        mylog.info("Finished generating photons.")

        p_bins = np.cumsum(photons["NumberOfPhotons"])
        p_bins = np.insert(p_bins, 0, [np.uint64(0)])
                        
        return cls(photons, parameters, cosmo, p_bins)
        
    def write_h5_file(self, photonfile):
        """
        Write the photons to the HDF5 file *photonfile*.
        """

        if parallel_capable:
            
            mpi_long = get_mpi_type("int64")
            mpi_double = get_mpi_type("float64")
        
            local_num_cells = len(self.photons["x"])
            sizes_c = comm.comm.gather(local_num_cells, root=0)
            
            local_num_photons = self.photons["NumberOfPhotons"].sum()
            sizes_p = comm.comm.gather(local_num_photons, root=0)
            
            if comm.rank == 0:
                num_cells = sum(sizes_c)
                num_photons = sum(sizes_p)        
                disps_c = [sum(sizes_c[:i]) for i in range(len(sizes_c))]
                disps_p = [sum(sizes_p[:i]) for i in range(len(sizes_p))]
                x = np.zeros(num_cells)
                y = np.zeros(num_cells)
                z = np.zeros(num_cells)
                vx = np.zeros(num_cells)
                vy = np.zeros(num_cells)
                vz = np.zeros(num_cells)
                dx = np.zeros(num_cells)
                n_ph = np.zeros(num_cells, dtype="uint64")
                e = np.zeros(num_photons)
            else:
                sizes_c = []
                sizes_p = []
                disps_c = []
                disps_p = []
                x = np.empty([])
                y = np.empty([])
                z = np.empty([])
                vx = np.empty([])
                vy = np.empty([])
                vz = np.empty([])
                dx = np.empty([])
                n_ph = np.empty([])
                e = np.empty([])
                                                
            comm.comm.Gatherv([self.photons["x"].ndarray_view(), local_num_cells, mpi_double],
                              [x, (sizes_c, disps_c), mpi_double], root=0)
            comm.comm.Gatherv([self.photons["y"].ndarray_view(), local_num_cells, mpi_double],
                              [y, (sizes_c, disps_c), mpi_double], root=0)
            comm.comm.Gatherv([self.photons["z"].ndarray_view(), local_num_cells, mpi_double],
                              [z, (sizes_c, disps_c), mpi_double], root=0)
            comm.comm.Gatherv([self.photons["vx"].ndarray_view(), local_num_cells, mpi_double],
                              [vx, (sizes_c, disps_c), mpi_double], root=0)
            comm.comm.Gatherv([self.photons["vy"].ndarray_view(), local_num_cells, mpi_double],
                              [vy, (sizes_c, disps_c), mpi_double], root=0)
            comm.comm.Gatherv([self.photons["vz"].ndarray_view(), local_num_cells, mpi_double],
                              [vz, (sizes_c, disps_c), mpi_double], root=0)
            comm.comm.Gatherv([self.photons["dx"].ndarray_view(), local_num_cells, mpi_double],
                              [dx, (sizes_c, disps_c), mpi_double], root=0)
            comm.comm.Gatherv([self.photons["NumberOfPhotons"], local_num_cells, mpi_long],
                              [n_ph, (sizes_c, disps_c), mpi_long], root=0)
            comm.comm.Gatherv([self.photons["Energy"].ndarray_view(), local_num_photons, mpi_double],
                              [e, (sizes_p, disps_p), mpi_double], root=0) 

        else:

            x = self.photons["x"].d
            y = self.photons["y"].d
            z = self.photons["z"].d
            vx = self.photons["vx"].d
            vy = self.photons["vy"].d
            vz = self.photons["vz"].d
            dx = self.photons["dx"].d
            n_ph = self.photons["NumberOfPhotons"]
            e = self.photons["Energy"].d
                                                
        if comm.rank == 0:
            
            f = h5py.File(photonfile, "w")

            # Scalars
       
            f.create_dataset("fid_area", data=float(self.parameters["FiducialArea"]))
            f.create_dataset("fid_exp_time", data=float(self.parameters["FiducialExposureTime"]))
            f.create_dataset("fid_redshift", data=self.parameters["FiducialRedshift"])
            f.create_dataset("hubble", data=self.parameters["HubbleConstant"])
            f.create_dataset("omega_matter", data=self.parameters["OmegaMatter"])
            f.create_dataset("omega_lambda", data=self.parameters["OmegaLambda"])
            f.create_dataset("fid_d_a", data=float(self.parameters["FiducialAngularDiameterDistance"]))
            f.create_dataset("dimension", data=self.parameters["Dimension"])
            f.create_dataset("width", data=float(self.parameters["Width"]))
                        
            # Arrays

            f.create_dataset("x", data=x)
            f.create_dataset("y", data=y)
            f.create_dataset("z", data=z)
            f.create_dataset("vx", data=vx)
            f.create_dataset("vy", data=vy)
            f.create_dataset("vz", data=vz)
            f.create_dataset("dx", data=dx)
            f.create_dataset("num_photons", data=n_ph)
            f.create_dataset("energy", data=e)

            f.close()

        comm.barrier()

    def project_photons(self, L, area_new=None, exp_time_new=None, 
                        redshift_new=None, dist_new=None,
                        absorb_model=None, psf_sigma=None,
                        sky_center=None, responses=None,
                        convolve_energies=False, no_shifting=False):
        r"""
        Projects photons onto an image plane given a line of sight.

        Parameters
        ----------

        L : array_like
            Normal vector to the plane of projection.
        area_new : float, optional
            New value for the effective area of the detector. If *responses*
            are specified the value of this keyword is ignored.
        exp_time_new : float, optional
            The new value for the exposure time.
        redshift_new : float, optional
            The new value for the cosmological redshift.
        dist_new : tuple, optional
            The new value for the angular diameter distance in the form
            (value, unit), used mainly for nearby sources. This may be optionally supplied
            instead of it being determined from the cosmology.
        absorb_model : 'yt.analysis_modules.photon_simulator.PhotonModel`, optional
            A model for galactic absorption.
        psf_sigma : float, optional
            Quick-and-dirty psf simulation using Gaussian smoothing with
            standard deviation *psf_sigma* in degrees. 
        sky_center : array_like, optional
            Center RA, Dec of the events in degrees.
        responses : list of strings, optional
            The names of the ARF and/or RMF files to convolve the photons with.
        convolve_energies : boolean, optional
            If this is set, the photon energies will be convolved with the RMF.
        no_shifting : boolean, optional
            If set, the photon energies will not be Doppler shifted.
            
        Examples
        --------
        >>> L = np.array([0.1,-0.2,0.3])
        >>> events = my_photons.project_photons(L, area_new="sim_arf.fits",
        ...                                     redshift_new=0.05,
        ...                                     psf_sigma=0.01)
        """

        if redshift_new is not None and dist_new is not None:
            mylog.error("You may specify a new redshift or distance, "+
                        "but not both!")
        
        if sky_center is None:
            sky_center = YTArray([30.,45.], "degree")
        else:
            sky_center = YTArray(sky_center, "degree")

        dx = self.photons["dx"].d
        nx = self.parameters["Dimension"]
        if psf_sigma is not None:
             psf_sigma = parse_value(psf_sigma, "degree")

        L = np.array(L)
        L /= np.sqrt(np.dot(L, L))
        vecs = np.identity(3)
        t = np.cross(L, vecs).sum(axis=1)
        ax = t.argmax()
        north = np.cross(L, vecs[ax,:]).ravel()
        orient = Orientation(L, north_vector=north)
        
        x_hat = orient.unit_vectors[0]
        y_hat = orient.unit_vectors[1]
        z_hat = orient.unit_vectors[2]

        n_ph = self.photons["NumberOfPhotons"]
        n_ph_tot = n_ph.sum()
        
        eff_area = None

        parameters = {}
        
        if responses is not None:
            responses = ensure_list(responses)
            parameters["ARF"] = responses[0]
            if len(responses) == 2:
                parameters["RMF"] = responses[1]
            area_new = parameters["ARF"]
            
        if (exp_time_new is None and area_new is None and
            redshift_new is None and dist_new is None):
            my_n_obs = n_ph_tot
            zobs = self.parameters["FiducialRedshift"]
            D_A = self.parameters["FiducialAngularDiameterDistance"]
        else:
            if exp_time_new is None:
                Tratio = 1.
            else:
                Tratio = parse_value(exp_time_new, "s")/self.parameters["FiducialExposureTime"]
            if area_new is None:
                Aratio = 1.
            elif isinstance(area_new, basestring):
                if comm.rank == 0:
                    mylog.info("Using energy-dependent effective area: %s" % (parameters["ARF"]))
                f = _astropy.pyfits.open(area_new)
                elo = f["SPECRESP"].data.field("ENERG_LO")
                ehi = f["SPECRESP"].data.field("ENERG_HI")
                eff_area = np.nan_to_num(f["SPECRESP"].data.field("SPECRESP"))
                if "RMF" in parameters:
                    weights = self._normalize_arf(parameters["RMF"])
                    eff_area *= weights
                else:
                    mylog.warning("You specified an ARF but not an RMF. This is ok if the "+
                                  "responses are normalized properly. If not, you may "+
                                  "get inconsistent results.")
                f.close()
                Aratio = eff_area.max()/self.parameters["FiducialArea"]
            else:
                mylog.info("Using constant effective area.")
                Aratio = parse_value(area_new, "cm**2")/self.parameters["FiducialArea"]
            if redshift_new is None and dist_new is None:
                Dratio = 1.
                zobs = self.parameters["FiducialRedshift"]
                D_A = self.parameters["FiducialAngularDiameterDistance"]
            else:
                if redshift_new is None:
                    zobs = 0.0
                    D_A = parse_value(dist_new, "Mpc")
                else:
                    zobs = redshift_new
                    D_A = self.cosmo.angular_diameter_distance(0.0,zobs).in_units("Mpc")
                fid_D_A = self.parameters["FiducialAngularDiameterDistance"]
                Dratio = fid_D_A*fid_D_A*(1.+self.parameters["FiducialRedshift"]**3) / \
                         (D_A*D_A*(1.+zobs)**3)
            fak = Aratio*Tratio*Dratio
            if fak > 1:
                raise ValueError("Spectrum scaling factor = %g, cannot be greater than unity." % (fak))
            my_n_obs = np.uint64(n_ph_tot*fak)

        n_obs_all = comm.mpi_allreduce(my_n_obs)
        if comm.rank == 0:
            mylog.info("Total number of photons to use: %d" % (n_obs_all))
        
        x = np.random.uniform(low=-0.5,high=0.5,size=my_n_obs)
        y = np.random.uniform(low=-0.5,high=0.5,size=my_n_obs)
        z = np.random.uniform(low=-0.5,high=0.5,size=my_n_obs)

        if not no_shifting:
            vz = self.photons["vx"]*z_hat[0] + \
                 self.photons["vy"]*z_hat[1] + \
                 self.photons["vz"]*z_hat[2]
            shift = -vz.in_cgs()/clight
            shift = np.sqrt((1.-shift)/(1.+shift))

        if my_n_obs == n_ph_tot:
            idxs = np.arange(my_n_obs,dtype='uint64')
        else:
            idxs = np.random.permutation(n_ph_tot)[:my_n_obs].astype("uint64")
        obs_cells = np.searchsorted(self.p_bins, idxs, side='right')-1
        delta = dx[obs_cells]

        x *= delta
        y *= delta
        z *= delta
        x += self.photons["x"][obs_cells]
        y += self.photons["y"][obs_cells]
        z += self.photons["z"][obs_cells]
        if no_shifting:
            eobs = self.photons["Energy"][idxs]
        else:
            eobs = self.photons["Energy"][idxs]*shift[obs_cells]

        xsky = x*x_hat[0] + y*x_hat[1] + z*x_hat[2]
        ysky = x*y_hat[0] + y*y_hat[1] + z*y_hat[2]
        eobs /= (1.+zobs)
        
        if absorb_model is None:
            not_abs = np.ones(eobs.shape, dtype='bool')
        else:
            mylog.info("Absorbing.")
            absorb_model.prepare()
            emid = absorb_model.emid
            aspec = absorb_model.get_spectrum()
            absorb = np.interp(eobs, emid, aspec, left=0.0, right=0.0)
            randvec = aspec.max()*np.random.random(eobs.shape)
            not_abs = randvec < absorb
        
        if eff_area is None:
            detected = np.ones(eobs.shape, dtype='bool')
        else:
            mylog.info("Applying energy-dependent effective area.")
            earf = 0.5*(elo+ehi)
            earea = np.interp(eobs, earf, eff_area, left=0.0, right=0.0)
            randvec = eff_area.max()*np.random.random(eobs.shape)
            detected = randvec < earea
        
        detected = np.logical_and(not_abs, detected)
                    
        events = {}

        dx_min = self.parameters["Width"]/self.parameters["Dimension"]
        dtheta = YTQuantity(np.rad2deg(dx_min/D_A), "degree")
        
        events["xpix"] = xsky[detected]/dx_min.v + 0.5*(nx+1)
        events["ypix"] = ysky[detected]/dx_min.v + 0.5*(nx+1)
        events["eobs"] = eobs[detected]

        events = comm.par_combine_object(events, datatype="dict", op="cat")

        if psf_sigma is not None:
            events["xpix"] += np.random.normal(sigma=psf_sigma/dtheta)
            events["ypix"] += np.random.normal(sigma=psf_sigma/dtheta)
        
        num_events = len(events["xpix"])
            
        if comm.rank == 0: mylog.info("Total number of observed photons: %d" % num_events)

        if "RMF" in parameters and convolve_energies:
            events, info = self._convolve_with_rmf(parameters["RMF"], events)
            for k, v in info.items(): parameters[k] = v
                
        if exp_time_new is None:
            parameters["ExposureTime"] = self.parameters["FiducialExposureTime"]
        else:
            parameters["ExposureTime"] = exp_time_new
        if area_new is None:
            parameters["Area"] = self.parameters["FiducialArea"]
        else:
            parameters["Area"] = area_new
        parameters["Redshift"] = zobs
        parameters["AngularDiameterDistance"] = D_A.in_units("Mpc")
        parameters["sky_center"] = sky_center
        parameters["pix_center"] = np.array([0.5*(nx+1)]*2)
        parameters["dtheta"] = dtheta
        
        return EventList(events, parameters)

    def _normalize_arf(self, respfile):
        rmf = _astropy.pyfits.open(respfile)
        table = rmf["MATRIX"]
        weights = np.array([w.sum() for w in table.data["MATRIX"]])
        rmf.close()
        return weights

    def _convolve_with_rmf(self, respfile, events):
        """
        Convolve the events with a RMF file.
        """
        mylog.warning("This routine has not been tested to work with all RMFs. YMMV.")
        mylog.info("Reading response matrix file (RMF): %s" % (respfile))
        
        hdulist = _astropy.pyfits.open(respfile)

        tblhdu = hdulist["MATRIX"]
        n_de = len(tblhdu.data["ENERG_LO"])
        mylog.info("Number of energy bins in RMF: %d" % (n_de))
        mylog.info("Energy limits: %g %g" % (min(tblhdu.data["ENERG_LO"]),
                                             max(tblhdu.data["ENERG_HI"])))

        tblhdu2 = hdulist["EBOUNDS"]
        n_ch = len(tblhdu2.data["CHANNEL"])
        mylog.info("Number of channels in RMF: %d" % (n_ch))
        
        eidxs = np.argsort(events["eobs"])

        phEE = events["eobs"][eidxs].ndarray_view()
        phXX = events["xpix"][eidxs]
        phYY = events["ypix"][eidxs]

        detectedChannels = []

        # run through all photon energies and find which bin they go in
        k = 0
        fcurr = 0
        last = len(phEE)-1

        pbar = get_pbar("Scattering energies with RMF:", n_de)
        
        for low,high in zip(tblhdu.data["ENERG_LO"],tblhdu.data["ENERG_HI"]):
            # weight function for probabilities from RMF
            weights = np.nan_to_num(tblhdu.data[k]["MATRIX"][:])
            weights /= weights.sum()
            # build channel number list associated to array value,
            # there are groups of channels in rmfs with nonzero probabilities
            trueChannel = []
            f_chan = np.nan_to_num(tblhdu.data["F_CHAN"][k])
            n_chan = np.nan_to_num(tblhdu.data["N_CHAN"][k])
            n_grp = np.nan_to_num(tblhdu.data["N_CHAN"][k])
            if not iterable(f_chan):
                f_chan = [f_chan]
                n_chan = [n_chan]
                n_grp  = [n_grp]
            for start,nchan in zip(f_chan, n_chan):
                end = start + nchan
                if start == end:
                    trueChannel.append(start)
                else:
                    for j in range(start,end):
                        trueChannel.append(j)
            if len(trueChannel) > 0:
                for q in range(fcurr,last):
                    if phEE[q] >= low and phEE[q] < high:
                        channelInd = np.random.choice(len(weights), p=weights)
                        fcurr +=1
                        detectedChannels.append(trueChannel[channelInd])
                    if phEE[q] >= high:
                        break
            pbar.update(k)
            k+=1
        pbar.finish()

        dchannel = np.array(detectedChannels)

        events["xpix"] = phXX
        events["ypix"] = phYY
        events["eobs"] = YTArray(phEE, "keV")
        events[tblhdu.header["CHANTYPE"]] = dchannel.astype(int)

        info = {"ChannelType" : tblhdu.header["CHANTYPE"],
                "Telescope" : tblhdu.header["TELESCOP"],
                "Instrument" : tblhdu.header["INSTRUME"]}

        info["Mission"] = tblhdu.header.get("MISSION","")

        return events, info

class EventList(object) :

    def __init__(self, events, parameters):

        self.events = events
        self.parameters = parameters
        self.num_events = events["xpix"].shape[0]
        self.wcs = _astropy.pywcs.WCS(naxis=2)
        self.wcs.wcs.crpix = parameters["pix_center"]
        self.wcs.wcs.crval = parameters["sky_center"].ndarray_view()
        self.wcs.wcs.cdelt = [-parameters["dtheta"].value, parameters["dtheta"].value]
        self.wcs.wcs.ctype = ["RA---TAN","DEC--TAN"]
        self.wcs.wcs.cunit = ["deg"]*2                                                
        x,y = self.wcs.wcs_pix2world(self.events["xpix"], self.events["ypix"], 1)
        self.events["xsky"] = YTArray(x, "degree")
        self.events["ysky"] = YTArray(y, "degree")

    def keys(self):
        return self.events.keys()

    def has_key(self, key):
        return key in self.keys()
    
    def items(self):
        return self.events.items()

    def values(self):
        return self.events.values()
    
    def __getitem__(self,key):                        
        return self.events[key]

    def __repr__(self):
        return self.events.__repr__()

    def __add__(self, other):
        keys1 = self.parameters.keys()
        keys2 = other.parameters.keys()
        keys1.sort()
        keys2.sort()
        if keys1 != keys2:
            raise RuntimeError("The two EventLists do not have the same parameters!")
        for k1, k2 in zip(keys1, keys2):
            v1 = self.parameters[k1]
            v2 = other.parameters[k2]
            if isinstance(v1, basestring) or isinstance(v2, basestring):
                check_equal = v1 == v2
            else:
                check_equal = np.allclose(v1, v2, rtol=0.0, atol=1.0e-10)
            if not check_equal:
                raise RuntimeError("The values for the parameter '%s' in the two EventLists" % k1 +
                                   " are not identical (%s vs. %s)!" % (v1, v2))
        events = {}
        for item1, item2 in zip(self.items(), other.items()):
            k1, v1 = item1
            k2, v2 = item2
            events[k1] = uconcatenate([v1,v2])
        return EventList(events, self.parameters)

    def filter_events(self, region):
        """                                                                                                                                 
        Filter events using a ds9 region. Requires the pyregion package.                                                                    
        Returns a new EventList.                                                                                                            
        """
        import pyregion
        import os
        if os.path.exists(region):
            reg = pyregion.open(region)
        else:
            reg = pyregion.parse(region)
        r = reg.as_imagecoord(header=self.wcs.to_header())
        f = r.get_filter()
        idxs = f.inside_x_y(self.events["xpix"], self.events["ypix"])
        if idxs.sum() == 0:
            raise RuntimeError("No events are inside this region!")
        new_events = {}
        for k, v in self.events.items():
            new_events[k] = v[idxs]
        return EventList(new_events, self.parameters)
   
    @classmethod
    def from_h5_file(cls, h5file):
        """
        Initialize an EventList from a HDF5 file with filename *h5file*.
        """
        events = {}
        parameters = {}
        
        f = h5py.File(h5file, "r")

        parameters["ExposureTime"] = YTQuantity(f["/exp_time"].value, "s")
        parameters["Area"] = YTQuantity(f["/area"].value, "cm**2")
        parameters["Redshift"] = f["/redshift"].value
        parameters["AngularDiameterDistance"] = YTQuantity(f["/d_a"].value, "Mpc")
        if "rmf" in f:
            parameters["RMF"] = f["/rmf"].value
        if "arf" in f:
            parameters["ARF"] = f["/arf"].value
        if "channel_type" in f:
            parameters["ChannelType"] = f["/channel_type"].value
        if "mission" in f:
            parameters["Mission"] = f["/mission"].value
        if "telescope" in f:
            parameters["Telescope"] = f["/telescope"].value
        if "instrument" in f:
            parameters["Instrument"] = f["/instrument"].value

        events["xpix"] = f["/xpix"][:]
        events["ypix"] = f["/ypix"][:]
        events["eobs"] = YTArray(f["/eobs"][:], "keV")
        if "pi" in f:
            events["PI"] = f["/pi"][:]
        if "pha" in f:
            events["PHA"] = f["/pha"][:]
        parameters["sky_center"] = YTArray(f["/sky_center"][:], "deg")
        parameters["dtheta"] = YTQuantity(f["/dtheta"].value, "deg")
        parameters["pix_center"] = f["/pix_center"][:]
        
        f.close()
        
        return cls(events, parameters)

    @classmethod
    def from_fits_file(cls, fitsfile):
        """
        Initialize an EventList from a FITS file with filename *fitsfile*.
        """
        hdulist = _astropy.pyfits.open(fitsfile)

        tblhdu = hdulist["EVENTS"]
        
        events = {}
        parameters = {}
        
        parameters["ExposureTime"] = YTQuantity(tblhdu.header["EXPOSURE"], "s")
        parameters["Area"] = YTQuantity(tblhdu.header["AREA"], "cm**2")
        parameters["Redshift"] = tblhdu.header["REDSHIFT"]
        parameters["AngularDiameterDistance"] = YTQuantity(tblhdu.header["D_A"], "Mpc")
        if "RMF" in tblhdu.header:
            parameters["RMF"] = tblhdu["RMF"]
        if "ARF" in tblhdu.header:
            parameters["ARF"] = tblhdu["ARF"]
        if "CHANTYPE" in tblhdu.header:
            parameters["ChannelType"] = tblhdu["CHANTYPE"]
        if "MISSION" in tblhdu.header:
            parameters["Mission"] = tblhdu["MISSION"]
        if "TELESCOP" in tblhdu.header:
            parameters["Telescope"] = tblhdu["TELESCOP"]
        if "INSTRUME" in tblhdu.header:
            parameters["Instrument"] = tblhdu["INSTRUME"]
        parameters["sky_center"] = YTArray([tblhdu["TCRVL2"],tblhdu["TCRVL3"]], "deg")
        parameters["pix_center"] = np.array([tblhdu["TCRVL2"],tblhdu["TCRVL3"]])
        parameters["dtheta"] = YTQuantity(tblhdu["TCRVL3"], "deg")
        events["xpix"] = tblhdu.data.field("X")
        events["ypix"] = tblhdu.data.field("Y")
        events["eobs"] = YTArray(tblhdu.data.field("ENERGY")/1000., "keV")
        if "PI" in tblhdu.columns.names:
            events["PI"] = tblhdu.data.field("PI")
        if "PHA" in tblhdu.columns.names:
            events["PHA"] = tblhdu.data.field("PHA")
        
        return cls(events, parameters)

    @parallel_root_only
    def write_fits_file(self, fitsfile, clobber=False):
        """
        Write events to a FITS binary table file with filename *fitsfile*.
        Set *clobber* to True if you need to overwrite a previous file.
        """
        pyfits = _astropy.pyfits
        
        cols = []

        col1 = pyfits.Column(name='ENERGY', format='E', unit='eV',
                             array=self.events["eobs"].in_units("eV").ndarray_view())
        col2 = pyfits.Column(name='X', format='D', unit='pixel',
                             array=self.events["xpix"])
        col3 = pyfits.Column(name='Y', format='D', unit='pixel',
                             array=self.events["ypix"])

        cols = [col1, col2, col3]

        if "ChannelType" in self.parameters:
             chantype = self.parameters["ChannelType"]
             if chantype == "PHA":
                  cunit="adu"
             elif chantype == "PI":
                  cunit="Chan"
             col4 = pyfits.Column(name=chantype.upper(), format='1J',
                                  unit=cunit, array=self.events[chantype])
             cols.append(col4)
            
        coldefs = pyfits.ColDefs(cols)
        tbhdu = pyfits.new_table(coldefs)
        tbhdu.update_ext_name("EVENTS")

        tbhdu.header["MTYPE1"] = "sky"
        tbhdu.header["MFORM1"] = "x,y"
        tbhdu.header["MTYPE2"] = "EQPOS"
        tbhdu.header["MFORM2"] = "RA,DEC"
        tbhdu.header["TCTYP2"] = "RA---TAN"
        tbhdu.header["TCTYP3"] = "DEC--TAN"
        tbhdu.header["TCRVL2"] = float(self.parameters["sky_center"][0])
        tbhdu.header["TCRVL3"] = float(self.parameters["sky_center"][1])
        tbhdu.header["TCDLT2"] = -float(self.parameters["dtheta"])
        tbhdu.header["TCDLT3"] = float(self.parameters["dtheta"])
        tbhdu.header["TCRPX2"] = self.parameters["pix_center"][0]
        tbhdu.header["TCRPX3"] = self.parameters["pix_center"][1]
        tbhdu.header["TLMIN2"] = 0.5
        tbhdu.header["TLMIN3"] = 0.5
        tbhdu.header["TLMAX2"] = 2.*self.parameters["pix_center"][0]-0.5
        tbhdu.header["TLMAX3"] = 2.*self.parameters["pix_center"][1]-0.5
        tbhdu.header["EXPOSURE"] = float(self.parameters["ExposureTime"])
        if isinstance(self.parameters["Area"], basestring):
            tbhdu.header["AREA"] = self.parameters["Area"]
        else:
            tbhdu.header["AREA"] = float(self.parameters["Area"])
        tbhdu.header["D_A"] = float(self.parameters["AngularDiameterDistance"])
        tbhdu.header["REDSHIFT"] = self.parameters["Redshift"]
        tbhdu.header["HDUVERS"] = "1.1.0"
        tbhdu.header["RADECSYS"] = "FK5"
        tbhdu.header["EQUINOX"] = 2000.0
        if "RMF" in self.parameters:
            tbhdu.header["RESPFILE"] = self.parameters["RMF"]
        if "ARF" in self.parameters:
            tbhdu.header["ANCRFILE"] = self.parameters["ARF"]
        if "ChannelType" in self.parameters:
            tbhdu.header["CHANTYPE"] = self.parameters["ChannelType"]
        if "Mission" in self.parameters:
            tbhdu.header["MISSION"] = self.parameters["Mission"]
        if "Telescope" in self.parameters:
            tbhdu.header["TELESCOP"] = self.parameters["Telescope"]
        if "Instrument" in self.parameters:
            tbhdu.header["INSTRUME"] = self.parameters["Instrument"]
            
        tbhdu.writeto(fitsfile, clobber=clobber)

    @parallel_root_only    
    def write_simput_file(self, prefix, clobber=False, emin=None, emax=None):
        r"""
        Write events to a SIMPUT file that may be read by the SIMX instrument
        simulator.

        Parameters
        ----------

        prefix : string
            The filename prefix.
        clobber : boolean, optional
            Set to True to overwrite previous files.
        e_min : float, optional
            The minimum energy of the photons to save in keV.
        e_max : float, optional
            The maximum energy of the photons to save in keV.
        """
        pyfits = _astropy.pyfits
        if isinstance(self.parameters["Area"], basestring):
             mylog.error("Writing SIMPUT files is only supported if you didn't convolve with an ARF.")
             raise TypeError("Writing SIMPUT files is only supported if you didn't convolve with an ARF.")
        
        if emin is None:
            emin = self.events["eobs"].min().value*0.95
        if emax is None:
            emax = self.events["eobs"].max().value*1.05

        idxs = np.logical_and(self.events["eobs"].ndarray_view() >= emin,
                              self.events["eobs"].ndarray_view() <= emax)
        flux = np.sum(self.events["eobs"][idxs].in_units("erg")) / \
               self.parameters["ExposureTime"]/self.parameters["Area"]
        
        col1 = pyfits.Column(name='ENERGY', format='E',
                             array=self["eobs"].ndarray_view())
        col2 = pyfits.Column(name='DEC', format='D',
                             array=self["ysky"].ndarray_view())
        col3 = pyfits.Column(name='RA', format='D',
                             array=self["xsky"].ndarray_view())

        coldefs = pyfits.ColDefs([col1, col2, col3])

        tbhdu = pyfits.new_table(coldefs)
        tbhdu.update_ext_name("PHLIST")

        tbhdu.header["HDUCLASS"] = "HEASARC/SIMPUT"
        tbhdu.header["HDUCLAS1"] = "PHOTONS"
        tbhdu.header["HDUVERS"] = "1.1.0"
        tbhdu.header["EXTVER"] = 1
        tbhdu.header["REFRA"] = 0.0
        tbhdu.header["REFDEC"] = 0.0
        tbhdu.header["TUNIT1"] = "keV"
        tbhdu.header["TUNIT2"] = "deg"
        tbhdu.header["TUNIT3"] = "deg"

        phfile = prefix+"_phlist.fits"

        tbhdu.writeto(phfile, clobber=clobber)

        col1 = pyfits.Column(name='SRC_ID', format='J', array=np.array([1]).astype("int32"))
        col2 = pyfits.Column(name='RA', format='D', array=np.array([float(self.parameters["sky_center"][0])]))
        col3 = pyfits.Column(name='DEC', format='D', array=np.array([float(self.parameters["sky_center"][1])]))
        col4 = pyfits.Column(name='E_MIN', format='D', array=np.array([float(emin)]))
        col5 = pyfits.Column(name='E_MAX', format='D', array=np.array([float(emax)]))
        col6 = pyfits.Column(name='FLUX', format='D', array=np.array([flux.value]))
        col7 = pyfits.Column(name='SPECTRUM', format='80A', array=np.array([phfile+"[PHLIST,1]"]))
        col8 = pyfits.Column(name='IMAGE', format='80A', array=np.array([phfile+"[PHLIST,1]"]))
                        
        coldefs = pyfits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8])
        
        wrhdu = pyfits.new_table(coldefs)
        wrhdu.update_ext_name("SRC_CAT")
                                
        wrhdu.header["HDUCLASS"] = "HEASARC"
        wrhdu.header["HDUCLAS1"] = "SIMPUT"
        wrhdu.header["HDUCLAS2"] = "SRC_CAT"
        wrhdu.header["HDUVERS"] = "1.1.0"
        wrhdu.header["RADECSYS"] = "FK5"
        wrhdu.header["EQUINOX"] = 2000.0
        wrhdu.header["TUNIT2"] = "deg"
        wrhdu.header["TUNIT3"] = "deg"
        wrhdu.header["TUNIT4"] = "keV"
        wrhdu.header["TUNIT5"] = "keV"
        wrhdu.header["TUNIT6"] = "erg/s/cm**2"

        simputfile = prefix+"_simput.fits"
                
        wrhdu.writeto(simputfile, clobber=clobber)

    @parallel_root_only
    def write_h5_file(self, h5file):
        """
        Write an EventList to the HDF5 file given by *h5file*.
        """
        f = h5py.File(h5file, "w")

        f.create_dataset("/exp_time", data=float(self.parameters["ExposureTime"]))
        area = self.parameters["Area"]
        if not isinstance(area, basestring):
            area = float(area)
        f.create_dataset("/area", data=area)
        f.create_dataset("/redshift", data=self.parameters["Redshift"])
        f.create_dataset("/d_a", data=float(self.parameters["AngularDiameterDistance"]))
        if "ARF" in self.parameters:
            f.create_dataset("/arf", data=self.parameters["ARF"])
        if "RMF" in self.parameters:
            f.create_dataset("/rmf", data=self.parameters["RMF"])
        if "ChannelType" in self.parameters:
            f.create_dataset("/channel_type", data=self.parameters["ChannelType"])
        if "Mission" in self.parameters:
            f.create_dataset("/mission", data=self.parameters["Mission"]) 
        if "Telescope" in self.parameters:
            f.create_dataset("/telescope", data=self.parameters["Telescope"])
        if "Instrument" in self.parameters:
            f.create_dataset("/instrument", data=self.parameters["Instrument"])
                            
        f.create_dataset("/xpix", data=self.events["xpix"])
        f.create_dataset("/ypix", data=self.events["ypix"])
        f.create_dataset("/xsky", data=self.events["xsky"].ndarray_view())
        f.create_dataset("/ysky", data=self.events["ysky"].ndarray_view())
        f.create_dataset("/eobs", data=self.events["eobs"].ndarray_view())
        if "PI" in self.events:
            f.create_dataset("/pi", data=self.events["PI"])                  
        if "PHA" in self.events:
            f.create_dataset("/pha", data=self.events["PHA"])                  
        f.create_dataset("/sky_center", data=self.parameters["sky_center"].ndarray_view())
        f.create_dataset("/pix_center", data=self.parameters["pix_center"])
        f.create_dataset("/dtheta", data=float(self.parameters["dtheta"]))

        f.close()

    @parallel_root_only
    def write_fits_image(self, imagefile, clobber=False,
                         emin=None, emax=None):
        r"""
        Generate a image by binning X-ray counts and write it to a FITS file.

        Parameters
        ----------

        imagefile : string
            The name of the image file to write.
        clobber : boolean, optional
            Set to True to overwrite a previous file.
        emin : float, optional
            The minimum energy of the photons to put in the image, in keV.
        emax : float, optional
            The maximum energy of the photons to put in the image, in keV.
        """
        if emin is None:
            mask_emin = np.ones((self.num_events), dtype='bool')
        else:
            mask_emin = self.events["eobs"].ndarray_view() > emin
        if emax is None:
            mask_emax = np.ones((self.num_events), dtype='bool')
        else:
            mask_emax = self.events["eobs"].ndarray_view() < emax

        mask = np.logical_and(mask_emin, mask_emax)

        nx = int(2*self.parameters["pix_center"][0]-1.)
        ny = int(2*self.parameters["pix_center"][1]-1.)
        
        xbins = np.linspace(0.5, float(nx)+0.5, nx+1, endpoint=True)
        ybins = np.linspace(0.5, float(ny)+0.5, ny+1, endpoint=True)

        H, xedges, yedges = np.histogram2d(self.events["xpix"][mask],
                                           self.events["ypix"][mask],
                                           bins=[xbins,ybins])
        
        hdu = _astropy.pyfits.PrimaryHDU(H.T)
        
        hdu.header["MTYPE1"] = "EQPOS"
        hdu.header["MFORM1"] = "RA,DEC"
        hdu.header["CTYPE1"] = "RA---TAN"
        hdu.header["CTYPE2"] = "DEC--TAN"
        hdu.header["CRPIX1"] = 0.5*(nx+1)
        hdu.header["CRPIX2"] = 0.5*(nx+1)
        hdu.header["CRVAL1"] = float(self.parameters["sky_center"][0])
        hdu.header["CRVAL2"] = float(self.parameters["sky_center"][1])
        hdu.header["CUNIT1"] = "deg"
        hdu.header["CUNIT2"] = "deg"
        hdu.header["CDELT1"] = -float(self.parameters["dtheta"])
        hdu.header["CDELT2"] = float(self.parameters["dtheta"])
        hdu.header["EXPOSURE"] = float(self.parameters["ExposureTime"])
        
        hdu.writeto(imagefile, clobber=clobber)
                                    
    @parallel_root_only
    def write_spectrum(self, specfile, energy_bins=False, emin=0.1,
                       emax=10.0, nchan=2000, clobber=False):
        r"""
        Bin event energies into a spectrum and write it to a FITS binary table. Can bin
        on energy or channel, the latter only if the photons were convolved with an
        RMF. In that case, the spectral binning will be determined by the RMF binning.

        Parameters
        ----------

        specfile : string
            The name of the FITS file to be written.
        energy_bins : boolean, optional
            Bin on energy or channel. 
        emin : float, optional
            The minimum energy of the spectral bins in keV. Only used if binning on energy.
        emax : float, optional
            The maximum energy of the spectral bins in keV. Only used if binning on energy.
        nchan : integer, optional
            The number of channels. Only used if binning on energy.
        """
        pyfits = _astropy.pyfits
        if energy_bins:
            spectype = "energy"
            espec = self.events["eobs"].d
            range = (emin, emax)
            spec, ee = np.histogram(espec, bins=nchan, range=range)
            bins = 0.5*(ee[1:]+ee[:-1])
        else:
            if "ChannelType" not in self.parameters:
                mylog.error("These events were not convolved with an RMF. Set energy_bins=True.")
                raise KeyError
            spectype = self.parameters["ChannelType"]
            espec = self.events[spectype]
            f = pyfits.open(self.parameters["RMF"])
            nchan = int(f[1].header["DETCHANS"])
            try:
                cmin = int(f[1].header["TLMIN4"])
            except KeyError:
                mylog.warning("Cannot determine minimum allowed value for channel. " +
                              "Setting to 0, which may be wrong.")
                cmin = int(0)
            try:
                cmax = int(f[1].header["TLMAX4"])
            except KeyError:
                mylog.warning("Cannot determine maximum allowed value for channel. " +
                              "Setting to DETCHANS, which may be wrong.")
                cmax = int(nchan)
            f.close()
            minlength = nchan
            if cmin == 1: minlength += 1
            spec = np.bincount(self.events[spectype],minlength=minlength)
            if cmin == 1: spec = spec[1:]
            bins = np.arange(nchan).astype("int32")+cmin
            
        col1 = pyfits.Column(name='CHANNEL', format='1J', array=bins)
        col2 = pyfits.Column(name=spectype.upper(), format='1D', array=bins.astype("float64"))
        col3 = pyfits.Column(name='COUNTS', format='1J', array=spec.astype("int32"))
        col4 = pyfits.Column(name='COUNT_RATE', format='1D', array=spec/float(self.parameters["ExposureTime"]))
        
        coldefs = pyfits.ColDefs([col1, col2, col3, col4])
        
        tbhdu = pyfits.new_table(coldefs)
        tbhdu.update_ext_name("SPECTRUM")

        if not energy_bins:
            tbhdu.header["DETCHANS"] = spec.shape[0]
            tbhdu.header["TOTCTS"] = spec.sum()
            tbhdu.header["EXPOSURE"] = float(self.parameters["ExposureTime"])
            tbhdu.header["LIVETIME"] = float(self.parameters["ExposureTime"])
            tbhdu.header["CONTENT"] = spectype
            tbhdu.header["HDUCLASS"] = "OGIP"
            tbhdu.header["HDUCLAS1"] = "SPECTRUM"
            tbhdu.header["HDUCLAS2"] = "TOTAL"
            tbhdu.header["HDUCLAS3"] = "TYPE:I"
            tbhdu.header["HDUCLAS4"] = "COUNT"
            tbhdu.header["HDUVERS"] = "1.1.0"
            tbhdu.header["HDUVERS1"] = "1.1.0"
            tbhdu.header["CHANTYPE"] = spectype
            tbhdu.header["BACKFILE"] = "none"
            tbhdu.header["CORRFILE"] = "none"
            tbhdu.header["POISSERR"] = True
            if "RMF" in self.parameters:
                 tbhdu.header["RESPFILE"] = self.parameters["RMF"]
            else:
                 tbhdu.header["RESPFILE"] = "none"
            if "ARF" in self.parameters:
                tbhdu.header["ANCRFILE"] = self.parameters["ARF"]
            else:        
                tbhdu.header["ANCRFILE"] = "none"
            if "Mission" in self.parameters:
                tbhdu.header["MISSION"] = self.parameters["Mission"]
            else:
                tbhdu.header["MISSION"] = "none"
            if "Telescope" in self.parameters:
                tbhdu.header["TELESCOP"] = self.parameters["Telescope"]
            else:
                tbhdu.header["TELESCOP"] = "none"
            if "Instrument" in self.parameters:
                tbhdu.header["INSTRUME"] = self.parameters["Instrument"]
            else:
                tbhdu.header["INSTRUME"] = "none"
            tbhdu.header["AREASCAL"] = 1.0
            tbhdu.header["CORRSCAL"] = 0.0
            tbhdu.header["BACKSCAL"] = 1.0
                                
        hdulist = pyfits.HDUList([pyfits.PrimaryHDU(), tbhdu])
        
        hdulist.writeto(specfile, clobber=clobber)
