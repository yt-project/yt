"""
Classes for generating lists of photons and detected events
The algorithms used here are based off of the method used by the
PHOX code (http://www.mpa-garching.mpg.de/~kdolag/Phox/),
developed by Veronica Biffi and Klaus Dolag. References for
PHOX may be found at:

Biffi et al 2012: http://adsabs.harvard.edu/abs/2012MNRAS.420.3545B
Biffi et al 2013: http://adsabs.harvard.edu/abs/2013MNRAS.428.1395B
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from numpy.testing import assert_allclose
from yt.funcs import *
from yt.utilities.physical_constants import clight, \
     cm_per_km, erg_per_keV
from yt.utilities.cosmology import Cosmology
from yt.utilities.orientation import Orientation
from yt.utilities.parallel_tools.parallel_analysis_interface import \
     communication_system, parallel_root_only, get_mpi_type, \
     op_names, parallel_capable

import h5py

try:
    import astropy.io.fits as pyfits
    import astropy.wcs as pywcs
except ImportError:
    pass

comm = communication_system.communicators[-1]
    
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
    
    @classmethod
    def from_file(cls, filename):
        r"""
        Initialize a PhotonList from the HDF5 file *filename*.
        """

        photons = {}
        parameters = {}
        
        f = h5py.File(filename, "r")

        parameters["FiducialExposureTime"] = f["/fid_exp_time"].value
        parameters["FiducialArea"] = f["/fid_area"].value
        parameters["FiducialRedshift"] = f["/fid_redshift"].value
        parameters["FiducialAngularDiameterDistance"] = f["/fid_d_a"].value
        parameters["Dimension"] = f["/dimension"].value
        parameters["Width"] = f["/width"].value
        parameters["HubbleConstant"] = f["/hubble"].value
        parameters["OmegaMatter"] = f["/omega_matter"].value
        parameters["OmegaLambda"] = f["/omega_lambda"].value

        num_cells = f["/x"][:].shape[0]
        start_c = comm.rank*num_cells/comm.size
        end_c = (comm.rank+1)*num_cells/comm.size
        
        photons["x"] = f["/x"][start_c:end_c]
        photons["y"] = f["/y"][start_c:end_c]
        photons["z"] = f["/z"][start_c:end_c]
        photons["dx"] = f["/dx"][start_c:end_c]
        photons["vx"] = f["/vx"][start_c:end_c]
        photons["vy"] = f["/vy"][start_c:end_c]
        photons["vz"] = f["/vz"][start_c:end_c]

        n_ph = f["/num_photons"][:]
        
        if comm.rank == 0:
            start_e = np.uint64(0)
        else:
            start_e = n_ph[:start_c].sum()
        end_e = start_e + np.uint64(n_ph[start_c:end_c].sum())

        photons["NumberOfPhotons"] = n_ph[start_c:end_c]

        p_bins = np.cumsum(photons["NumberOfPhotons"])
        p_bins = np.insert(p_bins, 0, [np.uint64(0)])
        
        photons["Energy"] = f["/energy"][start_e:end_e]
        
        f.close()

        cosmo = Cosmology(HubbleConstantNow=parameters["HubbleConstant"],
                          OmegaMatterNow=parameters["OmegaMatter"],
                          OmegaLambdaNow=parameters["OmegaLambda"])

        return cls(photons, parameters, cosmo, p_bins)
    
    @classmethod
    def from_scratch(cls, data_source, redshift, area,
                     exp_time, photon_model, parameters=None,
                     center=None, dist=None, cosmology=None):
        """
        Initialize a PhotonList from a photon model. The redshift, collecting area,
        exposure time, and cosmology are stored in the *parameters* dictionary which
        is passed to the *photon_model* function. 

        Parameters
        ----------

        data_source : `yt.data_objects.api.AMRData`
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
            Cosmological information. If not supplied, it assumes \LambdaCDM with
            the default yt parameters.

        Examples
        --------

        This is a simple example where a point source with a single line emission
        spectrum of photons is created. More complicated examples which actually
        create photons based on the fields in the dataset could be created. 

        >>> from scipy.stats import powerlaw
        >>> def line_func(source, photons, parameters):
        ...
        ...     pf = source.pf
        ... 
        ...     num_photons = parameters["num_photons"]
        ...     E0  = parameters["line_energy"] # Energies are in keV
        ...     sigE = parameters["line_sigma"] 
        ...
        ...     energies = norm.rvs(loc=E0, scale=sigE, size=num_photons)
        ...     
        ...     photons["x"] = np.zeros((1)) # Place everything in the center cell
        ...     photons["y"] = np.zeros((1))
        ...     photons["z"] = np.zeros((1))
        ...     photons["vx"] = np.zeros((1))
        ...     photons["vy"] = np.zeros((1))
        ...     photons["vz"] = 100.*np.ones((1))
        ...     photons["dx"] = source["dx"][0]*pf.units["kpc"]*np.ones((1)) 
        ...     photons["NumberOfPhotons"] = num_photons*np.ones((1))
        ...     photons["Energy"] = np.array(energies)
        >>>
        >>> redshift = 0.05
        >>> area = 6000.0
        >>> time = 2.0e5
        >>> parameters = {"num_photons" : 10000, "line_energy" : 5.0,
        ...               "line_sigma" : 0.1}
        >>> ddims = (128,128,128)
        >>> random_data = {"Density":np.random.random(ddims)}
        >>> pf = load_uniform_grid(random_data, ddims)
        >>> dd = pf.h.all_data
        >>> my_photons = PhotonList.from_user_model(dd, redshift, area,
        ...                                         time, line_func)

        """

        pf = data_source.pf

        if parameters is None:
             parameters = {}
        if cosmology is None:
            cosmo = Cosmology()
        else:
            cosmo = cosmology
        if dist is None:
            D_A = cosmo.AngularDiameterDistance(0.0,redshift)
        else:
            D_A = dist[0]*pf.units["mpc"]/pf.units[dist[1]]
            redshift = 0.0

        if center == "c":
            parameters["center"] = pf.domain_center
        elif center == "max":
            parameters["center"] = pf.h.find_max("Density")[-1]
        elif iterable(center):
            parameters["center"] = center
        elif center is None:
            parameters["center"] = data_source.get_field_parameter("center")
            
        parameters["FiducialExposureTime"] = exp_time
        parameters["FiducialArea"] = area
        parameters["FiducialRedshift"] = redshift
        parameters["FiducialAngularDiameterDistance"] = D_A
        parameters["HubbleConstant"] = cosmo.HubbleConstantNow
        parameters["OmegaMatter"] = cosmo.OmegaMatterNow
        parameters["OmegaLambda"] = cosmo.OmegaLambdaNow

        dimension = 0
        width = 0.0
        for i, ax in enumerate("xyz"):
            pos = data_source[ax]
            delta = data_source["d%s"%(ax)]
            le = np.min(pos-0.5*delta)
            re = np.max(pos+0.5*delta)
            width = 2.*max(width, re-parameters["center"][i], parameters["center"][i]-le)
            dimension = max(dimension, int(width/delta.min()))
        parameters["Dimension"] = dimension
        parameters["Width"] = width*pf.units["kpc"]
                
        photons = photon_model(data_source, parameters)
        
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
                x = np.zeros((num_cells))
                y = np.zeros((num_cells))
                z = np.zeros((num_cells))
                vx = np.zeros((num_cells))
                vy = np.zeros((num_cells))
                vz = np.zeros((num_cells))
                dx = np.zeros((num_cells))
                n_ph = np.zeros((num_cells), dtype="uint64")
                e = np.zeros((num_photons))
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
                                                
            comm.comm.Gatherv([self.photons["x"], local_num_cells, mpi_double],
                              [x, (sizes_c, disps_c), mpi_double], root=0)
            comm.comm.Gatherv([self.photons["y"], local_num_cells, mpi_double],
                              [y, (sizes_c, disps_c), mpi_double], root=0)
            comm.comm.Gatherv([self.photons["z"], local_num_cells, mpi_double],
                              [z, (sizes_c, disps_c), mpi_double], root=0)
            comm.comm.Gatherv([self.photons["vx"], local_num_cells, mpi_double],
                              [vx, (sizes_c, disps_c), mpi_double], root=0)
            comm.comm.Gatherv([self.photons["vy"], local_num_cells, mpi_double],
                              [vy, (sizes_c, disps_c), mpi_double], root=0)
            comm.comm.Gatherv([self.photons["vz"], local_num_cells, mpi_double],
                              [vz, (sizes_c, disps_c), mpi_double], root=0)
            comm.comm.Gatherv([self.photons["dx"], local_num_cells, mpi_double],
                              [dx, (sizes_c, disps_c), mpi_double], root=0)
            comm.comm.Gatherv([self.photons["NumberOfPhotons"], local_num_cells, mpi_long],
                              [n_ph, (sizes_c, disps_c), mpi_long], root=0)
            comm.comm.Gatherv([self.photons["Energy"], local_num_photons, mpi_double],
                              [e, (sizes_p, disps_p), mpi_double], root=0) 

        else:

            x = self.photons["x"]
            y = self.photons["y"]
            z = self.photons["z"]
            vx = self.photons["vx"]
            vy = self.photons["vy"]
            vz = self.photons["vz"]
            dx = self.photons["dx"]
            n_ph = self.photons["NumberOfPhotons"]
            e = self.photons["Energy"]
                                                
        if comm.rank == 0:
            
            f = h5py.File(photonfile, "w")

            # Scalars
       
            f.create_dataset("fid_area", data=self.parameters["FiducialArea"])
            f.create_dataset("fid_exp_time", data=self.parameters["FiducialExposureTime"])
            f.create_dataset("fid_redshift", data=self.parameters["FiducialRedshift"])
            f.create_dataset("hubble", data=self.parameters["HubbleConstant"])
            f.create_dataset("omega_matter", data=self.parameters["OmegaMatter"])
            f.create_dataset("omega_lambda", data=self.parameters["OmegaLambda"])
            f.create_dataset("fid_d_a", data=self.parameters["FiducialAngularDiameterDistance"])
            f.create_dataset("dimension", data=self.parameters["Dimension"])
            f.create_dataset("width", data=self.parameters["Width"])
                        
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

    def project_photons(self, L, area_new=None, texp_new=None, 
                        redshift_new=None, dist_new=None,
                        absorb_model=None, psf_sigma=None,
                        sky_center=None):
        r"""
        Projects photons onto an image plane given a line of sight.

        Parameters
        ----------

        L : array_like
            Normal vector to the plane of projection.
        area_new : float or filename, optional
            New value for the effective area of the detector. 
            Either a single float value or a standard ARF file
            containing the effective area as a function of energy.
        texp_new : float, optional
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
             sky_center = np.array([30.,45.])

        dx = self.photons["dx"]
        nx = self.parameters["Dimension"]
        if psf_sigma is not None:
             psf_sigma /= 3600.

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
        num_cells = len(n_ph)
        n_ph_tot = n_ph.sum()
        
        eff_area = None
        
        if (texp_new is None and area_new is None and
            redshift_new is None and dist_new is None):
            my_n_obs = n_ph_tot
            zobs = self.parameters["FiducialRedshift"]
            D_A = self.parameters["FiducialAngularDiameterDistance"]*1000.
        else:
            if texp_new is None:
                Tratio = 1.
            else:
                Tratio = texp_new/self.parameters["FiducialExposureTime"]
            if area_new is None:
                Aratio = 1.
            elif isinstance(area_new, basestring):
                if comm.rank == 0:
                    mylog.info("Using energy-dependent effective area.")
                f = pyfits.open(area_new)
                elo = f["SPECRESP"].data.field("ENERG_LO")
                ehi = f["SPECRESP"].data.field("ENERG_HI")
                eff_area = np.nan_to_num(f["SPECRESP"].data.field("SPECRESP"))
                f.close()
                Aratio = eff_area.max()/self.parameters["FiducialArea"]
            else:
                mylog.info("Using constant effective area.")
                Aratio = area_new/self.parameters["FiducialArea"]
            if redshift_new is None and dist_new is None:
                Dratio = 1.
                zobs = self.parameters["FiducialRedshift"]
                D_A = self.parameters["FiducialAngularDiameterDistance"]*1000.                    
            else:
                if redshift_new is None:
                    zobs = 0.0
                    D_A = dist[0]*self.pf.units["kpc"]/self.pf.units[dist[1]]
                else:
                    zobs = redshift_new
                    D_A = self.cosmo.AngularDiameterDistance(0.0,zobs)*1000.
                fid_D_A = self.parameters["FiducialAngularDiameterDistance"]*1000.
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
                    
        vz = self.photons["vx"]*z_hat[0] + \
             self.photons["vy"]*z_hat[1] + \
             self.photons["vz"]*z_hat[2]
        shift = -vz*cm_per_km/clight
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
        dtheta = np.rad2deg(dx_min/D_A)
        
        events["xpix"] = xsky[detected]/dx_min + 0.5*(nx+1) 
        events["ypix"] = ysky[detected]/dx_min + 0.5*(nx+1)
        events["eobs"] = eobs[detected]

        events = comm.par_combine_object(events, datatype="dict", op="cat")

        if psf_sigma is not None:
            events["xpix"] += np.random.normal(sigma=psf_sigma/dtheta)
            events["ypix"] += np.random.normal(sigma=psf_sigma/dtheta)
        
        num_events = len(events["xpix"])
            
        if comm.rank == 0: mylog.info("Total number of observed photons: %d" % (num_events))

        parameters = {}
        
        if texp_new is None:
            parameters["ExposureTime"] = self.parameters["FiducialExposureTime"]
        else:
            parameters["ExposureTime"] = texp_new
        if area_new is None:
            parameters["Area"] = self.parameters["FiducialArea"]
        else:
            parameters["Area"] = area_new
        parameters["Redshift"] = zobs
        parameters["AngularDiameterDistance"] = D_A/1000.
        if isinstance(area_new, basestring):
            parameters["ARF"] = area_new
        parameters["sky_center"] = np.array(sky_center)
        parameters["pix_center"] = np.array([0.5*(nx+1)]*2)
        parameters["dtheta"] = dtheta
        
        return EventList(events, parameters)

class EventList(object) :

    def __init__(self, events, parameters) :

        self.events = events
        self.parameters = parameters
        self.num_events = events["xpix"].shape[0]
        self.wcs = pywcs.WCS(naxis=2)
        self.wcs.wcs.crpix = parameters["pix_center"]
        self.wcs.wcs.crval = parameters["sky_center"]
        self.wcs.wcs.cdelt = [-parameters["dtheta"], parameters["dtheta"]]
        self.wcs.wcs.ctype = ["RA---TAN","DEC--TAN"]
        self.wcs.wcs.cunit = ["deg"]*2                                                
        
    def keys(self):
        return self.events.keys()

    def has_key(self, key):
        return key in self.keys()
    
    def items(self):
        return self.events.items()

    def values(self):
        return self.events.values()
    
    def __getitem__(self,key):

        if key == "xsky" or key == "ysky":
            if not self.has_key(key):
                (self.events["xsky"],
                 self.events["ysky"]) = \
                 self.wcs.wcs_pix2world(events["xpix"], events["ypix"],
                                        1, ra_dec_order=True)
                        
        return self.events[key]
        
    @classmethod
    def from_h5_file(cls, h5file):
        """
        Initialize an EventList from a HDF5 file with filename *h5file*.
        """
        events = {}
        parameters = {}
        
        f = h5py.File(h5file, "r")

        parameters["ExposureTime"] = f["/exp_time"].value
        parameters["Area"] = f["/area"].value
        parameters["Redshift"] = f["/redshift"].value
        parameters["AngularDiameterDistance"] = f["/d_a"].value
        if "rmf" in f:
            parameters["RMF"] = f["/rmf"].value
        if "arf" in f:
            parameters["ARF"] = f["/arf"].value
        if "channel_type" in f:
            parameters["ChannelType"] = f["/channel_type"].value
        if "telescope" in f:
            parameters["Telescope"] = f["/telescope"].value
        if "instrument" in f:
            parameters["Instrument"] = f["/instrument"].value

        events["xpix"] = f["/xpix"][:]
        events["ypix"] = f["/ypix"][:]
        events["eobs"] = f["/eobs"][:]
        if "pi" in f:
            events["PI"] = f["/pi"][:]
        if "pha" in f:
            events["PHA"] = f["/pha"][:]
        parameters["sky_center"] = f["/sky_center"][:]
        parameters["dtheta"] = f["/dtheta"].value
        parameters["pix_center"] = f["/pix_center"][:]
        
        f.close()
        
        return cls(events, parameters)

    @classmethod
    def from_fits_file(cls, fitsfile):
        """
        Initialize an EventList from a FITS file with filename *fitsfile*.
        """
        hdulist = pyfits.open(fitsfile)

        tblhdu = hdulist["EVENTS"]
        
        events = {}
        parameters = {}
        
        parameters["ExposureTime"] = tblhdu.header["EXPOSURE"]
        parameters["Area"] = tblhdu.header["AREA"]
        parameters["Redshift"] = tblhdu.header["REDSHIFT"]
        parameters["AngularDiameterDistance"] = tblhdu.header["D_A"]
        if "RMF" in tblhdu.header:
            parameters["RMF"] = tblhdu["RMF"]
        if "ARF" in tblhdu.header:
            parameters["ARF"] = tblhdu["ARF"]
        if "CHANTYPE" in tblhdu.header:
            parameters["ChannelType"] = tblhdu["CHANTYPE"]
        if "TELESCOP" in tblhdu.header:
            parameters["Telescope"] = tblhdu["TELESCOP"]
        if "INSTRUME" in tblhdu.header:
            parameters["Instrument"] = tblhdu["INSTRUME"]
        parameters["sky_center"] = np.array([tblhdu["TCRVL2"],tblhdu["TCRVL3"]])
        parameters["pix_center"] = np.array([tblhdu["TCRVL2"],tblhdu["TCRVL3"]])
        parameters["dtheta"] = tblhdu["TCRVL3"]
        events["xpix"] = tblhdu.data.field("X")
        events["ypix"] = tblhdu.data.field("Y")
        events["eobs"] = tblhdu.data.field("ENERGY")/1000. # Convert to keV
        if "PI" in tblhdu.columns.names:
            events["PI"] = tblhdu.data.field("PI")
        if "PHA" in tblhdu.columns.names:
            events["PHA"] = tblhdu.data.field("PHA")
        
        return cls(events, parameters)

    @classmethod
    def join_events(cls, events1, events2):
        """
        Join two sets of events, *events1* and *events2*.
        """
        events = {}
        for item1, item2 in zip(events1.items(), events2.items()):
            k1, v1 = item1
            k2, v2 = item2
            events[k1] = np.concatenate([v1,v2])
        
        return cls(events, events1.parameters)

    def convolve_with_response(self, respfile):
        """
        Convolve the events with a RMF file *respfile*.
        """
        mylog.warning("This routine has not been tested to work with all RMFs. YMMV.")
        if not "ARF" in self.parameters:
            mylog.warning("Photons have not been processed with an"+
                          " auxiliary response file. Spectral fitting"+
                          " may be inaccurate.")

        mylog.info("Reading response matrix file (RMF): %s" % (respfile))
        
        hdulist = pyfits.open(respfile)

        tblhdu = hdulist["MATRIX"]
        n_de = len(tblhdu.data["ENERG_LO"])
        mylog.info("Number of Energy Bins: %d" % (n_de))
        de = tblhdu.data["ENERG_HI"] - tblhdu.data["ENERG_LO"]

        mylog.info("Energy limits: %g %g" % (min(tblhdu.data["ENERG_LO"]),
                                             max(tblhdu.data["ENERG_HI"])))

        if "ARF" in self.parameters:
            f = pyfits.open(self.parameters["ARF"])
            elo = f["SPECRESP"].data.field("ENERG_LO")
            ehi = f["SPECRESP"].data.field("ENERG_HI")
            f.close()
            try:
                assert_allclose(elo, tblhdu.data["ENERG_LO"], rtol=1.0e-6)
                assert_allclose(ehi, tblhdu.data["ENERG_HI"], rtol=1.0e-6)
            except AssertionError:
                mylog.warning("Energy binning does not match for "+
                              "ARF and RMF. This may make spectral "+
                              "fitting difficult.")
                                              
        tblhdu2 = hdulist["EBOUNDS"]
        n_ch = len(tblhdu2.data["CHANNEL"])
        mylog.info("Number of Channels: %d" % (n_ch))
        
        eidxs = np.argsort(self.events["eobs"])

        phEE = self.events["eobs"][eidxs]
        phXX = self.events["xpix"][eidxs]
        phYY = self.events["ypix"][eidxs]

        detectedChannels = []
        pindex = 0

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

        self.events["xpix"] = phXX
        self.events["ypix"] = phYY
        self.events["eobs"] = phEE
        self.events[tblhdu.header["CHANTYPE"]] = dchannel.astype(int)
        self.parameters["RMF"] = respfile
        self.parameters["ChannelType"] = tblhdu.header["CHANTYPE"]
        self.parameters["Telescope"] = tblhdu.header["TELESCOP"]
        self.parameters["Instrument"] = tblhdu.header["INSTRUME"]
        
    @parallel_root_only
    def write_fits_file(self, fitsfile, clobber=False):
        """
        Write events to a FITS binary table file with filename *fitsfile*.
        Set *clobber* to True if you need to overwrite a previous file.
        """
        
        cols = []

        col1 = pyfits.Column(name='ENERGY', format='E', unit='eV',
                             array=self.events["eobs"]*1000.)
        col2 = pyfits.Column(name='X', format='D', unit='pixel',
                             array=self.events["xpix"])
        col3 = pyfits.Column(name='Y', format='D', unit='pixel',
                             array=self.events["ypix"])

        cols = [col1, col2, col3]

        if self.parameters.has_key("ChannelType"):
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

        tbhdu.header.update("MTYPE1", "sky")
        tbhdu.header.update("MFORM1", "x,y")        
        tbhdu.header.update("MTYPE2", "EQPOS")
        tbhdu.header.update("MFORM2", "RA,DEC")
        tbhdu.header.update("TCTYP2", "RA---TAN")
        tbhdu.header.update("TCTYP3", "DEC--TAN")
        tbhdu.header.update("TCRVL2", self.parameters["sky_center"][0])
        tbhdu.header.update("TCRVL3", self.parameters["sky_center"][1])
        tbhdu.header.update("TCDLT2", -self.parameters["dtheta"])
        tbhdu.header.update("TCDLT3", self.parameters["dtheta"])
        tbhdu.header.update("TCRPX2", self.parameters["pix_center"][0])
        tbhdu.header.update("TCRPX3", self.parameters["pix_center"][1])
        tbhdu.header.update("TLMIN2", 0.5)
        tbhdu.header.update("TLMIN3", 0.5)
        tbhdu.header.update("TLMAX2", 2.*self.parameters["pix_center"][0]-0.5)
        tbhdu.header.update("TLMAX3", 2.*self.parameters["pix_center"][1]-0.5)
        tbhdu.header.update("EXPOSURE", self.parameters["ExposureTime"])
        tbhdu.header.update("AREA", self.parameters["Area"])
        tbhdu.header.update("D_A", self.parameters["AngularDiameterDistance"])
        tbhdu.header.update("REDSHIFT", self.parameters["Redshift"])
        tbhdu.header.update("HDUVERS", "1.1.0")
        tbhdu.header.update("RADECSYS", "FK5")
        tbhdu.header.update("EQUINOX", 2000.0)
        if "RMF" in self.parameters:
            tbhdu.header.update("RMF", self.parameters["RMF"])
        if "ARF" in self.parameters:
            tbhdu.header.update("ARF", self.parameters["ARF"])
        if "ChannelType" in self.parameters:
            tbhdu.header.update("CHANTYPE", self.parameters["ChannelType"])
        if "Telescope" in self.parameters:
            tbhdu.header.update("TELESCOP", self.parameters["Telescope"])
        if "Instrument" in self.parameters:
            tbhdu.header.update("INSTRUME", self.parameters["Instrument"])
            
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

        if isinstance(self.parameters["Area"], basestring):
             mylog.error("Writing SIMPUT files is only supported if you didn't convolve with an ARF.")
             raise TypeError
        
        if emin is None:
            emin = self.events["eobs"].min()*0.95
        if emax is None:
            emax = self.events["eobs"].max()*1.05

        idxs = np.logical_and(self.events["eobs"] >= emin, self.events["eobs"] <= emax)
        flux = erg_per_keV*np.sum(self.events["eobs"][idxs]) / \
               self.parameters["ExposureTime"]/self.parameters["Area"]
        
        col1 = pyfits.Column(name='ENERGY', format='E',
                             array=self["eobs"])
        col2 = pyfits.Column(name='DEC', format='D',
                             array=self["xsky"])
        col3 = pyfits.Column(name='RA', format='D',
                             array=self["ysky"])

        coldefs = pyfits.ColDefs([col1, col2, col3])

        tbhdu = pyfits.new_table(coldefs)
        tbhdu.update_ext_name("PHLIST")

        tbhdu.header.update("HDUCLASS", "HEASARC/SIMPUT")
        tbhdu.header.update("HDUCLAS1", "PHOTONS")
        tbhdu.header.update("HDUVERS", "1.1.0")
        tbhdu.header.update("EXTVER", 1)
        tbhdu.header.update("REFRA", 0.0)
        tbhdu.header.update("REFDEC", 0.0)
        tbhdu.header.update("TUNIT1", "keV")
        tbhdu.header.update("TUNIT2", "deg")
        tbhdu.header.update("TUNIT3", "deg")                

        phfile = prefix+"_phlist.fits"

        tbhdu.writeto(phfile, clobber=clobber)

        col1 = pyfits.Column(name='SRC_ID', format='J', array=np.array([1]).astype("int32"))
        col2 = pyfits.Column(name='RA', format='D', array=np.array([self.parameters["sky_center"][0]]))
        col3 = pyfits.Column(name='DEC', format='D', array=np.array([self.parameters["sky_center"][1]]))
        col4 = pyfits.Column(name='E_MIN', format='D', array=np.array([emin]))
        col5 = pyfits.Column(name='E_MAX', format='D', array=np.array([emax]))
        col6 = pyfits.Column(name='FLUX', format='D', array=np.array([flux]))
        col7 = pyfits.Column(name='SPECTRUM', format='80A', array=np.array([phfile+"[PHLIST,1]"]))
        col8 = pyfits.Column(name='IMAGE', format='80A', array=np.array([phfile+"[PHLIST,1]"]))
                        
        coldefs = pyfits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8])
        
        wrhdu = pyfits.new_table(coldefs)
        wrhdu.update_ext_name("SRC_CAT")
                                
        wrhdu.header.update("HDUCLASS", "HEASARC")
        wrhdu.header.update("HDUCLAS1", "SIMPUT")
        wrhdu.header.update("HDUCLAS2", "SRC_CAT")        
        wrhdu.header.update("HDUVERS", "1.1.0")
        wrhdu.header.update("RADECSYS", "FK5")
        wrhdu.header.update("EQUINOX", 2000.0)
        wrhdu.header.update("TUNIT2", "deg")
        wrhdu.header.update("TUNIT3", "deg")
        wrhdu.header.update("TUNIT4", "keV")
        wrhdu.header.update("TUNIT5", "keV")
        wrhdu.header.update("TUNIT6", "erg/s/cm**2")

        simputfile = prefix+"_simput.fits"
                
        wrhdu.writeto(simputfile, clobber=clobber)

    @parallel_root_only
    def write_h5_file(self, h5file):
        """
        Write an EventList to the HDF5 file given by *h5file*.
        """
        f = h5py.File(h5file, "w")

        f.create_dataset("/exp_time", data=self.parameters["ExposureTime"])
        f.create_dataset("/area", data=self.parameters["Area"])
        f.create_dataset("/redshift", data=self.parameters["Redshift"])
        f.create_dataset("/d_a", data=self.parameters["AngularDiameterDistance"])        
        if "ARF" in self.parameters:
            f.create_dataset("/arf", data=self.parameters["ARF"])
        if "RMF" in self.parameters:
            f.create_dataset("/rmf", data=self.parameters["RMF"])
        if "ChannelType" in self.parameters:
            f.create_dataset("/channel_type", data=self.parameters["ChannelType"])
        if "Telescope" in self.parameters:
            f.create_dataset("/telescope", data=self.parameters["Telescope"])
        if "Instrument" in self.parameters:
            f.create_dataset("/instrument", data=self.parameters["Instrument"])
                            
        f.create_dataset("/xpix", data=self.events["xpix"])
        f.create_dataset("/ypix", data=self.events["ypix"])
        f.create_dataset("/eobs", data=self.events["eobs"])
        if "PI" in self.events:
            f.create_dataset("/pi", data=self.events["PI"])                  
        if "PHA" in self.events:
            f.create_dataset("/pha", data=self.events["PHA"])                  
        f.create_dataset("/sky_center", data=self.parameters["sky_center"])
        f.create_dataset("/pix_center", data=self.parameters["pix_center"])
        f.create_dataset("/dtheta", data=self.parameters["dtheta"])

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
            mask_emin = self.events["eobs"] > emin
        if emax is None:
            mask_emax = np.ones((self.num_events), dtype='bool')
        else:
            mask_emax = self.events["eobs"] < emax

        mask = np.logical_and(mask_emin, mask_emax)

        nx = int(2*self.parameters["pix_center"][0]-1.)
        ny = int(2*self.parameters["pix_center"][1]-1.)
        
        xbins = np.linspace(0.5, float(nx)+0.5, nx+1, endpoint=True)
        ybins = np.linspace(0.5, float(ny)+0.5, ny+1, endpoint=True)

        H, xedges, yedges = np.histogram2d(self.events["xpix"][mask],
                                           self.events["ypix"][mask],
                                           bins=[xbins,ybins])
        
        hdu = pyfits.PrimaryHDU(H.T)
        
        hdu.header.update("MTYPE1", "EQPOS")
        hdu.header.update("MFORM1", "RA,DEC")
        hdu.header.update("CTYPE1", "RA---TAN")
        hdu.header.update("CTYPE2", "DEC--TAN")
        hdu.header.update("CRPIX1", 0.5*(nx+1))
        hdu.header.update("CRPIX2", 0.5*(nx+1))                
        hdu.header.update("CRVAL1", self.parameters["sky_center"][0])
        hdu.header.update("CRVAL2", self.parameters["sky_center"][1])
        hdu.header.update("CUNIT1", "deg")
        hdu.header.update("CUNIT2", "deg")
        hdu.header.update("CDELT1", -self.parameters["dtheta"])
        hdu.header.update("CDELT2", self.parameters["dtheta"])
        hdu.header.update("EXPOSURE", self.parameters["ExposureTime"])
        
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

        if energy_bins:
            spectype = "energy"
            espec = self.events["eobs"]
            range = (emin, emax)
            spec, ee = np.histogram(espec, bins=nchan, range=range)
            bins = 0.5*(ee[1:]+ee[:-1])
        else:
            if not self.parameters.has_key("ChannelType"):
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
        col4 = pyfits.Column(name='COUNT_RATE', format='1D', array=spec/self.parameters["ExposureTime"])
        
        coldefs = pyfits.ColDefs([col1, col2, col3, col4])
        
        tbhdu = pyfits.new_table(coldefs)
        tbhdu.update_ext_name("SPECTRUM")

        if not energy_bins:
            tbhdu.header.update("DETCHANS", spec.shape[0])
            tbhdu.header.update("TOTCTS", spec.sum())
            tbhdu.header.update("EXPOSURE", self.parameters["ExposureTime"])
            tbhdu.header.update("LIVETIME", self.parameters["ExposureTime"])        
            tbhdu.header.update("CONTENT", spectype)
            tbhdu.header.update("HDUCLASS", "OGIP")
            tbhdu.header.update("HDUCLAS1", "SPECTRUM")
            tbhdu.header.update("HDUCLAS2", "TOTAL")
            tbhdu.header.update("HDUCLAS3", "TYPE:I")
            tbhdu.header.update("HDUCLAS4", "COUNT")
            tbhdu.header.update("HDUVERS", "1.1.0")
            tbhdu.header.update("HDUVERS1", "1.1.0")
            tbhdu.header.update("CHANTYPE", spectype)
            tbhdu.header.update("BACKFILE", "none")
            tbhdu.header.update("CORRFILE", "none")
            tbhdu.header.update("POISSERR", True)
            if self.parameters.has_key("RMF"):
                 tbhdu.header.update("RESPFILE", self.parameters["RMF"])
            else:
                 tbhdu.header.update("RESPFILE", "none")
            if self.parameters.has_key("ARF"):
                tbhdu.header.update("ANCRFILE", self.parameters["ARF"])
            else:        
                tbhdu.header.update("ANCRFILE", "none")
            if self.parameters.has_key("Telescope"):
                tbhdu.header.update("TELESCOP", self.parameters["Telescope"])
            else:
                tbhdu.header.update("TELESCOP", "none")
            if self.parameters.has_key("Instrument"):
                tbhdu.header.update("INSTRUME", self.parameters["Instrument"])
            else:
                tbhdu.header.update("INSTRUME", "none")
            tbhdu.header.update("AREASCAL", 1.0)
            tbhdu.header.update("CORRSCAL", 0.0)
            tbhdu.header.update("BACKSCAL", 1.0)
                                
        hdulist = pyfits.HDUList([pyfits.PrimaryHDU(), tbhdu])
        
        hdulist.writeto(specfile, clobber=clobber)
