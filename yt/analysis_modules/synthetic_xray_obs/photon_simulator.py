import numpy as np
from yt.funcs import *
from yt.utilities.physical_constants import mp, clight, cm_per_kpc, \
     cm_per_mpc, cm_per_km, K_per_keV, erg_per_keV
from yt.utilities.cosmology import Cosmology
from yt.utilities.orientation import Orientation
from yt.utilities.parallel_tools.parallel_analysis_interface import \
     communication_system, parallel_root_only, get_mpi_type, parallel_capable
from yt.data_objects.api import add_field

import os
import h5py

try:
    import pyfits
except ImportError:
    try:
        import astropy.io.fits as pyfits
    except ImportError:
        mylog.warning("You don't have pyFITS installed. " + 
                      "Writing to and reading from FITS files won't be available.")
    
N_TBIN = 10000
TMIN = 8.08e-2
TMAX = 50.
FOUR_PI = 4.*np.pi

class XRayPhotonList(object):

    def __init__(self, photons=None, comm=None, cosmo=None, p_bins=None):
        if photons is None: photons = {}
        self.photons = photons
        self.comm = comm
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
    def from_file(cls, filename, cosmology=None):
        """
        Initialize a XRayPhotonList from an HDF5 file given by filename.
        """
        photons = {}

        comm = communication_system.communicators[-1]
                
        f = h5py.File(filename, "r")

        photons["FiducialExposureTime"] = f["/fid_exp_time"].value
        photons["FiducialArea"] = f["/fid_area"].value
        photons["FiducialRedshift"] = f["/fid_redshift"].value
        photons["Hubble0"] = f["/hubble"].value
        photons["OmegaMatter"] = f["/omega_matter"].value
        photons["OmegaLambda"] = f["/omega_lambda"].value                    
        photons["FiducialAngularDiameterDistance"] = f["/fid_d_a"].value

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

        if cosmology is None:
            cosmo = Cosmology(HubbleConstantNow=71., OmegaMatterNow=0.27,
                              OmegaLambdaNow=0.73)
        else:
            cosmo = cosmology
                                        
        return cls(photons=photons, comm=comm, cosmo=cosmo, p_bins=p_bins)

    @classmethod
    def from_scratch(cls, data_source, redshift, eff_A,
                     exp_time, emission_model, center="c",
                     X_H=0.75, Zmet=0.3, cosmology=None):
        """
        Initialize a XRayPhotonList from a data container. 
        """
        pf = data_source.pf

        comm = communication_system.communicators[-1]
                
        vol_scale = pf.units["cm"]**(-3)/np.prod(pf.domain_width)

        if cosmology is None:
            cosmo = Cosmology(HubbleConstantNow=71., OmegaMatterNow=0.27,
                              OmegaLambdaNow=0.73)
        else:
            cosmo = cosmology
        
        D_A = cosmo.AngularDiameterDistance(0.0,redshift)*cm_per_mpc
        cosmo_fac = 1.0/(FOUR_PI*D_A*D_A*(1.+redshift)**3)

        num_cells = data_source["Temperature"].shape[0]
        start_c = comm.rank*num_cells/comm.size
        end_c = (comm.rank+1)*num_cells/comm.size

        kT = data_source["Temperature"][start_c:end_c].copy()/K_per_keV
        vol = data_source["CellVolume"][start_c:end_c].copy()
        dx = data_source["dx"][start_c:end_c].copy()
        EM = (data_source["Density"][start_c:end_c].copy()/mp)**2
        EM *= 0.5*(1.+X_H)*X_H*vol
        
        data_source.clear_data()
        
        x = data_source["x"][start_c:end_c].copy()
        y = data_source["y"][start_c:end_c].copy()
        z = data_source["z"][start_c:end_c].copy()

        data_source.clear_data()
                
        vx = data_source["x-velocity"][start_c:end_c].copy()
        vy = data_source["y-velocity"][start_c:end_c].copy()
        vz = data_source["z-velocity"][start_c:end_c].copy()
        
        data_source.clear_data()
                
        idxs = np.argsort(kT)
        dshape = idxs.shape

        if center == "c":
            src_ctr = pf.domain_center
        elif center == "max":
            src_ctr = pf.h.find_max("Density")[-1]
        elif iterable(center):
            src_ctr = center
                    
        kT_bins = np.linspace(TMIN, max(kT[idxs][-1], TMAX), num=N_TBIN+1)
        dkT = kT_bins[1]-kT_bins[0]
        kT_idxs = np.digitize(kT[idxs], kT_bins)
        kT_idxs = np.minimum(np.maximum(1, kT_idxs), N_TBIN) - 1
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

        emission_model.prepare()
        energy = emission_model.ebins
        de = emission_model.de
        emid = 0.5*(energy[1:]+energy[:-1])
        
        cell_em = EM[idxs]*vol_scale
        cell_vol = vol[idxs]*vol_scale
        cell_emd = cell_em/cell_vol
        
        number_of_photons = np.zeros(dshape, dtype='uint64')
        energies = []
                                
        pbar = get_pbar("Generating Photons", dshape[0])

        for i, ikT in enumerate(kT_idxs):
            
            ncells = int(bcounts[i])
            ibegin = bcell[i]
            iend = ecell[i]
            kT = kT_bins[ikT] + 0.5*dkT
            
            em_sum = cell_em[ibegin:iend].sum()
            vol_sum = cell_vol[ibegin:iend].sum()
            em_avg = em_sum/vol_sum
            
            tot_norm = cosmo_fac*em_sum/vol_scale
            
            spec = emission_model.get_spectrum(kT, Zmet)
            spec *= tot_norm
            cumspec = np.cumsum(spec)
            counts = cumspec[:]/cumspec[-1]
            tot_ph = cumspec[-1]*eff_A*exp_time
            
            for icell in xrange(ibegin, iend):
                
                cell_norm = tot_ph * (cell_emd[icell]/em_avg) * (cell_vol[icell]/vol_sum)                    
                cell_Nph = int(cell_norm) + int(np.modf(cell_norm)[0] >= np.random.random())
                
                if cell_Nph > 0:
                    number_of_photons[icell] = cell_Nph                    
                    randvec = np.random.uniform(low=counts[0], high=counts[-1], size=cell_Nph)
                    randvec.sort()
                    eidxs = np.searchsorted(counts, randvec)-1
                    cell_e = emid[eidxs]+de*(randvec-counts[eidxs])/(counts[eidxs+1]-counts[eidxs])
                    energies.append(cell_e)
                            
                pbar.update(icell)
            
        pbar.finish()

        del cell_emd, cell_vol, cell_em

        active_cells = number_of_photons > 0
        idxs = idxs[active_cells]
        
        photons = {}
        photons["x"] = (x[idxs]-src_ctr[0])*pf.units["kpc"]
        photons["y"] = (y[idxs]-src_ctr[1])*pf.units["kpc"]
        photons["z"] = (z[idxs]-src_ctr[2])*pf.units["kpc"]
        photons["vx"] = vx[idxs]/cm_per_km
        photons["vy"] = vy[idxs]/cm_per_km
        photons["vz"] = vz[idxs]/cm_per_km
        photons["dx"] = dx[idxs]*pf.units["kpc"]
        photons["NumberOfPhotons"] = number_of_photons[active_cells]
        photons["Energy"] = np.concatenate(energies)
                
        photons["FiducialExposureTime"] = exp_time
        photons["FiducialArea"] = eff_A
        photons["FiducialRedshift"] = redshift
        photons["Hubble0"] = cosmo.HubbleConstantNow
        photons["OmegaMatter"] = cosmo.OmegaMatterNow
        photons["OmegaLambda"] = cosmo.OmegaLambdaNow
        photons["FiducialAngularDiameterDistance"] = D_A/cm_per_mpc

        p_bins = np.cumsum(photons["NumberOfPhotons"])
        p_bins = np.insert(p_bins, 0, [np.uint64(0)])
        
        return cls(photons=photons, comm=comm, cosmo=cosmo, p_bins=p_bins)

    def write_h5_file(self, photonfile):

        if parallel_capable:
            
            mpi_long = get_mpi_type("int64")
            mpi_double = get_mpi_type("float64")
        
            local_num_cells = len(self.photons["x"])
            sizes_c = self.comm.comm.gather(local_num_cells, root=0)
            
            local_num_photons = self.photons["NumberOfPhotons"].sum()
            sizes_p = self.comm.comm.gather(local_num_photons, root=0)
            
            if self.comm.rank == 0:
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
                                                
            self.comm.comm.Gatherv([self.photons["x"], local_num_cells, mpi_double],
                                   [x, (sizes_c, disps_c), mpi_double], root=0)
            self.comm.comm.Gatherv([self.photons["y"], local_num_cells, mpi_double],
                                   [y, (sizes_c, disps_c), mpi_double], root=0)
            self.comm.comm.Gatherv([self.photons["z"], local_num_cells, mpi_double],
                                   [z, (sizes_c, disps_c), mpi_double], root=0)
            self.comm.comm.Gatherv([self.photons["vx"], local_num_cells, mpi_double],
                                   [vx, (sizes_c, disps_c), mpi_double], root=0)
            self.comm.comm.Gatherv([self.photons["vy"], local_num_cells, mpi_double],
                                   [vy, (sizes_c, disps_c), mpi_double], root=0)
            self.comm.comm.Gatherv([self.photons["vz"], local_num_cells, mpi_double],
                                   [vz, (sizes_c, disps_c), mpi_double], root=0)
            self.comm.comm.Gatherv([self.photons["dx"], local_num_cells, mpi_double],
                                   [dx, (sizes_c, disps_c), mpi_double], root=0)
            self.comm.comm.Gatherv([self.photons["NumberOfPhotons"], local_num_cells, mpi_long],
                                   [n_ph, (sizes_c, disps_c), mpi_long], root=0)
            self.comm.comm.Gatherv([self.photons["Energy"], local_num_photons, mpi_double],
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
                                                
        if self.comm.rank == 0:
            
            f = h5py.File(photonfile, "w")

            # Scalars
       
            f.create_dataset("fid_area", data=self.photons["FiducialArea"])
            f.create_dataset("fid_exp_time", data=self.photons["FiducialExposureTime"])
            f.create_dataset("fid_redshift", data=self.photons["FiducialRedshift"])
            f.create_dataset("omega_matter", data=self.photons["OmegaMatter"])
            f.create_dataset("omega_lambda", data=self.photons["OmegaLambda"])
            f.create_dataset("hubble", data=self.photons["Hubble0"])
            f.create_dataset("fid_d_a", data=self.photons["FiducialAngularDiameterDistance"])
        
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

        self.comm.barrier()

    def project_photons(self, L, area_new=None, texp_new=None, 
                        redshift_new=None, absorb_model=None, psf_sigma=None):
        """
        Projects photons onto an image plane given a line of sight. 
        """

        dx = self.photons["dx"]
        
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
        
        if texp_new is None and area_new is None and redshift_new is None:
            my_n_obs = n_ph_tot
        else:
            if texp_new is None:
                Tratio = 1.
            else:
                Tratio = texp_new/self.photons["FiducialExposureTime"]
            if area_new is None:
                Aratio = 1.
            elif isinstance(area_new, basestring):
                mylog.info("Using energy-dependent effective area.")
                f = pyfits.open(area_new)
                elo = f[1].data.field("ENERG_LO")
                ehi = f[1].data.field("ENERG_HI")
                eff_area = f[1].data.field("SPECRESP")
                f.close()
                Aratio = eff_area.max()/self.photons["FiducialArea"]
            else:
                mylog.info("Using constant effective area.")
                Aratio = area_new/self.photons["FiducialArea"]
            if redshift_new is None:
                Zratio = 1.
                zobs = self.photons["FiducialRedshift"]
                D_A = self.photons["FiducialAngularDiameterDistance"]*1000.                    
            else:
                zobs = redshift_new
                fid_D_A = self.photons["FiducialAngularDiameterDistance"]*1000.
                D_A = self.cosmo.AngularDiameterDistance(0.0,zobs)*1000.
                Zratio = fid_D_A*fid_D_A*(1.+self.photons["FiducialRedshift"]**3) / \
                         (D_A*D_A*(1.+zobs)**3)
            fak = Aratio*Tratio*Zratio
            if fak > 1:
                raise ValueError("Spectrum scaling factor = %g, cannot be greater than unity." % (fak))
            my_n_obs = np.uint64(n_ph_tot*fak)

        n_obs_all = self.comm.mpi_allreduce(my_n_obs)
        if self.comm.rank == 0: mylog.info("Total number of photons to use: %d" % (n_obs_all))
        
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
            energy = absorb_model.ebins
            de = absorb_model.de
            emid = 0.5*(energy[1:]+energy[:-1])
            aspec = absorb_model.get_spectrum()
            eidxs = np.searchsorted(emid, eobs)-1
            dx = (eobs-emid[eidxs])/de
            absorb = aspec[eidxs]*(1.-dx) + aspec[eidxs+1]*dx
            randvec = aspec.max()*np.random.random(eobs.shape)
            not_abs = randvec < absorb

        if eff_area is None:
            detected = np.ones(eobs.shape, dtype='bool')
        else:
            mylog.info("Applying energy-dependent effective area.")
            earf = 0.5*(elo+ehi)
            de = earf[1]-earf[0]            
            eidxs = np.searchsorted(earf, eobs)-1
            dx = (eobs-earf[eidxs])/de
            earea = eff_area[eidxs]*(1.-dx) + eff_area[eidxs+1]*dx
            randvec = eff_area.max()*np.random.random(eobs.shape)
            detected = randvec < earea
        
        detected = np.logical_and(not_abs, detected)
                    
        events = {}

        events["xsky"] = np.rad2deg(xsky[detected]/D_A)*3600.
        events["ysky"] = np.rad2deg(ysky[detected]/D_A)*3600.
        events["eobs"] = eobs[detected]

        if psf_sigma is not None:
            events["xsky"] += np.random.normal(sigma=psf_sigma)
            events["ysky"] += np.random.normal(sigma=psf_sigma)

        events = self.comm.par_combine_object(events, datatype="dict", op="cat")
        
        num_events = len(events["xsky"])
            
        if self.comm.rank == 0: mylog.info("Total number of observed photons: %d" % (num_events))
                        
        if texp_new is None:
            events["ExposureTime"] = self.photons["FiducialExposureTime"]
        else:
            events["ExposureTime"] = texp_new
        if area_new is None:
            events["Area"] = self.photons["FiducialArea"]
        else:
            events["Area"] = area_new
        events["Hubble0"] = self.photons["Hubble0"]
        events["Redshift"] = zobs
        events["OmegaMatter"] = self.photons["OmegaMatter"]
        events["OmegaLambda"] = self.photons["OmegaLambda"] 
        events["AngularDiameterDistance"] = D_A/1000.
                
        return XRayEventList(events)

class XRayEventList(object) :

    def __init__(self, events = None) :

        if events is None : events = {}
        self.events = events
        self.num_events = events["xsky"].shape[0]

    def keys(self):
        return self.events.keys()

    def items(self):
        return self.events.items()

    def values(self):
        return self.events.values()
    
    def __getitem__(self,key):
        return self.events[key]
        
    @classmethod
    def from_h5_file(cls, h5file):
        """
        Initialize a XRayEventList from a HDF5 file with filename h5file.
        """
        events = {}
        
        f = h5py.File(h5file, "r")

        events["ExposureTime"] = f["/exp_time"].value
        events["Area"] = f["/area"].value
        events["Hubble0"] = f["/hubble"].value
        events["Redshift"] = f["/redshift"].value
        events["OmegaMatter"] = f["/omega_matter"].value
        events["OmegaLambda"] = f["/omega_lambda"].value
        events["AngularDiameterDistance"] = f["/d_a"].value
        
        events["xsky"] = f["/xsky"][:]
        events["ysky"] = f["/ysky"][:]
        events["eobs"] = f["/eobs"][:]
        
        f.close()
        
        return cls(events)

    @classmethod
    def from_fits_file(cls, fitsfile):
        """
        Initialize a XRayEventList from a FITS file with filename fitsfile.
        """
        hdulist = pyfits.open(fitsfile)

        tblhdu = hdulist[1]

        events = {}
        
        events["ExposureTime"] = tblhdu.header["EXPOSURE"]
        events["Area"] = tblhdu.header["AREA"]
        events["Hubble0"] = tblhdu.header["HUBBLE"]
        events["Redshift"] = tblhdu.header["REDSHIFT"]
        events["OmegaMatter"] = tblhdu.header["OMEGA_M"]
        events["OmegaLambda"] = tblhdu.header["OMEGA_L"]
        events["AngularDiameterDistance"] = tblhdu.header["D_A"]
        
        events["xsky"] = tblhdu.data.field("POS_X")
        events["ysky"] = tblhdu.data.field("POS_Y")
        events["eobs"] = tblhdu.data.field("ENERGY")
        
        return cls(events)

    @classmethod
    def join_events(cls, events1, events2):
        events = {}
        for item1, item2 in zip(events1.items(), events2.items()):
            k1, v1 = item1
            k2, v2 = item2
            if isinstance(v1, np.ndarray):
                events[k1] = np.concatenate([v1,v2])
            else:
                events[k1] = v1            
        return cls(events)
    
    def convolve_with_response(self, respfile):
        pass

    @parallel_root_only
    def write_fits_file(self, fitsfile, clobber=False):
        """
        Write events to a FITS binary table file.
        """
        col1 = pyfits.Column(name='ENERGY', format='E',
                             array=self.events["eobs"])
        col2 = pyfits.Column(name='XSKY', format='D',
                             array=self.events["xsky"])
        col3 = pyfits.Column(name='YSKY', format='D',
                             array=self.events["ysky"])
        
        coldefs = pyfits.ColDefs([col1, col2, col3])
        
        tbhdu = pyfits.new_table(coldefs)

        tbhdu.header.update("EXPOSURE", self.events["ExposureTime"])
        tbhdu.header.update("AREA", self.events["Area"])
        tbhdu.header.update("HUBBLE", self.events["Hubble0"])
        tbhdu.header.update("OMEGA_M", self.events["OmegaMatter"])
        tbhdu.header.update("OMEGA_L", self.events["OmegaLambda"])
        tbhdu.header.update("D_A", self.evenets["AngularDiameterDistance"])
                
        tbhdu.writeto(fitsfile, clobber=clobber)

    @parallel_root_only    
    def write_simput_file(self, prefix, clobber=False, e_min=None, e_max=None):
        """
        Write events to a SIMPUT file.
        """
        if e_min is None:
            e_min = 0.0
        if e_max is None:
            e_max = 100.0

        idxs = np.logical_and(self.events["eobs"] >= e_min, self.events["eobs"] <= e_max)
        flux = erg_per_keV*np.sum(self.events["eobs"][idxs])/self.events["ExposureTime"]/self.events["Area"]
        
        col1 = pyfits.Column(name='ENERGY', format='E',
                             array=self.events["eobs"])
        col2 = pyfits.Column(name='DEC', format='D',
                             array=self.events["ysky"]/3600.)
        col3 = pyfits.Column(name='RA', format='D',
                             array=self.events["xsky"]/3600.)

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
        col2 = pyfits.Column(name='RA', format='D', array=np.array([0.0]))
        col3 = pyfits.Column(name='DEC', format='D', array=np.array([0.0]))
        col4 = pyfits.Column(name='E_MIN', format='D', array=np.array([e_min]))
        col5 = pyfits.Column(name='E_MAX', format='D', array=np.array([e_max]))
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
        Write a XRayEventList to the HDF5 file given by h5file.
        """
        f = h5py.File(h5file, "w")

        f.create_dataset("/exp_time", data=self.events["ExposureTime"])
        f.create_dataset("/area", data=self.events["Area"])
        f.create_dataset("/redshift", data=self.events["Redshift"])
        f.create_dataset("/hubble", data=self.events["Hubble0"])
        f.create_dataset("/omega_matter", data=self.events["OmegaMatter"])
        f.create_dataset("/omega_lambda", data=self.events["OmegaLambda"])        
        f.create_dataset("/d_a", data=self.events["AngularDiameterDistance"])        
        f.create_dataset("/xsky", data=self.events["xsky"])
        f.create_dataset("/ysky", data=self.events["ysky"])
        f.create_dataset("/eobs", data=self.events["eobs"])
                        
        f.close()

    @parallel_root_only
    def write_fits_image(self, imagefile, width, nx, center,
                         clobber=False, gzip_file=False,
                         emin=None, emax=None):
        """
        Generate a image by binning X-ray counts and write it to a FITS file.
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

        dx_pixel = width/nx
        xmin = -0.5*width
        xmax = -xmin
        xbins = np.linspace(xmin, xmax, nx+1, endpoint=True)
        
        H, xedges, yedges = np.histogram2d(self.events["xsky"][mask],
                                           self.events["ysky"][mask],
                                           bins=[xbins,xbins])
        
        hdu = pyfits.PrimaryHDU(H.T[::-1,::])

        hdu.header.update("MTYPE1", "EQPOS")
        hdu.header.update("MFORM1", "RA,DEC")
        hdu.header.update("CTYPE1", "RA---TAN")
        hdu.header.update("CTYPE2", "DEC--TAN")
        hdu.header.update("CRPIX1", 0.5*(nx+1))
        hdu.header.update("CRPIX2", 0.5*(nx+1))                
        hdu.header.update("CRVAL1", center[0])
        hdu.header.update("CRVAL2", center[1])
        hdu.header.update("CUNIT1", "deg")
        hdu.header.update("CUNIT2", "deg")
        hdu.header.update("CDELT1", -dx_pixel/3600.)
        hdu.header.update("CDELT2", dx_pixel/3600.)
        hdu.header.update("EXPOSURE", self.events["ExposureTime"])
        
        hdu.writeto(imagefile, clobber=clobber)

        if (gzip_file):
            clob = ""
            if (clobber) : clob="-f"
            os.system("gzip "+clob+" %s.fits" % (prefix))
                                    
    @parallel_root_only
    def write_spectrum(self, specfile, emin, emax, nchan, clobber=False):
        
        spec, ee = np.histogram(self.events["eobs"], bins=nchan, range=(emin, emax))

        de = ee[1]-ee[0]
        emid = 0.5*(ee[1:]+ee[:-1])
        
        col1 = pyfits.Column(name='ENERGY', format='1E', array=emid)
        col2 = pyfits.Column(name='COUNTS', format='1E', array=spec)
        
        coldefs = pyfits.ColDefs([col1, col2])
        
        tbhdu = pyfits.new_table(coldefs)
        
        tbhdu.writeto(specfile, clobber=clobber)
