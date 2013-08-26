import numpy as np
import os
from yt.funcs import *
import h5py
try:
    import pyfits
except:
    try:
        import astropy.io.fits as pyfits
    except:
        mylog.error("Your mom doesn't have pyFITS")
try:
    import xspec
except ImportError:
    mylog.warning("You don't have PyXSpec installed. Some models won't be available.")
from scipy.integrate import cumtrapz
from scipy import stats
from yt.utilities.physical_constants import hcgs, clight, erg_per_keV, amu_cgs

hc = 1.0e8*hcgs*clight/erg_per_keV

class PhotonModel(object):

    def __init__(self, emin, emax, nchan):
        self.emin = emin
        self.emax = emax
        self.nchan = nchan
        self.ebins = np.linspace(emin, emax, nchan+1)
        self.de = self.ebins[1]-self.ebins[0]
        
    def prepare(self):
        pass
    
    def get_spectrum(self):
        pass
                                                        
class XSpecThermalModel(PhotonModel):

    def __init__(self, model_name, emin, emax, nchan):
        self.model_name = model_name
        PhotonModel.__init__(self, emin, emax, nchan)
        
    def prepare(self):
        xspec.Xset.chatter = 0
        xspec.AllModels.setEnergies("%f %f %d lin" %
                                    (self.emin, self.emax, self.nchan))
        self.model = xspec.Model(self.model_name)
        
    def get_spectrum(self, kT, Zmet):
        m = getattr(self.model,self.model_name)
        m.kT = kT
        m.Abundanc = Zmet
        m.norm = 1.0
        m.Redshift = 0.0
        return 1.0e-14*np.array(self.model.values(0))
    
class XSpecAbsorbModel(PhotonModel):

    def __init__(self, model_name, nH, emin=0.01, emax=50.0, nchan=100000):
        self.model_name = model_name
        self.nH = nH
        PhotonModel.__init__(self, emin, emax, nchan)
        
    def prepare(self):
        xspec.Xset.chatter = 0
        xspec.AllModels.setEnergies("%f %f %d lin" %
                                    (self.emin, self.emax, self.nchan))
        self.model = xspec.Model(self.model_name+"*powerlaw")
        self.model.powerlaw.norm = self.nchan/(self.emax-self.emin)
        self.model.powerlaw.PhoIndex = 0.0

    def get_spectrum(self):
        m = getattr(self.model,self.model_name)
        m.nH = self.nH
        return np.array(self.model.values(0))

class TableApecModel(PhotonModel):

    def __init__(self, apec_root, emin, emax, nchan,
                 apec_vers="2.0.2", thermal_broad=False):
        self.apec_root = apec_root
        self.apec_prefix = "apec_v"+apec_vers
        self.cocofile = os.path.join(self.apec_root,
                                     self.apec_prefix+"_coco.fits")
        self.linefile = os.path.join(self.apec_root,
                                     self.apec_prefix+"_line.fits")
        PhotonModel.__init__(self, emin, emax, nchan)
        self.wvbins = hc/self.ebins[::-1]
        # H, He, and trace elements
        self.cosmic_elem = [1,2,3,4,5,9,11,15,17,19,21,22,23,24,25,27,29,30]
        # Non-trace metals
        self.metal_elem = [6,7,8,10,12,13,14,16,18,20,26,28]
        self.thermal_broad = thermal_broad
        self.A = np.array([0.0,1.00794,4.00262,6.941,9.012182,10.811,
                           12.0107,14.0067,15.9994,18.9984,20.1797,
                           22.9898,24.3050,26.9815,28.0855,30.9738,
                           32.0650,35.4530,39.9480,39.0983,40.0780,
                           44.9559,47.8670,50.9415,51.9961,54.9380,
                           55.8450,58.9332,58.6934,63.5460,65.3800])
        
    def prepare(self):
        try:
            self.line_handle = pyfits.open(self.linefile)
        except IOError:
            mylog.error("LINE file %s does not exist" % (self.linefile))
        try:
            self.coco_handle = pyfits.open(self.cocofile)
        except IOError:
            mylog.error("COCO file %s does not exist" % (self.cocofile))
        self.Tvals = self.line_handle[1].data.field("kT")
        self.dTvals = np.diff(self.Tvals)
        self.minlam = self.wvbins.min()
        self.maxlam = self.wvbins.max()

    def expand_E_grid(self, Eedges, Ncont, Econt, cont):
        # linearly interpolate onto wvedges
        # All energies should be in keV
        newx = np.array(Econt[:Ncont])
        newx = np.append(newx,Eedges)
        newxargs = newx.argsort()
        newx.sort()
        newy = np.interp(newx, Econt, cont)
        # simple integration
        ct = cumtrapz(newy,newx)
        ret = np.zeros((self.nchan))
        iinew = np.where(newxargs==Ncont)[0][0]
        for i in range(self.nchan):
            ind = i + 1 + Ncont
            iiold = iinew
            iinew = np.where(newxargs==ind)[0][0]
            ret[i] = ct[iinew-1]-ct[iiold-1]
        return ret

    def make_spectrum(self, element, tindex):
        
        tmpspec = np.zeros((self.nchan))
        
        i = np.where((self.line_handle[tindex].data.field('element')==element) &
                     (self.line_handle[tindex].data.field('lambda') > self.minlam) &
                     (self.line_handle[tindex].data.field('lambda') < self.maxlam))[0]

        vec = np.zeros((self.nchan))
        E0 = hc/self.line_handle[tindex].data.field('lambda')[i]
        amp = self.line_handle[tindex].data.field('epsilon')[i]
        if self.thermal_broad:
            vec = np.zeros((self.nchan))
            sigma = E0*np.sqrt(self.Tvals[tindex]*erg_per_keV/(self.A[element]*amu_cgs))/clight
            for E, sig, a in zip(E0, sigma, amp):
                cdf = stats.norm(E,sig).cdf(self.ebins)
                vec += np.diff(cdf)*a
        else:
            ie = np.searchsorted(self.ebins, E0, side='right')-1
            for i,a in zip(ie,amp): vec[i] += a
        tmpspec += vec

        ind = np.where((self.coco_handle[tindex].data.field('Z')==element) &
                       (self.coco_handle[tindex].data.field('rmJ')==0))[0]
        if len(ind)==0:
            return tmpspec
        else:
            ind=ind[0]
                                                    
        n_cont=self.coco_handle[tindex].data.field('N_Cont')[ind]
        e_cont=self.coco_handle[tindex].data.field('E_Cont')[ind][:n_cont]
        continuum = self.coco_handle[tindex].data.field('Continuum')[ind][:n_cont]
        
        tmpspec += self.expand_E_grid(self.ebins, n_cont, e_cont, continuum)
        
        n_pseudo=self.coco_handle[tindex].data.field('N_Pseudo')[ind]
        e_pseudo=self.coco_handle[tindex].data.field('E_Pseudo')[ind][:n_pseudo]
        pseudo = self.coco_handle[tindex].data.field('Pseudo')[ind][:n_pseudo]
        
        tmpspec += self.expand_E_grid(self.ebins, n_pseudo, e_pseudo, pseudo)
                                                        
        return tmpspec

    def get_spectrum(self, kT, Zmet):
        cspec_l = np.zeros((self.nchan))
        mspec_l = np.zeros((self.nchan))
        cspec_r = np.zeros((self.nchan))
        mspec_r = np.zeros((self.nchan))
        tindex = np.searchsorted(self.Tvals, kT)-1
        dT = (kT-self.Tvals[tindex])/self.dTvals[tindex]
        # First do H,He, and trace elements
        for elem in self.cosmic_elem:
            cspec_l += self.make_spectrum(elem, tindex+2)
            cspec_r += self.make_spectrum(elem, tindex+3)            
        # Next do the metals
        for elem in self.metal_elem:
            mspec_l += self.make_spectrum(elem, tindex+2)
            mspec_r += self.make_spectrum(elem, tindex+3)
        cosmic_spec = cspec_l*(1.-dT)+cspec_r*dT
        metal_spec = mspec_l*(1.-dT)+mspec_r*dT        
        return cosmic_spec+Zmet*metal_spec

class TableAbsorbModel(PhotonModel):

    def __init__(self, filename, nH):
        if not os.path.exists(filename):
            raise IOError("File does not exist: %s." % filename)
        self.filename = filename
        f = h5py.File(self.filename,"r")
        emin = f["energ_lo"][:].min()
        emax = f["energ_hi"][:].max()
        self.sigma = f["cross_section"][:]
        nchan = self.sigma.shape[0]
        f.close()
        PhotonModel.__init__(self, emin, emax, nchan)
        self.nH = nH
        
    def prepare(self):
        pass
        
    def get_spectrum(self):
        return np.exp(-self.sigma*self.nH)
