import numpy as np
from yt.funcs import *
import h5py

try:
    import xspec
except ImportError:
    mylog.warning("You don't have PyXSpec installed. Some models won't be available.")
                
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

    def __init__(self, filename):
        if not os.path.exists(filename):
            raise IOError("File does not exist: %s." % filename)
        self.filename = filename
        f = h5py.File(self.filename,"r")
        self.T_vals = f["kT"][:]
        self.Z_vals = f["Zmet"][:]
        emin        = f["emin"].value
        emax        = f["emax"].value
        self.spec_table = f["spectrum"][:,:,:]
        nchan = self.spec_table.shape[-1]
        f.close()
        self.dT = self.T_vals[1]-self.T_vals[0]
        self.dZ = self.Z_vals[1]-self.Z_vals[0]
        PhotonModel.__init__(self, emin, emax, nchan)
        
    def prepare(self):
        pass

    def get_spectrum(self, kT, Zmet):
        iz = np.searchsorted(self.Z_vals, Zmet)-1
        dz = (Zmet-self.Z_vals[iz])/self.dZ
        it = np.searchsorted(self.T_vals, kT)-1
        dt = (kT-self.T_vals[it])/self.dT
        spec = self.spec_table[it+1,iz+1,:]*dt*dz + \
               self.spec_table[it,iz,:]*(1.-dt)*(1.-dz) + \
               self.spec_table[it,iz+1,:]*(1.-dt)*dz + \
               self.spec_table[it+1,iz,:]*dt*(1.-dz)
        return spec
    
class TableAbsorbModel(PhotonModel):

    def __init__(self, filename):
        if not os.path.exists(filename):
            raise IOError("File does not exist: %s." % filename)
        self.filename = filename
        f = h5py.File(self.filename,"r")
        emin = f["emin"].value
        emax = f["emax"].value
        self.abs = f["spectrum"][:]
        nchan = self.abs.shape[0]
        f.close()
        PhotonModel.__init__(self, emin, emax, nchan)
                                                                    
    def prepare(self):
        pass

    def get_spectrum(self):
        return self.abs ** nH
    
