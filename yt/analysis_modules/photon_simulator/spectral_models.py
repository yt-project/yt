"""
Photon emission and absoprtion models for use with the
photon simulator.
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.on_demand_imports import _h5py as h5py
import numpy as np
import os

from yt.funcs import mylog
from yt.units.yt_array import YTArray, YTQuantity
from yt.utilities.on_demand_imports import _astropy
from yt.utilities.physical_constants import hcgs, clight
from yt.utilities.physical_ratios import erg_per_keV, amu_grams
from yt.analysis_modules.photon_simulator.utils import broaden_lines

hc = (hcgs*clight).in_units("keV*angstrom").v
cl = clight.v
K = 1.0/np.sqrt(2.*np.pi)

class SpectralModel(object):

    def __init__(self, emin, emax, nchan):
        self.emin = YTQuantity(emin, "keV")
        self.emax = YTQuantity(emax, "keV")
        self.nchan = nchan
        self.ebins = YTArray(np.linspace(self.emin, self.emax, nchan+1), "keV")
        self.de = np.diff(self.ebins)
        self.emid = 0.5*(self.ebins[1:]+self.ebins[:-1])

    def prepare_spectrum(self):
        pass

    def get_spectrum(self):
        pass

    def cleanup_spectrum(self):
        pass

class XSpecThermalModel(SpectralModel):
    r"""
    Initialize a thermal gas emission model from PyXspec.

    Parameters
    ----------
    model_name : string
        The name of the thermal emission model.
    emin : float
        The minimum energy for the spectral model.
    emax : float
        The maximum energy for the spectral model.
    nchan : integer
        The number of channels in the spectral model.
    settings : dictionary, optional
        A dictionary of key, value pairs (must both be strings)
        that can be used to set various options in XSPEC.

    Examples
    --------
    >>> mekal_model = XSpecThermalModel("mekal", 0.05, 50.0, 1000)
    """
    def __init__(self, model_name, emin, emax, nchan,
                 thermal_broad=False, settings=None):
        self.model_name = model_name
        self.thermal_broad = thermal_broad
        if settings is None: settings = {}
        self.settings = settings
        super(XSpecThermalModel, self).__init__(emin, emax, nchan)

    def prepare_spectrum(self, zobs):
        """
        Prepare the thermal model for execution given a redshift *zobs* for the spectrum.
        """
        import xspec
        xspec.Xset.chatter = 0
        xspec.AllModels.setEnergies("%f %f %d lin" %
                                    (self.emin.value, self.emax.value, self.nchan))
        self.model = xspec.Model(self.model_name)
        self.thermal_comp = getattr(self.model,self.model_name)
        if self.model_name == "bremss":
            self.norm = 3.02e-15
        else:
            self.norm = 1.0e-14
        self.thermal_comp.norm = 1.0
        self.thermal_comp.Redshift = zobs
        if self.thermal_broad:
            xspec.Xset.addModelString("APECTHERMAL","yes")
        for k,v in self.settings.items():
            xspec.Xset.addModelString(k,v)

    def get_spectrum(self, kT):
        """
        Get the thermal emission spectrum given a temperature *kT* in keV. 
        """
        self.thermal_comp.kT = kT
        self.thermal_comp.Abundanc = 0.0
        cosmic_spec = np.array(self.model.values(0))
        if self.model_name == "bremss":
            metal_spec = np.zeros(self.nchan)
        else:
            self.thermal_comp.Abundanc = 1.0
            metal_spec = np.array(self.model.values(0)) - cosmic_spec
        cosmic_spec *= self.norm
        metal_spec *= self.norm
        return YTArray(cosmic_spec, "cm**3/s"), YTArray(metal_spec, "cm**3/s")
    
    def cleanup_spectrum(self):
        del self.thermal_comp
        del self.model

class XSpecAbsorbModel(SpectralModel):
    r"""
    Initialize an absorption model from PyXspec.

    Parameters
    ----------
    model_name : string
        The name of the absorption model.
    nH : float
        The foreground column density *nH* in units of 10^22 cm^{-2}.
    emin : float, optional
        The minimum energy for the spectral model.
    emax : float, optional
        The maximum energy for the spectral model.
    nchan : integer, optional
        The number of channels in the spectral model.
    settings : dictionary, optional
        A dictionary of key, value pairs (must both be strings)
        that can be used to set various options in XSPEC.

    Examples
    --------
    >>> abs_model = XSpecAbsorbModel("wabs", 0.1)
    """
    def __init__(self, model_name, nH, emin=0.01, emax=50.0,
                 nchan=100000, settings=None):
        self.model_name = model_name
        self.nH = nH
        if settings is None: settings = {}
        self.settings = settings
        super(XSpecAbsorbModel, self).__init__(emin, emax, nchan)

    def prepare_spectrum(self):
        """
        Prepare the absorption model for execution given a redshift *zobs* for the spectrum.
        """
        import xspec
        xspec.Xset.chatter = 0
        xspec.AllModels.setEnergies("%f %f %d lin" %
                                    (self.emin.value, self.emax.value, self.nchan))
        self.model = xspec.Model(self.model_name+"*powerlaw")
        self.model.powerlaw.norm = self.nchan/(self.emax.value-self.emin.value)
        self.model.powerlaw.PhoIndex = 0.0
        for k,v in self.settings.items():
            xspec.Xset.addModelString(k,v)

    def get_spectrum(self):
        """
        Get the absorption spectrum.
        """
        m = getattr(self.model,self.model_name)
        m.nH = self.nH
        return np.array(self.model.values(0))

    def cleanup_spectrum(self):
        del self.model

class TableApecModel(SpectralModel):
    r"""
    Initialize a thermal gas emission model from the AtomDB APEC tables
    available at http://www.atomdb.org. This code borrows heavily from Python
    routines used to read the APEC tables developed by Adam Foster at the
    CfA (afoster@cfa.harvard.edu).

    Parameters
    ----------
    apec_root : string
        The directory root where the APEC model files are stored.
    emin : float
        The minimum energy for the spectral model.
    emax : float
        The maximum energy for the spectral model.
    nchan : integer
        The number of channels in the spectral model.
    apec_vers : string, optional
        The version identifier string for the APEC files, e.g.
        "2.0.2"
    thermal_broad : boolean, optional
        Whether to apply thermal broadening to spectral lines. Only should
        be used if you are attemping to simulate a high-spectral resolution
        detector.

    Examples
    --------
    >>> apec_model = TableApecModel("$SPECTRAL_DATA/spectral/", 0.05, 50.0,
    ...                             1000, apec_vers="3.0", thermal_broad=True)
    """
    def __init__(self, apec_root, emin, emax, nchan,
                 apec_vers="2.0.2", thermal_broad=False):
        self.apec_root = apec_root
        self.apec_prefix = "apec_v"+apec_vers
        self.cocofile = os.path.join(self.apec_root,
                                     self.apec_prefix+"_coco.fits")
        self.linefile = os.path.join(self.apec_root,
                                     self.apec_prefix+"_line.fits")
        super(TableApecModel, self).__init__(emin, emax, nchan)
        self.wvbins = hc/self.ebins[::-1].d
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

    def prepare_spectrum(self, zobs):
        """
        Prepare the thermal model for execution.
        """
        try:
            self.line_handle = _astropy.pyfits.open(self.linefile)
        except IOError:
            mylog.error("LINE file %s does not exist" % self.linefile)
            raise IOError("LINE file %s does not exist" % self.linefile)
        try:
            self.coco_handle = _astropy.pyfits.open(self.cocofile)
        except IOError:
            mylog.error("COCO file %s does not exist" % self.cocofile)
            raise IOError("COCO file %s does not exist" % self.cocofile)

        self.Tvals = self.line_handle[1].data.field("kT")
        self.dTvals = np.diff(self.Tvals)
        self.minlam = self.wvbins.min()
        self.maxlam = self.wvbins.max()
        self.scale_factor = 1.0/(1.+zobs)

    def _make_spectrum(self, kT, element, tindex):

        tmpspec = np.zeros(self.nchan)

        line_data = self.line_handle[tindex].data
        coco_data = self.coco_handle[tindex].data

        i = np.where((line_data.field('element') == element) &
                     (line_data.field('lambda') > self.minlam) &
                     (line_data.field('lambda') < self.maxlam))[0]

        E0 = hc/line_data.field('lambda')[i].astype("float64")*self.scale_factor
        amp = line_data.field('epsilon')[i].astype("float64")
        ebins = self.ebins.d
        de = self.de.d
        emid = self.emid.d
        if self.thermal_broad:
            sigma = E0*np.sqrt(2.*kT*erg_per_keV/(self.A[element]*amu_grams))/cl
            vec = broaden_lines(E0, sigma, amp, ebins)
        else:
            vec = np.histogram(E0, ebins, weights=amp)[0]
        tmpspec += vec

        ind = np.where((coco_data.field('Z') == element) &
                       (coco_data.field('rmJ') == 0))[0]
        if len(ind) == 0:
            return tmpspec
        else:
            ind = ind[0]

        n_cont = coco_data.field('N_Cont')[ind]
        e_cont = coco_data.field('E_Cont')[ind][:n_cont]
        continuum = coco_data.field('Continuum')[ind][:n_cont]

        tmpspec += np.interp(emid, e_cont*self.scale_factor, continuum)*de/self.scale_factor

        n_pseudo = coco_data.field('N_Pseudo')[ind]
        e_pseudo = coco_data.field('E_Pseudo')[ind][:n_pseudo]
        pseudo = coco_data.field('Pseudo')[ind][:n_pseudo]

        tmpspec += np.interp(emid, e_pseudo*self.scale_factor, pseudo)*de/self.scale_factor

        return tmpspec*self.scale_factor

    def get_spectrum(self, kT):
        """
        Get the thermal emission spectrum given a temperature *kT* in keV. 
        """
        cspec_l = np.zeros(self.nchan)
        mspec_l = np.zeros(self.nchan)
        cspec_r = np.zeros(self.nchan)
        mspec_r = np.zeros(self.nchan)
        tindex = np.searchsorted(self.Tvals, kT)-1
        if tindex >= self.Tvals.shape[0]-1 or tindex < 0:
            return YTArray(cspec_l, "cm**3/s"), YTArray(mspec_l, "cm**3/s")
        dT = (kT-self.Tvals[tindex])/self.dTvals[tindex]
        # First do H,He, and trace elements
        for elem in self.cosmic_elem:
            cspec_l += self._make_spectrum(kT, elem, tindex+2)
            cspec_r += self._make_spectrum(kT, elem, tindex+3)
        # Next do the metals
        for elem in self.metal_elem:
            mspec_l += self._make_spectrum(kT, elem, tindex+2)
            mspec_r += self._make_spectrum(kT, elem, tindex+3)
        cosmic_spec = YTArray(cspec_l*(1.-dT)+cspec_r*dT, "cm**3/s")
        metal_spec = YTArray(mspec_l*(1.-dT)+mspec_r*dT, "cm**3/s")
        return cosmic_spec, metal_spec

class TableAbsorbModel(SpectralModel):
    r"""
    Initialize an absorption model from a table stored in an HDF5 file.

    Parameters
    ----------
    filename : string
        The name of the table file.
    nH : float
        The foreground column density *nH* in units of 10^22 cm^{-2}.

    Examples
    --------
    >>> abs_model = XSpecAbsorbModel("abs_table.h5", 0.1)
    """
    def __init__(self, filename, nH):
        if not os.path.exists(filename):
            raise IOError("File does not exist: %s." % filename)
        self.filename = filename
        f = h5py.File(self.filename,"r")
        emin = f["energy"][:].min()
        emax = f["energy"][:].max()
        self.sigma = YTArray(f["cross_section"][:], "cm**2")
        nchan = self.sigma.shape[0]
        f.close()
        super(TableAbsorbModel, self).__init__(emin, emax, nchan)
        self.nH = YTQuantity(nH*1.0e22, "cm**-2")

    def prepare_spectrum(self):
        """
        Prepare the absorption model for execution.
        """
        pass

    def get_spectrum(self):
        """
        Get the absorption spectrum.
        """
        return np.exp(-self.sigma*self.nH)
