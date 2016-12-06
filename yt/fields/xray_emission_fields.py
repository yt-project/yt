"""
Integrator classes to deal with interpolation and integration of input spectral
bins.  Currently only supports Cloudy and APEC-style data.



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

from yt.fields.derived_field import DerivedField
from yt.funcs import mylog, only_on_root, download_file
from yt.utilities.exceptions import YTFieldNotFound
from yt.utilities.exceptions import YTException
from yt.utilities.linear_interpolators import \
    UnilinearFieldInterpolator, BilinearFieldInterpolator
from yt.utilities.physical_constants import mp
from yt.units.yt_array import YTArray, YTQuantity
from yt.utilities.physical_ratios import \
    primordial_H_mass_fraction, erg_per_keV

def _get_data_file(table_type, data_dir=None):
    data_file = "%s_emissivity.h5" % table_type
    data_url = "http://yt-project.org/data"
    if data_dir is None:
        if "YT_DEST" in os.environ and \
            os.path.isdir(os.path.join(os.environ["YT_DEST"], "data")):
            # Try in default path
            data_dir = os.path.join(os.environ["YT_DEST"], "data")
        else:
            # Try in current working directory
            data_dir = "."
    data_path = os.path.join(data_dir, data_file)
    if not os.path.exists(data_path):
        mylog.info("Attempting to download supplementary data from %s to %s." %
                   (data_url, data_dir))
        fn = download_file(os.path.join(data_url, data_file), data_path)
        if fn != data_path:
            mylog.error("Failed to download supplementary data.")
            raise IOError("Failed to find emissivity data file %s!" % data_file)
    return data_path

class EnergyBoundsException(YTException):
    def __init__(self, lower, upper):
        self.lower = lower
        self.upper = upper

    def __str__(self):
        return "Energy bounds are %e to %e keV." % \
          (self.lower, self.upper)

class XrayEmissivityIntegrator(object):
    r"""Class for making X-ray emissivity fields. Uses hdf5 data tables
    generated from Cloudy and AtomDB/APEC.

    Initialize an XrayEmissivityIntegrator object.

    Parameters
    ----------
    table_type: string
        The type of data to use when computing the emissivity values. If "cloudy",
        a file called "cloudy_emissivity.h5" is used, for photoionized
        plasmas. If, "apec", a file called "apec_emissivity.h5" is used for 
        collisionally ionized plasmas. These files contain emissivity tables 
        for primordial elements and for metals at solar metallicity for the 
        energy range 0.1 to 100 keV.
    data_dir : string, optional
        The location to look for the data table in. If not supplied, the file
        will be looked for in the 
    """
    def __init__(self, table_type, data_dir=None):

        filename = _get_data_file(table_type, data_dir=data_dir)
        only_on_root(mylog.info, "Loading emissivity data from %s." % filename)
        in_file = h5py.File(filename, "r")
        if "info" in in_file.attrs:
            only_on_root(mylog.info, in_file.attrs["info"])
        only_on_root(mylog.info, "X-ray emissivity data version: %s." % \
                     in_file.attrs["version"])

        for field in ["emissivity_primordial", "emissivity_metals",
                      "log_nH", "log_T"]:
            if field in in_file:
                setattr(self, field, in_file[field][:])
        self.ebin = YTArray(in_file["ebin"], "keV")
        in_file.close()
        self.dE = np.diff(self.ebin)
        self.emid = 0.5*(self.ebin[1:]+self.ebin[:-1])

    def get_interpolator(self, data, e_min, e_max):
        e_min = YTQuantity(e_min, "keV")
        e_max = YTQuantity(e_max, "keV")
        if (e_min - self.ebin[0]) / e_min < -1e-3 or \
          (e_max - self.ebin[-1]) / e_max > 1e-3:
            raise EnergyBoundsException(self.ebin[0], self.ebin[-1])
        e_is, e_ie = np.digitize([e_min, e_max], self.ebin)
        e_is = np.clip(e_is - 1, 0, self.ebin.size - 1)
        e_ie = np.clip(e_ie, 0, self.ebin.size - 1)

        my_dE = self.dE[e_is: e_ie].copy()
        # clip edge bins if the requested range is smaller
        my_dE[0] -= e_min - self.ebin[e_is]
        my_dE[-1] -= self.ebin[e_ie] - e_max

        interp_data = (data[..., e_is:e_ie]*my_dE).sum(axis=-1)
        if len(data.shape) == 2:
            emiss = UnilinearFieldInterpolator(np.log10(interp_data),
                                               [self.log_T[0],  self.log_T[-1]],
                                               "log_T", truncate=True)
        else:
            emiss = BilinearFieldInterpolator(np.log10(interp_data),
                                              [self.log_nH[0], self.log_nH[-1],
                                               self.log_T[0],  self.log_T[-1]],
                                              ["log_nH", "log_T"], truncate=True)

        return emiss

def add_xray_emissivity_field(ds, e_min, e_max, 
                              metallicity=("gas", "metallicity"), 
                              table_type="cloudy"):
    r"""Create X-ray emissivity fields for a given energy range.

    Parameters
    ----------
    e_min : float
        The minimum energy in keV for the energy band.
    e_min : float
        The maximum energy in keV for the energy band.
    metallicity : field or float, optional
        Either the name of a metallicity field or a single floating-point
        number specifying a spatially constant metallicity. Must be in
        solar units. If set to None, no metals will be assumed. Default: 
        ("gas", "metallicity")
    table_type : string, optional
        The type of emissivity table to be used when creating the fields. 
        Options are "cloudy" or "apec". Default: "cloudy"

    This will create three fields:

    "xray_emissivity_{e_min}_{e_max}_keV" (erg s^-1 cm^-3)
    "xray_luminosity_{e_min}_{e_max}_keV" (erg s^-1)
    "xray_photon_emissivity_{e_min}_{e_max}_keV" (photons s^-1 cm^-3)

    Examples
    --------

    >>> import yt
    >>> ds = yt.load("sloshing_nomag2_hdf5_plt_cnt_0100")
    >>> yt.add_xray_emissivity_field(ds, 0.5, 2)
    >>> p = yt.ProjectionPlot(ds, 'x', "xray_emissivity_0.5_2_keV")
    >>> p.save()
    """
    if not isinstance(metallicity, float) and metallicity is not None:
        try:
            metallicity = ds._get_field_info(metallicity)
        except YTFieldNotFound:
            raise RuntimeError("Your dataset does not have a %s field! " % metallicity +
                               "Perhaps you should specify a constant metallicity?")

    my_si = XrayEmissivityIntegrator(table_type)

    em_0 = my_si.get_interpolator(my_si.emissivity_primordial, e_min, e_max)
    if metallicity is not None:
        em_Z = my_si.get_interpolator(my_si.emissivity_metals, e_min, e_max)

    energy_erg = my_si.emid.v * erg_per_keV
    emp_0 = my_si.get_interpolator((my_si.emissivity_primordial[..., :] / energy_erg),
                                   e_min, e_max)
    if metallicity is not None:
        emp_Z = my_si.get_interpolator((my_si.emissivity_metals[..., :] / energy_erg),
                                       e_min, e_max)

    try:
        ds._get_field_info("H_number_density")
    except YTFieldNotFound:
        mylog.warning("Could not find a field for \"H_number_density\". "
                      "Assuming primordial H mass fraction.")
        def _nh(field, data):
            return primordial_H_mass_fraction*data["gas", "density"]/mp
        ds.add_field(("gas", "H_number_density"), function=_nh, units="cm**-3")

    def _emissivity_field(field, data):
        dd = {"log_nH": np.log10(data["gas", "H_number_density"]),
              "log_T": np.log10(data["gas", "temperature"])}

        my_emissivity = np.power(10, em_0(dd))
        if metallicity is not None:
            if isinstance(metallicity, DerivedField):
                my_Z = data[metallicity]
            else:
                my_Z = metallicity
            my_emissivity += my_Z * np.power(10, em_Z(dd))

        return data["gas","H_number_density"]**2 * \
            YTArray(my_emissivity, "erg*cm**3/s")

    emiss_name = "xray_emissivity_%s_%s_keV" % (e_min, e_max)
    ds.add_field(("gas", emiss_name), function=_emissivity_field,
                 display_name=r"\epsilon_{X} (%s-%s keV)" % (e_min, e_max),
                 units="erg/cm**3/s")

    def _luminosity_field(field, data):
        return data[emiss_name] * data["cell_volume"]

    lum_name = "xray_luminosity_%s_%s_keV" % (e_min, e_max)
    ds.add_field(("gas", lum_name), function=_luminosity_field,
                 display_name=r"\rm{L}_{X} (%s-%s keV)" % (e_min, e_max),
                 units="erg/s")

    def _photon_emissivity_field(field, data):
        dd = {"log_nH": np.log10(data["gas", "H_number_density"]),
              "log_T": np.log10(data["gas", "temperature"])}

        my_emissivity = np.power(10, emp_0(dd))
        if metallicity is not None:
            if isinstance(metallicity, DerivedField):
                my_Z = data[metallicity]
            else:
                my_Z = metallicity
            my_emissivity += my_Z * np.power(10, emp_Z(dd))

        return data["gas", "H_number_density"]**2 * \
            YTArray(my_emissivity, "photons*cm**3/s")

    phot_name = "xray_photon_emissivity_%s_%s_keV" % (e_min, e_max)
    ds.add_field(("gas", phot_name), function=_photon_emissivity_field,
                 display_name=r"\epsilon_{X} (%s-%s keV)" % (e_min, e_max),
                 units="photons/cm**3/s")

    return emiss_name, lum_name, phot_name
