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

from yt.config import ytcfg
from yt.fields.derived_field import DerivedField
from yt.funcs import \
    mylog, \
    only_on_root, \
    issue_deprecation_warning, \
    parse_h5_attr
from yt.utilities.exceptions import YTFieldNotFound
from yt.utilities.exceptions import YTException
from yt.utilities.linear_interpolators import \
    UnilinearFieldInterpolator, BilinearFieldInterpolator
from yt.units.yt_array import YTArray, YTQuantity
from yt.utilities.cosmology import Cosmology

data_version = {"cloudy": 2,
                "apec": 2}

data_url = "http://yt-project.org/data"

def _get_data_file(table_type, data_dir=None):
    data_file = "%s_emissivity_v%d.h5" % (table_type, data_version[table_type])
    if data_dir is None:
        supp_data_dir = ytcfg.get("yt", "supp_data_dir")
        data_dir = supp_data_dir if os.path.exists(supp_data_dir) else "."
    data_path = os.path.join(data_dir, data_file)
    if not os.path.exists(data_path):
        msg = "Failed to find emissivity data file %s! " % data_file + \
            "Please download from http://yt-project.org/data!"
        mylog.error(msg)
        raise IOError(msg)
    return data_path

class EnergyBoundsException(YTException):
    def __init__(self, lower, upper):
        self.lower = lower
        self.upper = upper

    def __str__(self):
        return "Energy bounds are %e to %e keV." % \
          (self.lower, self.upper)

class ObsoleteDataException(YTException):
    def __init__(self, table_type):
        data_file = "%s_emissivity_v%d.h5" % (table_type, data_version[table_type])
        self.msg = "X-ray emissivity data is out of date.\n"
        self.msg += "Download the latest data from %s/%s." % (data_url, data_file)

    def __str__(self):
        return self.msg

class XrayEmissivityIntegrator(object):
    r"""Class for making X-ray emissivity fields. Uses hdf5 data tables
    generated from Cloudy and AtomDB/APEC.

    Initialize an XrayEmissivityIntegrator object.

    Parameters
    ----------
    table_type : string
        The type of data to use when computing the emissivity values. If "cloudy",
        a file called "cloudy_emissivity.h5" is used, for photoionized
        plasmas. If, "apec", a file called "apec_emissivity.h5" is used for 
        collisionally ionized plasmas. These files contain emissivity tables 
        for primordial elements and for metals at solar metallicity for the 
        energy range 0.1 to 100 keV.
    redshift : float, optional
        The cosmological redshift of the source of the field. Default: 0.0.
    data_dir : string, optional
        The location to look for the data table in. If not supplied, the file
        will be looked for in the location of the YT_DEST environment variable
        or in the current working directory.
    use_metals : boolean, optional
        If set to True, the emissivity will include contributions from metals.
        Default: True
    """
    def __init__(self, table_type, redshift=0.0, data_dir=None, use_metals=True):

        filename = _get_data_file(table_type, data_dir=data_dir)
        only_on_root(mylog.info, "Loading emissivity data from %s." % filename)
        in_file = h5py.File(filename, "r")
        if "info" in in_file.attrs:
            only_on_root(mylog.info, parse_h5_attr(in_file, "info"))
        if parse_h5_attr(in_file, "version") != data_version[table_type]:
            raise ObsoleteDataException(table_type)
        else:
            only_on_root(mylog.info, "X-ray '%s' emissivity data version: %s." % \
                         (table_type, parse_h5_attr(in_file, "version")))

        self.log_T = in_file["log_T"][:]
        self.emissivity_primordial = in_file["emissivity_primordial"][:]
        if "log_nH" in in_file:
            self.log_nH = in_file["log_nH"][:]
        if use_metals:
            self.emissivity_metals = in_file["emissivity_metals"][:]
        self.ebin = YTArray(in_file["E"], "keV")
        in_file.close()
        self.dE = np.diff(self.ebin)
        self.emid = 0.5*(self.ebin[1:]+self.ebin[:-1]).to("erg")
        self.redshift = redshift

    def get_interpolator(self, data_type, e_min, e_max, energy=True):
        data = getattr(self, "emissivity_%s" % data_type)
        if not energy:
            data = data[..., :] / self.emid.v
        e_min = YTQuantity(e_min, "keV")*(1.0+self.redshift)
        e_max = YTQuantity(e_max, "keV")*(1.0+self.redshift)
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
        if data.ndim == 2:
            emiss = UnilinearFieldInterpolator(np.log10(interp_data),
                                               [self.log_T[0],  self.log_T[-1]],
                                               "log_T", truncate=True)
        else:
            emiss = BilinearFieldInterpolator(np.log10(interp_data),
                                              [self.log_nH[0], self.log_nH[-1],
                                               self.log_T[0],  self.log_T[-1]],
                                              ["log_nH", "log_T"], truncate=True)

        return emiss

def add_xray_emissivity_field(ds, e_min, e_max, redshift=0.0,
                              metallicity=("gas", "metallicity"), 
                              table_type="cloudy", data_dir=None,
                              cosmology=None, **kwargs):
    r"""Create X-ray emissivity fields for a given energy range.

    Parameters
    ----------
    e_min : float
        The minimum energy in keV for the energy band.
    e_min : float
        The maximum energy in keV for the energy band.
    redshift : float, optional
        The cosmological redshift of the source of the field. Default: 0.0.
    metallicity : str or tuple of str or float, optional
        Either the name of a metallicity field or a single floating-point
        number specifying a spatially constant metallicity. Must be in
        solar units. If set to None, no metals will be assumed. Default: 
        ("gas", "metallicity")
    table_type : string, optional
        The type of emissivity table to be used when creating the fields. 
        Options are "cloudy" or "apec". Default: "cloudy"
    data_dir : string, optional
        The location to look for the data table in. If not supplied, the file
        will be looked for in the location of the YT_DEST environment variable
        or in the current working directory.
    cosmology : :class:`~yt.utilities.cosmology.Cosmology`, optional
        If set and redshift > 0.0, this cosmology will be used when computing the
        cosmological dependence of the emission fields. If not set, yt's default
        LCDM cosmology will be used.

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
    # The next several if constructs are for backwards-compatibility
    if "constant_metallicity" in kwargs:
        issue_deprecation_warning("The \"constant_metallicity\" parameter is deprecated. Set "
                                  "the \"metallicity\" parameter to a constant float value instead.")
        metallicity = kwargs["constant_metallicity"]
    if "with_metals" in kwargs:
        issue_deprecation_warning("The \"with_metals\" parameter is deprecated. Use the "
                                  "\"metallicity\" parameter to choose a constant or "
                                  "spatially varying metallicity.")
        if kwargs["with_metals"] and isinstance(metallicity, float):
            raise RuntimeError("\"with_metals=True\", but you specified a constant metallicity!")
        if not kwargs["with_metals"] and not isinstance(metallicity, float):
            raise RuntimeError("\"with_metals=False\", but you didn't specify a constant metallicity!")
    if not isinstance(metallicity, float) and metallicity is not None:
        try:
            metallicity = ds._get_field_info(*metallicity)
        except YTFieldNotFound:
            raise RuntimeError("Your dataset does not have a {} field! ".format(metallicity) +
                               "Perhaps you should specify a constant metallicity instead?")

    my_si = XrayEmissivityIntegrator(table_type, data_dir=data_dir, redshift=redshift)

    em_0 = my_si.get_interpolator("primordial", e_min, e_max)
    emp_0 = my_si.get_interpolator("primordial", e_min, e_max, energy=False)
    if metallicity is not None:
        em_Z = my_si.get_interpolator("metals", e_min, e_max)
        emp_Z = my_si.get_interpolator("metals", e_min, e_max, energy=False)

    def _emissivity_field(field, data):
        with np.errstate(all='ignore'):
            dd = {"log_nH": np.log10(data["gas", "H_nuclei_density"]),
                  "log_T": np.log10(data["gas", "temperature"])}

        my_emissivity = np.power(10, em_0(dd))
        if metallicity is not None:
            if isinstance(metallicity, DerivedField):
                my_Z = data[metallicity.name]
            else:
                my_Z = metallicity
            my_emissivity += my_Z * np.power(10, em_Z(dd))

        my_emissivity[np.isnan(my_emissivity)] = 0

        return data["gas","H_nuclei_density"]**2 * \
            YTArray(my_emissivity, "erg*cm**3/s")

    emiss_name = "xray_emissivity_%s_%s_keV" % (e_min, e_max)
    ds.add_field(("gas", emiss_name), function=_emissivity_field,
                 display_name=r"\epsilon_{X} (%s-%s keV)" % (e_min, e_max),
                 sampling_type="cell", units="erg/cm**3/s")

    def _luminosity_field(field, data):
        return data[emiss_name] * data["cell_volume"]

    lum_name = "xray_luminosity_%s_%s_keV" % (e_min, e_max)
    ds.add_field(("gas", lum_name), function=_luminosity_field,
                 display_name=r"\rm{L}_{X} (%s-%s keV)" % (e_min, e_max),
                 sampling_type="cell", units="erg/s")

    def _photon_emissivity_field(field, data):
        dd = {"log_nH": np.log10(data["gas", "H_nuclei_density"]),
              "log_T": np.log10(data["gas", "temperature"])}

        my_emissivity = np.power(10, emp_0(dd))
        if metallicity is not None:
            if isinstance(metallicity, DerivedField):
                my_Z = data[metallicity.name]
            else:
                my_Z = metallicity
            my_emissivity += my_Z * np.power(10, emp_Z(dd))

        return data["gas", "H_nuclei_density"]**2 * \
            YTArray(my_emissivity, "photons*cm**3/s")

    phot_name = "xray_photon_emissivity_%s_%s_keV" % (e_min, e_max)
    ds.add_field(("gas", phot_name), function=_photon_emissivity_field,
                 display_name=r"\epsilon_{X} (%s-%s keV)" % (e_min, e_max),
                 sampling_type="cell", units="photons/cm**3/s")

    fields = [emiss_name, lum_name, phot_name]

    if redshift > 0.0:

        if cosmology is None:
            if hasattr(ds, "cosmology"):
                cosmology = ds.cosmology
            else:
                cosmology = Cosmology()

        D_L = cosmology.luminosity_distance(0.0, redshift)
        angular_scale = 1.0/cosmology.angular_scale(0.0, redshift)
        dist_fac = 1.0/(4.0*np.pi*D_L*D_L*angular_scale*angular_scale)

        ei_name = "xray_intensity_%s_%s_keV" % (e_min, e_max)
        def _intensity_field(field, data):
            I = dist_fac*data[emiss_name]
            return I.in_units("erg/cm**3/s/arcsec**2")
        ds.add_field(("gas", ei_name), function=_intensity_field,
                     display_name=r"I_{X} (%s-%s keV)" % (e_min, e_max),
                     sampling_type="cell", units="erg/cm**3/s/arcsec**2")

        i_name = "xray_photon_intensity_%s_%s_keV" % (e_min, e_max)
        def _photon_intensity_field(field, data):
            I = (1.0+redshift)*dist_fac*data[phot_name]
            return I.in_units("photons/cm**3/s/arcsec**2")
        ds.add_field(("gas", i_name), function=_photon_intensity_field,
                     display_name=r"I_{X} (%s-%s keV)" % (e_min, e_max),
                     sampling_type="cell", units="photons/cm**3/s/arcsec**2")

        fields += [ei_name, i_name]

    [mylog.info("Adding %s field." % field) for field in fields]

    return fields
