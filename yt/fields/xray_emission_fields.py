import os

import numpy as np

from yt.config import ytcfg
from yt.fields.derived_field import DerivedField
from yt.funcs import mylog, only_on_root, parse_h5_attr
from yt.units.yt_array import YTArray, YTQuantity
from yt.utilities.cosmology import Cosmology
from yt.utilities.exceptions import YTException, YTFieldNotFound
from yt.utilities.linear_interpolators import (
    BilinearFieldInterpolator,
    UnilinearFieldInterpolator,
)
from yt.utilities.on_demand_imports import _h5py as h5py

data_version = {"cloudy": 2, "apec": 3}

data_url = "http://yt-project.org/data"


def _get_data_file(table_type, data_dir=None):
    data_file = "%s_emissivity_v%d.h5" % (table_type, data_version[table_type])
    if data_dir is None:
        supp_data_dir = ytcfg.get("yt", "supp_data_dir")
        data_dir = supp_data_dir if os.path.exists(supp_data_dir) else "."
    data_path = os.path.join(data_dir, data_file)
    if not os.path.exists(data_path):
        msg = "Failed to find emissivity data file {}! Please download from {}".format(
            data_file,
            data_url,
        )
        mylog.error(msg)
        raise OSError(msg)
    return data_path


class EnergyBoundsException(YTException):
    def __init__(self, lower, upper):
        self.lower = lower
        self.upper = upper

    def __str__(self):
        return f"Energy bounds are {self.lower:e} to {self.upper:e} keV."


class ObsoleteDataException(YTException):
    def __init__(self, table_type):
        data_file = "%s_emissivity_v%d.h5" % (table_type, data_version[table_type])
        self.msg = "X-ray emissivity data is out of date.\n"
        self.msg += f"Download the latest data from {data_url}/{data_file}."

    def __str__(self):
        return self.msg


class XrayEmissivityIntegrator:
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
        only_on_root(mylog.info, "Loading emissivity data from %s", filename)
        in_file = h5py.File(filename, mode="r")
        if "info" in in_file.attrs:
            only_on_root(mylog.info, parse_h5_attr(in_file, "info"))
        if parse_h5_attr(in_file, "version") != data_version[table_type]:
            raise ObsoleteDataException(table_type)
        else:
            only_on_root(
                mylog.info,
                "X-ray '%s' emissivity data version: %s."
                % (table_type, parse_h5_attr(in_file, "version")),
            )

        self.log_T = in_file["log_T"][:]
        self.emissivity_primordial = in_file["emissivity_primordial"][:]
        if "log_nH" in in_file:
            self.log_nH = in_file["log_nH"][:]
        if use_metals:
            self.emissivity_metals = in_file["emissivity_metals"][:]
        self.ebin = YTArray(in_file["E"], "keV")
        in_file.close()
        self.dE = np.diff(self.ebin)
        self.emid = 0.5 * (self.ebin[1:] + self.ebin[:-1]).to("erg")
        self.redshift = redshift

    def get_interpolator(self, data_type, e_min, e_max, energy=True):
        data = getattr(self, f"emissivity_{data_type}")
        if not energy:
            data = data[..., :] / self.emid.v
        e_min = YTQuantity(e_min, "keV") * (1.0 + self.redshift)
        e_max = YTQuantity(e_max, "keV") * (1.0 + self.redshift)
        if (e_min - self.ebin[0]) / e_min < -1e-3 or (
            e_max - self.ebin[-1]
        ) / e_max > 1e-3:
            raise EnergyBoundsException(self.ebin[0], self.ebin[-1])
        e_is, e_ie = np.digitize([e_min, e_max], self.ebin)
        e_is = np.clip(e_is - 1, 0, self.ebin.size - 1)
        e_ie = np.clip(e_ie, 0, self.ebin.size - 1)

        my_dE = self.dE[e_is:e_ie].copy()
        # clip edge bins if the requested range is smaller
        my_dE[0] -= e_min - self.ebin[e_is]
        my_dE[-1] -= self.ebin[e_ie] - e_max

        interp_data = (data[..., e_is:e_ie] * my_dE).sum(axis=-1)
        if data.ndim == 2:
            emiss = UnilinearFieldInterpolator(
                np.log10(interp_data),
                [self.log_T[0], self.log_T[-1]],
                "log_T",
                truncate=True,
            )
        else:
            emiss = BilinearFieldInterpolator(
                np.log10(interp_data),
                [self.log_nH[0], self.log_nH[-1], self.log_T[0], self.log_T[-1]],
                ["log_nH", "log_T"],
                truncate=True,
            )

        return emiss


def add_xray_emissivity_field(
    ds,
    e_min,
    e_max,
    redshift=0.0,
    metallicity=("gas", "metallicity"),
    table_type="cloudy",
    data_dir=None,
    cosmology=None,
    dist=None,
    ftype="gas",
):
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
    dist : (value, unit) tuple or :class:`~yt.units.yt_array.YTQuantity`, optional
        The distance to the source, used for making intensity fields. You should
        only use this if your source is nearby (not cosmological). Default: None
    ftype : string, optional
        The field type to use when creating the fields, default "gas"

    This will create at least three fields:

    "xray_emissivity_{e_min}_{e_max}_keV" (erg s^-1 cm^-3)
    "xray_luminosity_{e_min}_{e_max}_keV" (erg s^-1)
    "xray_photon_emissivity_{e_min}_{e_max}_keV" (photons s^-1 cm^-3)

    and if a redshift or distance is specified it will create two others:

    "xray_intensity_{e_min}_{e_max}_keV" (erg s^-1 cm^-3 arcsec^-2)
    "xray_photon_intensity_{e_min}_{e_max}_keV" (photons s^-1 cm^-3 arcsec^-2)

    These latter two are really only useful when making projections.

    Examples
    --------

    >>> import yt
    >>> ds = yt.load("sloshing_nomag2_hdf5_plt_cnt_0100")
    >>> yt.add_xray_emissivity_field(ds, 0.5, 2)
    >>> p = yt.ProjectionPlot(
    ...     ds, "x", ("gas", "xray_emissivity_0.5_2_keV"), table_type="apec"
    ... )
    >>> p.save()
    """
    if not isinstance(metallicity, float) and metallicity is not None:
        try:
            metallicity = ds._get_field_info(*metallicity)
        except YTFieldNotFound as e:
            raise RuntimeError(
                f"Your dataset does not have a {metallicity} field! "
                + "Perhaps you should specify a constant metallicity instead?"
            ) from e

    if table_type == "cloudy":
        # Cloudy wants to scale by nH**2
        other_n = "H_nuclei_density"
    else:
        # APEC wants to scale by nH*ne
        other_n = "El_number_density"

    def _norm_field(field, data):
        return data[ftype, "H_nuclei_density"] * data[ftype, other_n]

    ds.add_field(
        (ftype, "norm_field"), _norm_field, units="cm**-6", sampling_type="local"
    )

    my_si = XrayEmissivityIntegrator(table_type, data_dir=data_dir, redshift=redshift)

    em_0 = my_si.get_interpolator("primordial", e_min, e_max)
    emp_0 = my_si.get_interpolator("primordial", e_min, e_max, energy=False)
    if metallicity is not None:
        em_Z = my_si.get_interpolator("metals", e_min, e_max)
        emp_Z = my_si.get_interpolator("metals", e_min, e_max, energy=False)

    def _emissivity_field(field, data):
        with np.errstate(all="ignore"):
            dd = {
                "log_nH": np.log10(data[ftype, "H_nuclei_density"]),
                "log_T": np.log10(data[ftype, "temperature"]),
            }

        my_emissivity = np.power(10, em_0(dd))
        if metallicity is not None:
            if isinstance(metallicity, DerivedField):
                my_Z = data[metallicity.name].to_value("Zsun")
            else:
                my_Z = metallicity
            my_emissivity += my_Z * np.power(10, em_Z(dd))

        my_emissivity[np.isnan(my_emissivity)] = 0

        return data[ftype, "norm_field"] * YTArray(my_emissivity, "erg*cm**3/s")

    emiss_name = (ftype, f"xray_emissivity_{e_min}_{e_max}_keV")
    ds.add_field(
        emiss_name,
        function=_emissivity_field,
        display_name=rf"\epsilon_{{X}} ({e_min}-{e_max} keV)",
        sampling_type="local",
        units="erg/cm**3/s",
    )

    def _luminosity_field(field, data):
        return data[emiss_name] * data[ftype, "mass"] / data[ftype, "density"]

    lum_name = (ftype, f"xray_luminosity_{e_min}_{e_max}_keV")
    ds.add_field(
        lum_name,
        function=_luminosity_field,
        display_name=rf"\rm{{L}}_{{X}} ({e_min}-{e_max} keV)",
        sampling_type="local",
        units="erg/s",
    )

    def _photon_emissivity_field(field, data):
        dd = {
            "log_nH": np.log10(data[ftype, "H_nuclei_density"]),
            "log_T": np.log10(data[ftype, "temperature"]),
        }

        my_emissivity = np.power(10, emp_0(dd))
        if metallicity is not None:
            if isinstance(metallicity, DerivedField):
                my_Z = data[metallicity.name].to_value("Zsun")
            else:
                my_Z = metallicity
            my_emissivity += my_Z * np.power(10, emp_Z(dd))

        return data[ftype, "norm_field"] * YTArray(my_emissivity, "photons*cm**3/s")

    phot_name = (ftype, f"xray_photon_emissivity_{e_min}_{e_max}_keV")
    ds.add_field(
        phot_name,
        function=_photon_emissivity_field,
        display_name=rf"\epsilon_{{X}} ({e_min}-{e_max} keV)",
        sampling_type="local",
        units="photons/cm**3/s",
    )

    fields = [emiss_name, lum_name, phot_name]

    if redshift > 0.0 or dist is not None:

        if dist is None:
            if cosmology is None:
                if hasattr(ds, "cosmology"):
                    cosmology = ds.cosmology
                else:
                    cosmology = Cosmology()
            D_L = cosmology.luminosity_distance(0.0, redshift)
            angular_scale = 1.0 / cosmology.angular_scale(0.0, redshift)
            dist_fac = ds.quan(
                1.0 / (4.0 * np.pi * D_L * D_L * angular_scale * angular_scale).v,
                "rad**-2",
            )
        else:
            redshift = 0.0  # Only for local sources!
            try:
                # normal behaviour, if dist is a YTQuantity
                dist = ds.quan(dist.value, dist.units)
            except AttributeError as e:
                try:
                    dist = ds.quan(*dist)
                except (RuntimeError, TypeError):
                    raise TypeError(
                        "dist should be a YTQuantity or a (value, unit) tuple!"
                    ) from e

            angular_scale = dist / ds.quan(1.0, "radian")
            dist_fac = ds.quan(
                1.0 / (4.0 * np.pi * dist * dist * angular_scale * angular_scale).v,
                "rad**-2",
            )

        ei_name = (ftype, f"xray_intensity_{e_min}_{e_max}_keV")

        def _intensity_field(field, data):
            I = dist_fac * data[emiss_name]
            return I.in_units("erg/cm**3/s/arcsec**2")

        ds.add_field(
            ei_name,
            function=_intensity_field,
            display_name=rf"I_{{X}} ({e_min}-{e_max} keV)",
            sampling_type="local",
            units="erg/cm**3/s/arcsec**2",
        )

        i_name = (ftype, f"xray_photon_intensity_{e_min}_{e_max}_keV")

        def _photon_intensity_field(field, data):
            I = (1.0 + redshift) * dist_fac * data[phot_name]
            return I.in_units("photons/cm**3/s/arcsec**2")

        ds.add_field(
            i_name,
            function=_photon_intensity_field,
            display_name=rf"I_{{X}} ({e_min}-{e_max} keV)",
            sampling_type="local",
            units="photons/cm**3/s/arcsec**2",
        )

        fields += [ei_name, i_name]

    for field in fields:
        mylog.info("Adding ('%s','%s') field.", field[0], field[1])

    return fields
