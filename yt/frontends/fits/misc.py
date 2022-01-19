import base64
import os
from io import BytesIO

import numpy as np

from yt.fields.derived_field import ValidateSpatial
from yt.units.yt_array import YTArray, YTQuantity
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.on_demand_imports import _astropy


def _make_counts(emin, emax):
    def _counts(field, data):
        e = data[("all", "event_energy")].in_units("keV")
        mask = np.logical_and(e >= emin, e < emax)
        x = data[("all", "event_x")][mask]
        y = data[("all", "event_y")][mask]
        z = np.ones(x.shape)
        pos = np.array([x, y, z]).transpose()
        img = data.deposit(pos, method="count")
        if data.has_field_parameter("sigma"):
            sigma = data.get_field_parameter("sigma")
        else:
            sigma = None
        if sigma is not None and sigma > 0.0:
            kern = _astropy.conv.Gaussian2DKernel(x_stddev=sigma)
            img[:, :, 0] = _astropy.conv.convolve(img[:, :, 0], kern)
        return data.ds.arr(img, "counts/pixel")

    return _counts


def setup_counts_fields(ds, ebounds, ftype="gas"):
    r"""
    Create deposited image fields from X-ray count data in energy bands.

    Parameters
    ----------
    ds : ~yt.data_objects.static_output.Dataset
        The FITS events file dataset to add the counts fields to.
    ebounds : list of tuples
        A list of tuples, one for each field, with (emin, emax) as the
        energy bounds for the image.
    ftype : string, optional
        The field type of the resulting field. Defaults to "gas".

    Examples
    --------
    >>> ds = yt.load("evt.fits")
    >>> ebounds = [(0.1, 2.0), (2.0, 3.0)]
    >>> setup_counts_fields(ds, ebounds)
    """
    for (emin, emax) in ebounds:
        cfunc = _make_counts(emin, emax)
        fname = f"counts_{emin}-{emax}"
        mylog.info("Creating counts field %s.", fname)
        ds.add_field(
            (ftype, fname),
            sampling_type="cell",
            function=cfunc,
            units="counts/pixel",
            validators=[ValidateSpatial()],
            display_name=f"Counts ({emin}-{emax} keV)",
        )


def create_spectral_slabs(filename, slab_centers, slab_width, **kwargs):
    r"""
    Given a dictionary of spectral slab centers and a width in
    spectral units, extract data from a spectral cube at these slab
    centers and return a `FITSDataset` instance containing the different
    slabs as separate yt fields. Useful for extracting individual
    lines from a spectral cube and separating them out as different fields.

    Requires the SpectralCube (https://spectral-cube.readthedocs.io/en/latest/)
    library.

    All keyword arguments will be passed on to the `FITSDataset` constructor.

    Parameters
    ----------
    filename : string
        The spectral cube FITS file to extract the data from.
    slab_centers : dict of (float, string) tuples or YTQuantities
        The centers of the slabs, where the keys are the names
        of the new fields and the values are (float, string) tuples or
        YTQuantities, specifying a value for each center and its unit.
    slab_width : YTQuantity or (float, string) tuple
        The width of the slab along the spectral axis.

    Examples
    --------
    >>> slab_centers = {
    ...     "13CN": (218.03117, "GHz"),
    ...     "CH3CH2CHO": (218.284256, "GHz"),
    ...     "CH3NH2": (218.40956, "GHz"),
    ... }
    >>> slab_width = (0.05, "GHz")
    >>> ds = create_spectral_slabs(
    ...     "intensity_cube.fits", slab_centers, slab_width, nan_mask=0.0
    ... )
    """
    from spectral_cube import SpectralCube

    from yt.frontends.fits.api import FITSDataset
    from yt.visualization.fits_image import FITSImageData

    cube = SpectralCube.read(filename)
    if not isinstance(slab_width, YTQuantity):
        slab_width = YTQuantity(slab_width[0], slab_width[1])
    slab_data = {}
    field_units = cube.header.get("bunit", "dimensionless")
    for k, v in slab_centers.items():
        if not isinstance(v, YTQuantity):
            slab_center = YTQuantity(v[0], v[1])
        else:
            slab_center = v
        mylog.info("Adding slab field %s at %g %s", k, slab_center.v, slab_center.units)
        slab_lo = (slab_center - 0.5 * slab_width).to_astropy()
        slab_hi = (slab_center + 0.5 * slab_width).to_astropy()
        subcube = cube.spectral_slab(slab_lo, slab_hi)
        slab_data[k] = YTArray(subcube.filled_data[:, :, :], field_units)
    width = subcube.header["naxis3"] * cube.header["cdelt3"]
    w = subcube.wcs.copy()
    w.wcs.crpix[-1] = 0.5
    w.wcs.crval[-1] = -0.5 * width
    fid = FITSImageData(slab_data, wcs=w)
    for hdu in fid:
        hdu.header.pop("RESTFREQ", None)
        hdu.header.pop("RESTFRQ", None)
    ds = FITSDataset(fid, **kwargs)
    return ds


def ds9_region(ds, reg, obj=None, field_parameters=None):
    r"""
    Create a data container from a ds9 region file. Requires the pyregion
    package (https://pyregion.readthedocs.io/en/latest/) to be installed.

    Parameters
    ----------
    ds : FITSDataset
        The Dataset to create the region from.
    reg : string
        The filename of the ds9 region, or a region string to be parsed.
    obj : data container, optional
        The data container that will be used to create the new region.
        Defaults to ds.all_data.
    field_parameters : dictionary, optional
        A set of field parameters to apply to the region.

    Examples
    --------

    >>> ds = yt.load("m33_hi.fits")
    >>> circle_region = ds9_region(ds, "circle.reg")
    >>> print(circle_region.quantities.extrema("flux"))
    """
    import pyregion

    from yt.frontends.fits.api import EventsFITSDataset

    if os.path.exists(reg):
        r = pyregion.open(reg)
    else:
        r = pyregion.parse(reg)
    reg_name = reg
    header = ds.wcs_2d.to_header()
    # The FITS header only contains WCS-related keywords
    header["NAXIS1"] = nx = ds.domain_dimensions[ds.lon_axis]
    header["NAXIS2"] = ny = ds.domain_dimensions[ds.lat_axis]
    filter = r.get_filter(header=header)
    mask = filter.mask((ny, nx)).transpose()
    if isinstance(ds, EventsFITSDataset):
        prefix = "event_"
    else:
        prefix = ""

    def _reg_field(field, data):
        i = data[prefix + "xyz"[ds.lon_axis]].d.astype("int") - 1
        j = data[prefix + "xyz"[ds.lat_axis]].d.astype("int") - 1
        new_mask = mask[i, j]
        ret = np.zeros(data[prefix + "x"].shape)
        ret[new_mask] = 1.0
        return ret

    ds.add_field(("gas", reg_name), sampling_type="cell", function=_reg_field)
    if obj is None:
        obj = ds.all_data()
    if field_parameters is not None:
        for k, v in field_parameters.items():
            obj.set_field_parameter(k, v)
    return obj.cut_region([f"obj['{reg_name}'] > 0"])


class PlotWindowWCS:
    r"""
    Use AstroPy's WCSAxes class to plot celestial coordinates on the axes of a
    on-axis PlotWindow plot. See
    http://docs.astropy.org/en/stable/visualization/wcsaxes/ for more details
    on how it works under the hood. This functionality requires a version of
    AstroPy >= 1.3.

    Parameters
    ----------
    pw : on-axis PlotWindow instance
        The PlotWindow instance to add celestial axes to.
    """

    def __init__(self, pw):
        WCSAxes = _astropy.wcsaxes.WCSAxes

        if pw.oblique:
            raise NotImplementedError("WCS axes are not implemented for oblique plots.")
        if not hasattr(pw.ds, "wcs_2d"):
            raise NotImplementedError("WCS axes are not implemented for this dataset.")
        if pw.data_source.axis != pw.ds.spec_axis:
            raise NotImplementedError("WCS axes are not implemented for this axis.")
        self.plots = {}
        self.pw = pw
        for f in pw.plots:
            rect = pw.plots[f]._get_best_layout()[1]
            fig = pw.plots[f].figure
            ax = fig.axes[0]
            wcs_ax = WCSAxes(fig, rect, wcs=pw.ds.wcs_2d, frameon=False)
            fig.add_axes(wcs_ax)
            wcs = pw.ds.wcs_2d.wcs
            xax = pw.ds.coordinates.x_axis[pw.data_source.axis]
            yax = pw.ds.coordinates.y_axis[pw.data_source.axis]
            xlabel = f"{wcs.ctype[xax].split('-')[0]} ({wcs.cunit[xax]})"
            ylabel = f"{wcs.ctype[yax].split('-')[0]} ({wcs.cunit[yax]})"
            fp = pw._font_properties
            wcs_ax.coords[0].set_axislabel(xlabel, fontproperties=fp, minpad=0.5)
            wcs_ax.coords[1].set_axislabel(ylabel, fontproperties=fp, minpad=0.4)
            wcs_ax.coords[0].ticklabels.set_fontproperties(fp)
            wcs_ax.coords[1].ticklabels.set_fontproperties(fp)
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            wcs_ax.set_xlim(pw.xlim[0].value, pw.xlim[1].value)
            wcs_ax.set_ylim(pw.ylim[0].value, pw.ylim[1].value)
            wcs_ax.coords.frame._update_cache = []
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            self.plots[f] = fig

    def keys(self):
        return self.plots.keys()

    def values(self):
        return self.plots.values()

    def items(self):
        return self.plots.items()

    def __getitem__(self, key):
        for k in self.keys():
            if k[1] == key:
                return self.plots[k]

    def show(self):
        return self

    def save(self, name=None, mpl_kwargs=None):
        if mpl_kwargs is None:
            mpl_kwargs = {}
        mpl_kwargs["bbox_inches"] = "tight"
        self.pw.save(name=name, mpl_kwargs=mpl_kwargs)

    def _repr_html_(self):
        from yt.visualization._mpl_imports import FigureCanvasAgg

        ret = ""
        for v in self.plots.values():
            canvas = FigureCanvasAgg(v)
            f = BytesIO()
            canvas.print_figure(f)
            f.seek(0)
            img = base64.b64encode(f.read()).decode()
            ret += (
                r'<img style="max-width:100%%;max-height:100%%;" '
                r'src="data:image/png;base64,%s"><br>' % img
            )
        return ret
