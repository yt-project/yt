import re
import sys
from functools import partial
from numbers import Number as numeric_type
from typing import Tuple

import numpy as np
from more_itertools import first, mark_ends

from yt._maintenance.deprecation import issue_deprecation_warning
from yt.data_objects.construction_data_containers import YTCoveringGrid
from yt.data_objects.image_array import ImageArray
from yt.fields.derived_field import DerivedField
from yt.funcs import fix_axis, is_sequence, iter_fields, mylog, validate_moment
from yt.units import dimensions
from yt.units.unit_object import Unit  # type: ignore
from yt.units.yt_array import YTArray, YTQuantity
from yt.utilities.math_utils import compute_stddev_image
from yt.utilities.on_demand_imports import _astropy
from yt.utilities.parallel_tools.parallel_analysis_interface import parallel_root_only
from yt.visualization.fixed_resolution import FixedResolutionBuffer, ParticleImageBuffer
from yt.visualization.particle_plots import ParticleAxisAlignedDummyDataSource
from yt.visualization.volume_rendering.off_axis_projection import off_axis_projection


class UnitfulHDU:
    def __init__(self, hdu):
        self.hdu = hdu
        self.header = self.hdu.header
        self.name = self.header["BTYPE"]
        self.units = self.header["BUNIT"]
        self.shape = self.hdu.shape

    @property
    def data(self):
        return YTArray(self.hdu.data, self.units)

    def __repr__(self):
        im_shape = " x ".join(str(s) for s in self.shape)
        return f"FITSImage: {self.name} ({im_shape}, {self.units})"


class FITSImageData:
    def __init__(
        self,
        data,
        fields=None,
        length_unit=None,
        width=None,
        img_ctr=None,
        wcs=None,
        current_time=None,
        time_unit=None,
        mass_unit=None,
        velocity_unit=None,
        magnetic_unit=None,
        ds=None,
        unit_header=None,
        **kwargs,
    ):
        r"""Initialize a FITSImageData object.

        FITSImageData contains a collection of FITS ImageHDU instances and
        WCS information, along with units for each of the images. FITSImageData
        instances can be constructed from ImageArrays, NumPy arrays, dicts
        of such arrays, FixedResolutionBuffers, and YTCoveringGrids. The latter
        two are the most powerful because WCS information can be constructed
        automatically from their coordinates.

        Parameters
        ----------
        data : FixedResolutionBuffer or a YTCoveringGrid. Or, an
            ImageArray, an numpy.ndarray, or dict of such arrays
            The data to be made into a FITS image or images.
        fields : single string or list of strings, optional
            The field names for the data. If *fields* is none and *data* has
            keys, it will use these for the fields. If *data* is just a
            single array one field name must be specified.
        length_unit : string
            The units of the WCS coordinates and the length unit of the file.
            Defaults to the length unit of the dataset, if there is one, or
            "cm" if there is not.
        width : float or YTQuantity
            The width of the image. Either a single value or iterable of values.
            If a float, assumed to be in *units*. Only used if this information
            is not already provided by *data*.
        img_ctr : array_like or YTArray
            The center coordinates of the image. If a list or NumPy array,
            it is assumed to be in *units*. Only used if this information
            is not already provided by *data*.
        wcs : `~astropy.wcs.WCS` instance, optional
            Supply an AstroPy WCS instance. Will override automatic WCS
            creation from FixedResolutionBuffers and YTCoveringGrids.
        current_time : float, tuple, or YTQuantity, optional
            The current time of the image(s). If not specified, one will
            be set from the dataset if there is one. If a float, it will
            be assumed to be in *time_unit* units.
        time_unit : string
            The default time units of the file. Defaults to "s".
        mass_unit : string
            The default time units of the file. Defaults to "g".
        velocity_unit : string
            The default velocity units of the file. Defaults to "cm/s".
        magnetic_unit : string
            The default magnetic units of the file. Defaults to "gauss".
        ds : `~yt.static_output.Dataset` instance, optional
            The dataset associated with the image(s), typically used
            to transfer metadata to the header(s). Does not need to be
            specified if *data* has a dataset as an attribute.

        Examples
        --------

        >>> # This example uses a FRB.
        >>> ds = load("sloshing_nomag2_hdf5_plt_cnt_0150")
        >>> prj = ds.proj(2, "kT", weight_field=("gas", "density"))
        >>> frb = prj.to_frb((0.5, "Mpc"), 800)
        >>> # This example just uses the FRB and puts the coords in kpc.
        >>> f_kpc = FITSImageData(
        ...     frb, fields="kT", length_unit="kpc", time_unit=(1.0, "Gyr")
        ... )
        >>> # This example specifies a specific WCS.
        >>> from astropy.wcs import WCS
        >>> w = WCS(naxis=self.dimensionality)
        >>> w.wcs.crval = [30.0, 45.0]  # RA, Dec in degrees
        >>> w.wcs.cunit = ["deg"] * 2
        >>> nx, ny = 800, 800
        >>> w.wcs.crpix = [0.5 * (nx + 1), 0.5 * (ny + 1)]
        >>> w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        >>> scale = 1.0 / 3600.0  # One arcsec per pixel
        >>> w.wcs.cdelt = [-scale, scale]
        >>> f_deg = FITSImageData(frb, fields="kT", wcs=w)
        >>> f_deg.writeto("temp.fits")
        """

        if fields is not None:
            fields = list(iter_fields(fields))

        if ds is None:
            ds = getattr(data, "ds", None)

        self.fields = []
        self.field_units = {}

        if unit_header is None:
            self._set_units(
                ds, [length_unit, mass_unit, time_unit, velocity_unit, magnetic_unit]
            )
        else:
            self._set_units_from_header(unit_header)

        wcs_unit = str(self.length_unit.units)

        self._fix_current_time(ds, current_time)

        if width is None:
            width = 1.0
        if isinstance(width, tuple):
            if ds is None:
                width = YTQuantity(width[0], width[1])
            else:
                width = ds.quan(width[0], width[1])
        if img_ctr is None:
            img_ctr = np.zeros(3)

        exclude_fields = [
            "x",
            "y",
            "z",
            "px",
            "py",
            "pz",
            "pdx",
            "pdy",
            "pdz",
            "weight_field",
        ]

        if isinstance(data, _astropy.pyfits.PrimaryHDU):
            data = _astropy.pyfits.HDUList([data])

        if isinstance(data, _astropy.pyfits.HDUList):
            self.hdulist = data
            for hdu in data:
                self.fields.append(hdu.header["btype"])
                self.field_units[hdu.header["btype"]] = hdu.header["bunit"]

            self.shape = self.hdulist[0].shape
            self.dimensionality = len(self.shape)
            wcs_names = [key for key in self.hdulist[0].header if "WCSNAME" in key]
            for name in wcs_names:
                if name == "WCSNAME":
                    key = " "
                else:
                    key = name[-1]
                w = _astropy.pywcs.WCS(
                    header=self.hdulist[0].header, key=key, naxis=self.dimensionality
                )
                setattr(self, "wcs" + key.strip().lower(), w)

            return

        self.hdulist = _astropy.pyfits.HDUList()

        stddev = False
        if hasattr(data, "keys"):
            img_data = data
            if fields is None:
                fields = list(img_data.keys())
            if hasattr(data, "data_source"):
                stddev = getattr(data.data_source, "moment", 1) == 2
        elif isinstance(data, np.ndarray):
            if fields is None:
                mylog.warning(
                    "No field name given for this array. Calling it 'image_data'."
                )
                fn = "image_data"
                fields = [fn]
            else:
                fn = fields[0]
            img_data = {fn: data}

        for fd in fields:
            if isinstance(fd, tuple):
                fname = fd[1]
            elif isinstance(fd, DerivedField):
                fname = fd.name[1]
            else:
                fname = fd
            if stddev:
                fname += "_stddev"
            self.fields.append(fname)

        # Sanity checking names
        s = set()
        duplicates = {f for f in self.fields if f in s or s.add(f)}
        if len(duplicates) > 0:
            for i, fd in enumerate(self.fields):
                if fd in duplicates:
                    if isinstance(fields[i], tuple):
                        ftype, fname = fields[i]
                    elif isinstance(fields[i], DerivedField):
                        ftype, fname = fields[i].name
                    else:
                        raise RuntimeError(
                            f"Cannot distinguish between fields with same name {fd}!"
                        )
                    self.fields[i] = f"{ftype}_{fname}"

        for is_first, _is_last, (i, (name, field)) in mark_ends(
            enumerate(zip(self.fields, fields))
        ):
            if name not in exclude_fields:
                this_img = img_data[field]
                if hasattr(img_data[field], "units"):
                    has_code_unit = False
                    for atom in this_img.units.expr.atoms():
                        if str(atom).startswith("code"):
                            has_code_unit = True
                    if has_code_unit:
                        mylog.warning(
                            "Cannot generate an image with code "
                            "units. Converting to units in CGS."
                        )
                        funits = this_img.units.get_base_equivalent("cgs")
                        this_img.convert_to_units(funits)
                    else:
                        funits = this_img.units
                    self.field_units[name] = str(funits)
                else:
                    self.field_units[name] = "dimensionless"
                mylog.info("Making a FITS image of field %s", name)
                if isinstance(this_img, ImageArray):
                    if i == 0:
                        self.shape = this_img.shape[::-1]
                    this_img = np.asarray(this_img)
                else:
                    if i == 0:
                        self.shape = this_img.shape
                    this_img = np.asarray(this_img.T)
                if is_first:
                    hdu = _astropy.pyfits.PrimaryHDU(this_img)
                else:
                    hdu = _astropy.pyfits.ImageHDU(this_img)
                hdu.name = name
                hdu.header["btype"] = name
                hdu.header["bunit"] = re.sub("()", "", self.field_units[name])
                for unit in ("length", "time", "mass", "velocity", "magnetic"):
                    if unit == "magnetic":
                        short_unit = "bf"
                    else:
                        short_unit = unit[0]
                    key = f"{short_unit}unit"
                    value = getattr(self, f"{unit}_unit")
                    if value is not None:
                        hdu.header[key] = float(value.value)
                        hdu.header.comments[key] = f"[{value.units}]"
                hdu.header["time"] = float(self.current_time.value)
                if hasattr(self, "current_redshift"):
                    hdu.header["HUBBLE"] = self.hubble_constant
                    hdu.header["REDSHIFT"] = self.current_redshift
                self.hdulist.append(hdu)

        self.dimensionality = len(self.shape)

        if wcs is None:
            w = _astropy.pywcs.WCS(
                header=self.hdulist[0].header, naxis=self.dimensionality
            )
            # FRBs and covering grids are special cases where
            # we have coordinate information, so we take advantage
            # of this and construct the WCS object
            if isinstance(img_data, FixedResolutionBuffer):
                dx = (img_data.bounds[1] - img_data.bounds[0]).to_value(wcs_unit)
                dy = (img_data.bounds[3] - img_data.bounds[2]).to_value(wcs_unit)
                dx /= self.shape[0]
                dy /= self.shape[1]
                xctr = 0.5 * (img_data.bounds[1] + img_data.bounds[0]).to_value(
                    wcs_unit
                )
                yctr = 0.5 * (img_data.bounds[3] + img_data.bounds[2]).to_value(
                    wcs_unit
                )
                center = [xctr, yctr]
                cdelt = [dx, dy]
            elif isinstance(img_data, YTCoveringGrid):
                cdelt = img_data.dds.to_value(wcs_unit)
                center = 0.5 * (img_data.left_edge + img_data.right_edge).to_value(
                    wcs_unit
                )
            else:
                # If img_data is just an array we use the width and img_ctr
                # parameters to determine the cell widths
                if not is_sequence(width):
                    width = [width] * self.dimensionality
                if isinstance(width[0], YTQuantity):
                    cdelt = [
                        wh.to_value(wcs_unit) / n for wh, n in zip(width, self.shape)
                    ]
                else:
                    cdelt = [float(wh) / n for wh, n in zip(width, self.shape)]
                center = img_ctr[: self.dimensionality]
            w.wcs.crpix = 0.5 * (np.array(self.shape) + 1)
            w.wcs.crval = center
            w.wcs.cdelt = cdelt
            w.wcs.ctype = ["linear"] * self.dimensionality
            w.wcs.cunit = [wcs_unit] * self.dimensionality
            self.set_wcs(w)
        else:
            self.set_wcs(wcs)

    def _fix_current_time(self, ds, current_time):
        if ds is None:
            registry = None
        else:
            registry = ds.unit_registry
        tunit = Unit(self.time_unit, registry=registry)
        if current_time is None:
            if ds is not None:
                current_time = ds.current_time
            else:
                self.current_time = YTQuantity(0.0, "s")
                return
        elif isinstance(current_time, numeric_type):
            current_time = YTQuantity(current_time, tunit)
        elif isinstance(current_time, tuple):
            current_time = YTQuantity(current_time[0], current_time[1])
        self.current_time = current_time.to(tunit)

    def _set_units(self, ds, base_units):
        if ds is not None:
            if getattr(ds, "cosmological_simulation", False):
                self.hubble_constant = ds.hubble_constant
                self.current_redshift = ds.current_redshift
        attrs = (
            "length_unit",
            "mass_unit",
            "time_unit",
            "velocity_unit",
            "magnetic_unit",
        )
        cgs_units = ("cm", "g", "s", "cm/s", "gauss")
        for unit, attr, cgs_unit in zip(base_units, attrs, cgs_units):
            if unit is None:
                if ds is not None:
                    u = getattr(ds, attr, None)
                elif attr == "velocity_unit":
                    u = self.length_unit / self.time_unit
                elif attr == "magnetic_unit":
                    u = np.sqrt(
                        4.0
                        * np.pi
                        * self.mass_unit
                        / (self.time_unit**2 * self.length_unit)
                    )
                else:
                    u = cgs_unit
            else:
                u = unit

            if isinstance(u, str):
                uq = YTQuantity(1.0, u)
            elif isinstance(u, numeric_type):
                uq = YTQuantity(u, cgs_unit)
            elif isinstance(u, YTQuantity):
                uq = u.copy()
            elif isinstance(u, tuple):
                uq = YTQuantity(u[0], u[1])
            else:
                uq = None

            if uq is not None:
                atoms = {str(a) for a in uq.units.expr.atoms()}
                if hasattr(self, "hubble_constant"):
                    # Don't store cosmology units
                    if str(uq.units).startswith("cm") or "h" in atoms or "a" in atoms:
                        uq.convert_to_cgs()
                if any(a.startswith("code") for a in atoms):
                    # Don't store code units
                    mylog.warning(
                        "Cannot use code units of '%s' "
                        "when creating a FITSImageData instance! "
                        "Converting to a cgs equivalent.",
                        uq.units,
                    )
                    uq.convert_to_cgs()

            if attr == "length_unit" and uq.value != 1.0:
                mylog.warning("Converting length units from %s to %s.", uq, uq.units)
                uq = YTQuantity(1.0, uq.units)

            setattr(self, attr, uq)

    def _set_units_from_header(self, header):
        if "hubble" in header:
            self.hubble_constant = header["HUBBLE"]
            self.current_redshift = header["REDSHIFT"]
        for unit in ["length", "time", "mass", "velocity", "magnetic"]:
            if unit == "magnetic":
                key = "BFUNIT"
            else:
                key = unit[0].upper() + "UNIT"
            if key not in header:
                continue
            u = header.comments[key].strip("[]")
            uq = YTQuantity(header[key], u)
            setattr(self, unit + "_unit", uq)

    def set_wcs(self, wcs, wcsname=None, suffix=None):
        """
        Set the WCS coordinate information for all images
        with a WCS object *wcs*.
        """
        if wcsname is None:
            wcs.wcs.name = "yt"
        else:
            wcs.wcs.name = wcsname
        h = wcs.to_header()
        if suffix is None:
            self.wcs = wcs
        else:
            setattr(self, "wcs" + suffix, wcs)
        for img in self.hdulist:
            for k, v in h.items():
                kk = k
                if suffix is not None:
                    kk += suffix
                img.header[kk] = v

    def change_image_name(self, old_name, new_name):
        """
        Change the name of a FITS image.

        Parameters
        ----------
        old_name : string
            The old name of the image.
        new_name : string
            The new name of the image.
        """
        idx = self.fields.index(old_name)
        self.hdulist[idx].name = new_name
        self.hdulist[idx].header["BTYPE"] = new_name
        self.field_units[new_name] = self.field_units.pop(old_name)
        self.fields[idx] = new_name

    def convolve(self, field, kernel, **kwargs):
        """
        Convolve an image with a kernel, either a simple
        Gaussian kernel or one provided by AstroPy. Currently,
        this only works for 2D images.

        All keyword arguments are passed to
        :meth:`~astropy.convolution.convolve`.

        Parameters
        ----------
        field : string
            The name of the field to convolve.
        kernel : float, YTQuantity, (value, unit) tuple, or AstroPy Kernel object
            The kernel to convolve the image with. If this is an AstroPy Kernel
            object, the image will be convolved with it. Otherwise, it is
            assumed that the kernel is a Gaussian and that this value is
            the standard deviation. If a float, it is assumed that the units
            are pixels, but a (value, unit) tuple or YTQuantity can be supplied
            to specify the standard deviation in physical units.

        Examples
        --------
        >>> fid = FITSSlice(ds, "z", ("gas", "density"))
        >>> fid.convolve("density", (3.0, "kpc"))
        """
        if self.dimensionality == 3:
            raise RuntimeError("Convolution currently only works for 2D FITSImageData!")
        conv = _astropy.conv
        if field not in self.keys():
            raise KeyError(f"{field} not an image!")
        idx = self.fields.index(field)
        if not isinstance(kernel, conv.Kernel):
            if not isinstance(kernel, numeric_type):
                unit = str(self.wcs.wcs.cunit[0])
                pix_scale = YTQuantity(self.wcs.wcs.cdelt[0], unit)
                if isinstance(kernel, tuple):
                    stddev = YTQuantity(kernel[0], kernel[1]).to(unit)
                else:
                    stddev = kernel.to(unit)
                kernel = stddev / pix_scale
            kernel = conv.Gaussian2DKernel(x_stddev=kernel)
        self.hdulist[idx].data = conv.convolve(self.hdulist[idx].data, kernel, **kwargs)

    def update_header(self, field, key, value):
        """
        Update the FITS header for *field* with a
        *key*, *value* pair. If *field* == "all", all
        headers will be updated.
        """
        if field == "all":
            for img in self.hdulist:
                img.header[key] = value
        else:
            if field not in self.keys():
                raise KeyError(f"{field} not an image!")
            idx = self.fields.index(field)
            self.hdulist[idx].header[key] = value

    def update_all_headers(self, key, value):
        issue_deprecation_warning(
            "update_all_headers is deprecated. "
            "Use update_header('all', key, value) instead.",
            since="3.3",
            removal="4.2",
        )
        self.update_header("all", key, value)

    def keys(self):
        return self.fields

    def has_key(self, key):
        return key in self.fields

    def values(self):
        return [self[k] for k in self.fields]

    def items(self):
        return [(k, self[k]) for k in self.fields]

    def __getitem__(self, item):
        return UnitfulHDU(self.hdulist[item])

    def __repr__(self):
        return str([self[k] for k in self.keys()])

    def info(self, output=None):
        """
        Summarize the info of the HDUs in this `FITSImageData`
        instance.

        Note that this function prints its results to the console---it
        does not return a value.

        Parameters
        ----------
        output : file, boolean, optional
            A file-like object to write the output to.  If `False`, does not
            output to a file and instead returns a list of tuples representing
            the FITSImageData info.  Writes to ``sys.stdout`` by default.
        """
        hinfo = self.hdulist.info(output=False)
        num_cols = len(hinfo[0])
        if output is None:
            output = sys.stdout
        if num_cols == 8:
            header = "No.    Name      Ver    Type      Cards   Dimensions   Format     Units"
            format = "{:3d}  {:10}  {:3} {:11}  {:5d}   {}   {}   {}"
        else:
            header = (
                "No.    Name         Type      Cards   Dimensions   Format     Units"
            )
            format = "{:3d}  {:10}  {:11}  {:5d}   {}   {}   {}"
        if self.hdulist._file is None:
            name = "(No file associated with this FITSImageData)"
        else:
            name = self.hdulist._file.name
        results = [f"Filename: {name}", header]
        for line in hinfo:
            units = self.field_units[self.hdulist[line[0]].header["btype"]]
            summary = tuple(list(line[:-1]) + [units])
            if output:
                results.append(format.format(*summary))
            else:
                results.append(summary)

        if output:
            output.write("\n".join(results))
            output.write("\n")
            output.flush()
        else:
            return results[2:]

    @parallel_root_only
    def writeto(self, fileobj, fields=None, overwrite=False, **kwargs):
        r"""
        Write all of the fields or a subset of them to a FITS file.

        Parameters
        ----------
        fileobj : string
            The name of the file to write to.
        fields : list of strings, optional
            The fields to write to the file. If not specified
            all of the fields in the buffer will be written.
        overwrite : boolean
            Whether or not to overwrite a previously existing file.
            Default: False
        **kwargs
            Additional keyword arguments are passed to
            :meth:`~astropy.io.fits.HDUList.writeto`.
        """
        if fields is None:
            hdus = self.hdulist
        else:
            hdus = _astropy.pyfits.HDUList()
            for field in fields:
                hdus.append(self.hdulist[field])
        hdus.writeto(fileobj, overwrite=overwrite, **kwargs)

    def to_glue(self, label="yt", data_collection=None):
        """
        Takes the data in the FITSImageData instance and exports it to
        Glue (http://glueviz.org) for interactive analysis. Optionally
        add a *label*. If you are already within the Glue environment, you
        can pass a *data_collection* object, otherwise Glue will be started.
        """
        from glue.core import Data, DataCollection
        from glue.core.coordinates import coordinates_from_header

        try:
            from glue.app.qt.application import GlueApplication
        except ImportError:
            from glue.qt.glue_application import GlueApplication

        image = Data(label=label)
        image.coords = coordinates_from_header(self.wcs.to_header())
        for k in self.fields:
            image.add_component(self[k].data, k)
        if data_collection is None:
            dc = DataCollection([image])
            app = GlueApplication(dc)
            app.start()
        else:
            data_collection.append(image)

    def to_aplpy(self, **kwargs):
        """
        Use APLpy (http://aplpy.github.io) for plotting. Returns an
        `aplpy.FITSFigure` instance. All keyword arguments are passed
        to the `aplpy.FITSFigure` constructor.
        """
        import aplpy

        return aplpy.FITSFigure(self.hdulist, **kwargs)

    def get_data(self, field):
        """
        Return the data array of the image corresponding to *field*
        with units attached. Deprecated.
        """
        return self[field].data

    def set_unit(self, field, units):
        """
        Set the units of *field* to *units*.
        """
        if field not in self.keys():
            raise KeyError(f"{field} not an image!")
        idx = self.fields.index(field)
        new_data = YTArray(self.hdulist[idx].data, self.field_units[field]).to(units)
        self.hdulist[idx].data = new_data.v
        self.hdulist[idx].header["bunit"] = units
        self.field_units[field] = units

    def pop(self, key):
        """
        Remove a field with name *key*
        and return it as a new FITSImageData
        instance.
        """
        if key not in self.keys():
            raise KeyError(f"{key} not an image!")
        idx = self.fields.index(key)
        im = self.hdulist.pop(idx)
        self.field_units.pop(key)
        self.fields.remove(key)
        f = _astropy.pyfits.PrimaryHDU(im.data, header=im.header)
        return FITSImageData(f, current_time=f.header["TIME"], unit_header=f.header)

    def close(self):
        self.hdulist.close()

    @classmethod
    def from_file(cls, filename):
        """
        Generate a FITSImageData instance from one previously written to
        disk.

        Parameters
        ----------
        filename : string
            The name of the file to open.
        """
        f = _astropy.pyfits.open(filename, lazy_load_hdus=False)
        return cls(f, current_time=f[0].header["TIME"], unit_header=f[0].header)

    @classmethod
    def from_images(cls, image_list):
        """
        Generate a new FITSImageData instance from a list of FITSImageData
        instances.

        Parameters
        ----------
        image_list : list of FITSImageData instances
            The images to be combined.
        """
        image_list = image_list if isinstance(image_list, list) else [image_list]
        first_image = first(image_list)

        w = first_image.wcs
        img_shape = first_image.shape
        data = []
        for fid in image_list:
            assert_same_wcs(w, fid.wcs)
            if img_shape != fid.shape:
                raise RuntimeError("Images do not have the same shape!")
            for hdu in fid.hdulist:
                if len(data) == 0:
                    data.append(_astropy.pyfits.PrimaryHDU(hdu.data, header=hdu.header))
                else:
                    data.append(_astropy.pyfits.ImageHDU(hdu.data, header=hdu.header))
        data = _astropy.pyfits.HDUList(data)
        return cls(
            data,
            current_time=first_image.current_time,
            unit_header=first_image[0].header,
        )

    def create_sky_wcs(
        self,
        sky_center,
        sky_scale,
        ctype=None,
        crota=None,
        cd=None,
        pc=None,
        wcsname="celestial",
        replace_old_wcs=True,
    ):
        """
        Takes a Cartesian WCS and converts it to one in a
        sky-based coordinate system.

        Parameters
        ----------
        sky_center : iterable of floats
            Reference coordinates of the WCS in degrees.
        sky_scale : tuple or YTQuantity
            Conversion between an angle unit and a length unit,
            e.g. (3.0, "arcsec/kpc")
        ctype : list of strings, optional
            The type of the coordinate system to create. Default:
            A "tangential" projection.
        crota : 2-element ndarray, optional
            Rotation angles between cartesian coordinates and
            the celestial coordinates.
        cd : 2x2-element ndarray, optional
            Dimensioned coordinate transformation matrix.
        pc : 2x2-element ndarray, optional
            Coordinate transformation matrix.
        wcsname : string, optional
            The name of the WCS to be stored in the FITS header.
        replace_old_wcs : boolean, optional
            Whether or not to overwrite the default WCS of the
            FITSImageData instance. If false, a second WCS will
            be added to the header. Default: True.
        """
        if ctype is None:
            ctype = ["RA---TAN", "DEC--TAN"]
        old_wcs = self.wcs
        naxis = old_wcs.naxis
        crval = [sky_center[0], sky_center[1]]
        if isinstance(sky_scale, YTQuantity):
            scaleq = sky_scale
        else:
            scaleq = YTQuantity(sky_scale[0], sky_scale[1])
        if scaleq.units.dimensions != dimensions.angle / dimensions.length:
            raise RuntimeError(
                f"sky_scale {sky_scale} not in correct dimensions of angle/length!"
            )
        deltas = old_wcs.wcs.cdelt
        units = [str(unit) for unit in old_wcs.wcs.cunit]
        new_dx = (YTQuantity(-deltas[0], units[0]) * scaleq).in_units("deg")
        new_dy = (YTQuantity(deltas[1], units[1]) * scaleq).in_units("deg")
        new_wcs = _astropy.pywcs.WCS(naxis=naxis)
        cdelt = [new_dx.v, new_dy.v]
        cunit = ["deg"] * 2
        if naxis == 3:
            crval.append(old_wcs.wcs.crval[2])
            cdelt.append(old_wcs.wcs.cdelt[2])
            ctype.append(old_wcs.wcs.ctype[2])
            cunit.append(old_wcs.wcs.cunit[2])
        new_wcs.wcs.crpix = old_wcs.wcs.crpix
        new_wcs.wcs.cdelt = cdelt
        new_wcs.wcs.crval = crval
        new_wcs.wcs.cunit = cunit
        new_wcs.wcs.ctype = ctype
        if crota is not None:
            new_wcs.wcs.crota = crota
        if cd is not None:
            new_wcs.wcs.cd = cd
        if pc is not None:
            new_wcs.wcs.cd = pc
        if replace_old_wcs:
            self.set_wcs(new_wcs, wcsname=wcsname)
        else:
            self.set_wcs(new_wcs, wcsname=wcsname, suffix="a")


class FITSImageBuffer(FITSImageData):
    pass


def sanitize_fits_unit(unit):
    if unit == "Mpc":
        mylog.info("Changing FITS file length unit to kpc.")
        unit = "kpc"
    elif unit == "au":
        unit = "AU"
    return unit


# This list allows one to determine which axes are the
# correct axes of the image in a right-handed coordinate
# system depending on which axis is sliced or projected
axis_wcs = [[1, 2], [2, 0], [0, 1]]


def construct_image(ds, axis, data_source, center, image_res, width, length_unit):
    if width is None:
        width = ds.domain_width[axis_wcs[axis]]
        unit = ds.get_smallest_appropriate_unit(width[0])
        mylog.info(
            "Making an image of the entire domain, "
            "so setting the center to the domain center."
        )
    else:
        width = ds.coordinates.sanitize_width(axis, width, None)
        unit = str(width[0].units)
    if is_sequence(image_res):
        nx, ny = image_res
    else:
        nx, ny = image_res, image_res
    dx = width[0] / nx
    dy = width[1] / ny
    crpix = [0.5 * (nx + 1), 0.5 * (ny + 1)]
    if unit == "unitary":
        unit = ds.get_smallest_appropriate_unit(ds.domain_width.max())
    elif unit == "code_length":
        unit = ds.get_smallest_appropriate_unit(ds.quan(1.0, "code_length"))
    unit = sanitize_fits_unit(unit)
    if length_unit is None:
        length_unit = unit
    if any(char.isdigit() for char in length_unit) and "*" in length_unit:
        length_unit = length_unit.split("*")[-1]
    cunit = [length_unit] * 2
    ctype = ["LINEAR"] * 2
    cdelt = [dx.in_units(length_unit), dy.in_units(length_unit)]
    if is_sequence(axis):
        crval = center.in_units(length_unit)
    else:
        crval = [center[idx].in_units(length_unit) for idx in axis_wcs[axis]]
    if hasattr(data_source, "to_frb"):
        if is_sequence(axis):
            frb = data_source.to_frb(width[0], (nx, ny), height=width[1])
        else:
            frb = data_source.to_frb(width[0], (nx, ny), center=center, height=width[1])
    elif isinstance(data_source, ParticleAxisAlignedDummyDataSource):
        axes = axis_wcs[axis]
        bounds = (
            center[axes[0]] - width[0] / 2,
            center[axes[0]] + width[0] / 2,
            center[axes[1]] - width[1] / 2,
            center[axes[1]] + width[1] / 2,
        )
        frb = ParticleImageBuffer(
            data_source, bounds, (nx, ny), periodic=all(ds.periodicity)
        )
    else:
        frb = None
    w = _astropy.pywcs.WCS(naxis=2)
    w.wcs.crpix = crpix
    w.wcs.cdelt = cdelt
    w.wcs.crval = crval
    w.wcs.cunit = cunit
    w.wcs.ctype = ctype
    return w, frb, length_unit


def assert_same_wcs(wcs1, wcs2):
    from numpy.testing import assert_allclose

    assert wcs1.naxis == wcs2.naxis
    for i in range(wcs1.naxis):
        assert wcs1.wcs.cunit[i] == wcs2.wcs.cunit[i]
        assert wcs1.wcs.ctype[i] == wcs2.wcs.ctype[i]
    assert_allclose(wcs1.wcs.crpix, wcs2.wcs.crpix)
    assert_allclose(wcs1.wcs.cdelt, wcs2.wcs.cdelt)
    assert_allclose(wcs1.wcs.crval, wcs2.wcs.crval)
    crota1 = getattr(wcs1.wcs, "crota", None)
    crota2 = getattr(wcs2.wcs, "crota", None)
    if crota1 is None or crota2 is None:
        assert crota1 == crota2
    else:
        assert_allclose(wcs1.wcs.crota, wcs2.wcs.crota)
    cd1 = getattr(wcs1.wcs, "cd", None)
    cd2 = getattr(wcs2.wcs, "cd", None)
    if cd1 is None or cd2 is None:
        assert cd1 == cd2
    else:
        assert_allclose(wcs1.wcs.cd, wcs2.wcs.cd)
    pc1 = getattr(wcs1.wcs, "pc", None)
    pc2 = getattr(wcs2.wcs, "pc", None)
    if pc1 is None or pc2 is None:
        assert pc1 == pc2
    else:
        assert_allclose(wcs1.wcs.pc, wcs2.wcs.pc)


class FITSSlice(FITSImageData):
    r"""
    Generate a FITSImageData of an on-axis slice.

    Parameters
    ----------
    ds : :class:`~yt.data_objects.static_output.Dataset`
        The dataset object.
    axis : character or integer
        The axis of the slice. One of "x","y","z", or 0,1,2.
    fields : string or list of strings
        The fields to slice
    image_res : an int or 2-tuple of ints
        Specify the resolution of the resulting image. A single value will be
        used for both axes, whereas a tuple of values will be used for the
        individual axes. Default: 512
    center : A sequence of floats, a string, or a tuple.
        The coordinate of the center of the image. If set to 'c', 'center' or
        left blank, the plot is centered on the middle of the domain. If set
        to 'max' or 'm', the center will be located at the maximum of the
        ('gas', 'density') field. Centering on the max or min of a specific
        field is supported by providing a tuple such as ("min","temperature")
        or ("max","dark_matter_density"). Units can be specified by passing in
        *center* as a tuple containing a coordinate and string unit name or by
        passing in a YTArray. If a list or unitless array is supplied, code
        units are assumed.
    width : tuple or a float.
        Width can have four different formats to support variable
        x and y widths.  They are:

        ==================================     =======================
        format                                 example
        ==================================     =======================
        (float, string)                        (10,'kpc')
        ((float, string), (float, string))     ((10,'kpc'),(15,'kpc'))
        float                                  0.2
        (float, float)                         (0.2, 0.3)
        ==================================     =======================

        For example, (10, 'kpc') specifies a width that is 10 kiloparsecs
        wide in the x and y directions, ((10,'kpc'),(15,'kpc')) specifies a
        width that is 10 kiloparsecs wide along the x axis and 15
        kiloparsecs wide along the y axis.  In the other two examples, code
        units are assumed, for example (0.2, 0.3) specifies a width that has an
        x width of 0.2 and a y width of 0.3 in code units.
    length_unit : string, optional
        the length units that the coordinates are written in. The default
        is to use the default length unit of the dataset.
    """

    def __init__(
        self,
        ds,
        axis,
        fields,
        image_res=512,
        center="c",
        width=None,
        length_unit=None,
        **kwargs,
    ):
        fields = list(iter_fields(fields))
        axis = fix_axis(axis, ds)
        center, dcenter = ds.coordinates.sanitize_center(center, axis)
        slc = ds.slice(axis, center[axis], **kwargs)
        w, frb, lunit = construct_image(
            ds, axis, slc, dcenter, image_res, width, length_unit
        )
        super().__init__(frb, fields=fields, length_unit=lunit, wcs=w)


class FITSProjection(FITSImageData):
    r"""
    Generate a FITSImageData of an on-axis projection.

    Parameters
    ----------
    ds : :class:`~yt.data_objects.static_output.Dataset`
        The dataset object.
    axis : character or integer
        The axis along which to project. One of "x","y","z", or 0,1,2.
    fields : string or list of strings
        The fields to project
    image_res : an int or 2-tuple of ints
        Specify the resolution of the resulting image. A single value will be
        used for both axes, whereas a tuple of values will be used for the
        individual axes. Default: 512
    center : A sequence of floats, a string, or a tuple.
        The coordinate of the center of the image. If set to 'c', 'center' or
        left blank, the plot is centered on the middle of the domain. If set
        to 'max' or 'm', the center will be located at the maximum of the
        ('gas', 'density') field. Centering on the max or min of a specific
        field is supported by providing a tuple such as ("min","temperature")
        or ("max","dark_matter_density"). Units can be specified by passing in
        *center* as a tuple containing a coordinate and string unit name or by
        passing in a YTArray. If a list or unitless array is supplied, code
        units are assumed.
    width : tuple or a float.
        Width can have four different formats to support variable
        x and y widths.  They are:

        ==================================     =======================
        format                                 example
        ==================================     =======================
        (float, string)                        (10,'kpc')
        ((float, string), (float, string))     ((10,'kpc'),(15,'kpc'))
        float                                  0.2
        (float, float)                         (0.2, 0.3)
        ==================================     =======================

        For example, (10, 'kpc') specifies a width that is 10 kiloparsecs
        wide in the x and y directions, ((10,'kpc'),(15,'kpc')) specifies a
        width that is 10 kiloparsecs wide along the x axis and 15
        kiloparsecs wide along the y axis.  In the other two examples, code
        units are assumed, for example (0.2, 0.3) specifies a width that has an
        x width of 0.2 and a y width of 0.3 in code units.
    weight_field : string
        The field used to weight the projection.
    length_unit : string, optional
        the length units that the coordinates are written in. The default
        is to use the default length unit of the dataset.
    moment : integer, optional
        for a weighted projection, moment = 1 (the default) corresponds to a
        weighted average. moment = 2 corresponds to a weighted standard
        deviation.
    """

    def __init__(
        self,
        ds,
        axis,
        fields,
        image_res=512,
        center="c",
        width=None,
        weight_field=None,
        length_unit=None,
        *,
        moment=1,
        **kwargs,
    ):
        fields = list(iter_fields(fields))
        axis = fix_axis(axis, ds)
        center, dcenter = ds.coordinates.sanitize_center(center, axis)
        prj = ds.proj(
            fields[0], axis, weight_field=weight_field, moment=moment, **kwargs
        )
        w, frb, lunit = construct_image(
            ds, axis, prj, dcenter, image_res, width, length_unit
        )
        super().__init__(frb, fields=fields, length_unit=lunit, wcs=w)


class FITSParticleProjection(FITSImageData):
    r"""
    Generate a FITSImageData of an on-axis projection of a
    particle field.

    Parameters
    ----------
    ds : :class:`~yt.data_objects.static_output.Dataset`
        The dataset object.
    axis : character or integer
        The axis along which to project. One of "x","y","z", or 0,1,2.
    fields : string or list of strings
        The fields to project
    image_res : an int or 2-tuple of ints
        Specify the resolution of the resulting image. A single value will be
        used for both axes, whereas a tuple of values will be used for the
        individual axes. Default: 512
    center : A sequence of floats, a string, or a tuple.
        The coordinate of the center of the image. If set to 'c', 'center' or
        left blank, the plot is centered on the middle of the domain. If set
        to 'max' or 'm', the center will be located at the maximum of the
        ('gas', 'density') field. Centering on the max or min of a specific
        field is supported by providing a tuple such as ("min","temperature")
        or ("max","dark_matter_density"). Units can be specified by passing in
        *center* as a tuple containing a coordinate and string unit name or by
        passing in a YTArray. If a list or unitless array is supplied, code
        units are assumed.
    width : tuple or a float.
        Width can have four different formats to support variable
        x and y widths.  They are:

        ==================================     =======================
        format                                 example
        ==================================     =======================
        (float, string)                        (10,'kpc')
        ((float, string), (float, string))     ((10,'kpc'),(15,'kpc'))
        float                                  0.2
        (float, float)                         (0.2, 0.3)
        ==================================     =======================

        For example, (10, 'kpc') specifies a width that is 10 kiloparsecs
        wide in the x and y directions, ((10,'kpc'),(15,'kpc')) specifies a
        width that is 10 kiloparsecs wide along the x axis and 15
        kiloparsecs wide along the y axis.  In the other two examples, code
        units are assumed, for example (0.2, 0.3) specifies a width that has an
        x width of 0.2 and a y width of 0.3 in code units.
    depth : A tuple or a float
         A tuple containing the depth to project through and the string
         key of the unit: (width, 'unit').  If set to a float, code units
         are assumed. Defaults to the entire domain.
    weight_field : string
        The field used to weight the projection.
    length_unit : string, optional
        the length units that the coordinates are written in. The default
        is to use the default length unit of the dataset.
    deposition : string, optional
        Controls the order of the interpolation of the particles onto the
        mesh. "ngp" is 0th-order "nearest-grid-point" method (the default),
        "cic" is 1st-order "cloud-in-cell".
    density : boolean, optional
        If True, the quantity to be projected will be divided by the area of
        the cells, to make a projected density of the quantity. Default: False
    field_parameters : dictionary
         A dictionary of field parameters than can be accessed by derived
         fields.
    data_source : yt.data_objects.data_containers.YTSelectionContainer, optional
        If specified, this will be the data source used for selecting regions
        to project.
    """

    def __init__(
        self,
        ds,
        axis,
        fields,
        image_res=512,
        center="c",
        width=None,
        depth=(1, "1"),
        weight_field=None,
        length_unit=None,
        deposition="ngp",
        density=False,
        field_parameters=None,
        data_source=None,
    ):
        fields = list(iter_fields(fields))
        axis = fix_axis(axis, ds)
        center, dcenter = ds.coordinates.sanitize_center(center, axis)
        width = ds.coordinates.sanitize_width(axis, width, depth)
        width[-1].convert_to_units(width[0].units)

        if field_parameters is None:
            field_parameters = {}

        ps = ParticleAxisAlignedDummyDataSource(
            center,
            ds,
            axis,
            width,
            fields,
            weight_field,
            field_parameters=field_parameters,
            data_source=data_source,
            deposition=deposition,
            density=density,
        )
        w, frb, lunit = construct_image(
            ds, axis, ps, dcenter, image_res, width, length_unit
        )
        super().__init__(frb, fields=fields, length_unit=lunit, wcs=w)


class FITSOffAxisSlice(FITSImageData):
    r"""
    Generate a FITSImageData of an off-axis slice.

    Parameters
    ----------
    ds : :class:`~yt.data_objects.static_output.Dataset`
        The dataset object.
    normal : a sequence of floats
        The vector normal to the projection plane.
    fields : string or list of strings
        The fields to slice
    image_res : an int or 2-tuple of ints
        Specify the resolution of the resulting image. A single value will be
        used for both axes, whereas a tuple of values will be used for the
        individual axes. Default: 512
    center : A sequence of floats, a string, or a tuple.
        The coordinate of the center of the image. If set to 'c', 'center' or
        left blank, the plot is centered on the middle of the domain. If set
        to 'max' or 'm', the center will be located at the maximum of the
        ('gas', 'density') field. Centering on the max or min of a specific
        field is supported by providing a tuple such as ("min","temperature")
        or ("max","dark_matter_density"). Units can be specified by passing in
        *center* as a tuple containing a coordinate and string unit name or by
        passing in a YTArray. If a list or unitless array is supplied, code
        units are assumed.
    width : tuple or a float.
        Width can have four different formats to support variable
        x and y widths.  They are:

        ==================================     =======================
        format                                 example
        ==================================     =======================
        (float, string)                        (10,'kpc')
        ((float, string), (float, string))     ((10,'kpc'),(15,'kpc'))
        float                                  0.2
        (float, float)                         (0.2, 0.3)
        ==================================     =======================

        For example, (10, 'kpc') specifies a width that is 10 kiloparsecs
        wide in the x and y directions, ((10,'kpc'),(15,'kpc')) specifies a
        width that is 10 kiloparsecs wide along the x axis and 15
        kiloparsecs wide along the y axis.  In the other two examples, code
        units are assumed, for example (0.2, 0.3) specifies a width that has an
        x width of 0.2 and a y width of 0.3 in code units.
    north_vector : a sequence of floats
        A vector defining the 'up' direction in the plot.  This
        option sets the orientation of the slicing plane.  If not
        set, an arbitrary grid-aligned north-vector is chosen.
    length_unit : string, optional
        the length units that the coordinates are written in. The default
        is to use the default length unit of the dataset.
    """

    def __init__(
        self,
        ds,
        normal,
        fields,
        image_res=512,
        center="c",
        width=None,
        north_vector=None,
        length_unit=None,
    ):
        fields = list(iter_fields(fields))
        center, dcenter = ds.coordinates.sanitize_center(center, 4)
        cut = ds.cutting(normal, center, north_vector=north_vector)
        center = ds.arr([0.0] * 2, "code_length")
        w, frb, lunit = construct_image(
            ds, normal, cut, center, image_res, width, length_unit
        )
        super().__init__(frb, fields=fields, length_unit=lunit, wcs=w)


class FITSOffAxisProjection(FITSImageData):
    r"""
    Generate a FITSImageData of an off-axis projection.

    Parameters
    ----------
    ds : :class:`~yt.data_objects.static_output.Dataset`
        This is the dataset object corresponding to the
        simulation output to be plotted.
    normal : a sequence of floats
        The vector normal to the projection plane.
    fields : string, list of strings
        The name of the field(s) to be plotted.
    image_res : an int or 2-tuple of ints
        Specify the resolution of the resulting image. A single value will be
        used for both axes, whereas a tuple of values will be used for the
        individual axes. Default: 512
    center : A sequence of floats, a string, or a tuple.
        The coordinate of the center of the image. If set to 'c', 'center' or
        left blank, the plot is centered on the middle of the domain. If set
        to 'max' or 'm', the center will be located at the maximum of the
        ('gas', 'density') field. Centering on the max or min of a specific
        field is supported by providing a tuple such as ("min","temperature")
        or ("max","dark_matter_density"). Units can be specified by passing in
        *center* as a tuple containing a coordinate and string unit name or by
        passing in a YTArray. If a list or unitless array is supplied, code
        units are assumed.
    width : tuple or a float.
        Width can have four different formats to support variable
        x and y widths.  They are:

        ==================================     =======================
        format                                 example
        ==================================     =======================
        (float, string)                        (10,'kpc')
        ((float, string), (float, string))     ((10,'kpc'),(15,'kpc'))
        float                                  0.2
        (float, float)                         (0.2, 0.3)
        ==================================     =======================

        For example, (10, 'kpc') specifies a width that is 10 kiloparsecs
        wide in the x and y directions, ((10,'kpc'),(15,'kpc')) specifies a
        width that is 10 kiloparsecs wide along the x axis and 15
        kiloparsecs wide along the y axis.  In the other two examples, code
        units are assumed, for example (0.2, 0.3) specifies a width that has an
        x width of 0.2 and a y width of 0.3 in code units.
    depth : A tuple or a float
        A tuple containing the depth to project through and the string
        key of the unit: (width, 'unit').  If set to a float, code units
        are assumed
    weight_field : string
         The name of the weighting field.  Set to None for no weight.
    north_vector : a sequence of floats
        A vector defining the 'up' direction in the plot.  This
        option sets the orientation of the slicing plane.  If not
        set, an arbitrary grid-aligned north-vector is chosen.
    method : string
        The method of projection.  Valid methods are:

        "integrate" with no weight_field specified : integrate the requested
        field along the line of sight.

        "integrate" with a weight_field specified : weight the requested
        field by the weighting field and integrate along the line of sight.

        "sum" : This method is the same as integrate, except that it does not
        multiply by a path length when performing the integration, and is
        just a straight summation of the field along the given axis. WARNING:
        This should only be used for uniform resolution grid datasets, as other
        datasets may result in unphysical images.
    data_source : yt.data_objects.data_containers.YTSelectionContainer, optional
        If specified, this will be the data source used for selecting regions
        to project.
    length_unit : string, optional
        the length units that the coordinates are written in. The default
        is to use the default length unit of the dataset.
    moment : integer, optional
        for a weighted projection, moment = 1 (the default) corresponds to a
        weighted average. moment = 2 corresponds to a weighted standard
        deviation.
    """

    def __init__(
        self,
        ds,
        normal,
        fields,
        center="c",
        width=(1.0, "unitary"),
        weight_field=None,
        image_res=512,
        data_source=None,
        north_vector=None,
        depth=(1.0, "unitary"),
        method="integrate",
        length_unit=None,
        *,
        moment=1,
    ):
        validate_moment(moment, weight_field)
        center, dcenter = ds.coordinates.sanitize_center(center, 4)
        buf = {}
        width = ds.coordinates.sanitize_width(normal, width, depth)
        wd = tuple(el.in_units("code_length").v for el in width)
        if not is_sequence(image_res):
            image_res = (image_res, image_res)
        res = (image_res[0], image_res[1])
        if data_source is None:
            source = ds.all_data()
        else:
            source = data_source
        fields = source._determine_fields(list(iter_fields(fields)))
        stddev_str = "_stddev" if moment == 2 else ""
        for item in fields:

            key = (item[0], item[1] + stddev_str)

            buf[key] = off_axis_projection(
                source,
                center,
                normal,
                wd,
                res,
                item,
                north_vector=north_vector,
                method=method,
                weight=weight_field,
            ).swapaxes(0, 1)

            if moment == 2:

                def _sq_field(field, data, item: Tuple[str, str]):
                    return data[item] ** 2

                fd = ds._get_field_info(*item)

                field_sq = (item[0], f"tmp_{item[1]}_squared")

                ds.add_field(
                    field_sq,
                    partial(_sq_field, item=item),
                    sampling_type=fd.sampling_type,
                    units=f"({fd.units})*({fd.units})",
                )

                buff2 = off_axis_projection(
                    source,
                    center,
                    normal,
                    wd,
                    res,
                    field_sq,
                    north_vector=north_vector,
                    method=method,
                    weight=weight_field,
                ).swapaxes(0, 1)

                buf[key] = compute_stddev_image(buff2, buf[key])

                ds.field_info.pop(field_sq)

        center = ds.arr([0.0] * 2, "code_length")
        w, not_an_frb, lunit = construct_image(
            ds, normal, buf, center, image_res, width, length_unit
        )
        super().__init__(buf, fields=list(buf.keys()), wcs=w, length_unit=lunit, ds=ds)
