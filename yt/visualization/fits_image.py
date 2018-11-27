"""
FITSImageData Class
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
from yt.extern.six import string_types
import numpy as np
from yt.fields.derived_field import DerivedField
from yt.funcs import mylog, iterable, fix_axis, ensure_list, \
    issue_deprecation_warning
from yt.visualization.fixed_resolution import FixedResolutionBuffer
from yt.data_objects.construction_data_containers import YTCoveringGrid
from yt.utilities.on_demand_imports import _astropy
from yt.units.yt_array import YTQuantity, YTArray
from yt.units import dimensions
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_root_only
from yt.visualization.volume_rendering.off_axis_projection import off_axis_projection
import re
import sys

class UnitfulHDU(object):
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
        im_shape = " x ".join([str(s) for s in self.shape])
        return "FITSImage: %s (%s, %s)" % (self.name, im_shape, self.units)

class FITSImageData(object):

    def __init__(self, data, fields=None, units=None, width=None, wcs=None):
        r""" Initialize a FITSImageData object.

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
        units : string
            The units of the WCS coordinates. Defaults to "cm".
        width : float or YTQuantity
            The width of the image. Either a single value or iterable of values.
            If a float, assumed to be in *units*. Only used if this information 
            is not already provided by *data*.
        wcs : `astropy.wcs.WCS` instance, optional
            Supply an AstroPy WCS instance. Will override automatic WCS
            creation from FixedResolutionBuffers and YTCoveringGrids.

        Examples
        --------

        >>> # This example uses a FRB.
        >>> ds = load("sloshing_nomag2_hdf5_plt_cnt_0150")
        >>> prj = ds.proj(2, "kT", weight_field="density")
        >>> frb = prj.to_frb((0.5, "Mpc"), 800)
        >>> # This example just uses the FRB and puts the coords in kpc.
        >>> f_kpc = FITSImageData(frb, fields="kT", units="kpc")
        >>> # This example specifies a specific WCS.
        >>> from astropy.wcs import WCS
        >>> w = WCS(naxis=self.dimensionality)
        >>> w.wcs.crval = [30., 45.] # RA, Dec in degrees
        >>> w.wcs.cunit = ["deg"]*2
        >>> nx, ny = 800, 800
        >>> w.wcs.crpix = [0.5*(nx+1), 0.5*(ny+1)]
        >>> w.wcs.ctype = ["RA---TAN","DEC--TAN"]
        >>> scale = 1./3600. # One arcsec per pixel
        >>> w.wcs.cdelt = [-scale, scale]
        >>> f_deg = FITSImageData(frb, fields="kT", wcs=w)
        >>> f_deg.writeto("temp.fits")
        """

        self.fields = []
        self.field_units = {}

        if units is None:
            units = "cm"
        if width is None:
            width = 1.0

        exclude_fields = ['x','y','z','px','py','pz',
                          'pdx','pdy','pdz','weight_field']

        if isinstance(data, _astropy.pyfits.HDUList):
            self.hdulist = data
            for hdu in data:
                self.fields.append(hdu.header["btype"])
                self.field_units[hdu.header["btype"]] = hdu.header['bunit']

            self.shape = self.hdulist[0].shape
            self.dimensionality = len(self.shape)
            wcs_names = [key for key in self.hdulist[0].header 
                         if "WCSNAME" in key]
            for name in wcs_names:
                if name == "WCSNAME":
                    key = ' '
                else:
                    key = name[-1]
                w = _astropy.pywcs.WCS(header=self.hdulist[0].header,
                                       key=key, naxis=self.dimensionality)
                setattr(self, "wcs"+key.strip().lower(), w)

            return

        self.hdulist = _astropy.pyfits.HDUList()

        if isinstance(fields, string_types):
            fields = [fields]

        if hasattr(data, 'keys'):
            img_data = data
            if fields is None:
                fields = list(img_data.keys())
        elif isinstance(data, np.ndarray):
            if fields is None:
                mylog.warning("No field name given for this array. "
                              "Calling it 'image_data'.")
                fn = 'image_data'
                fields = [fn]
            else:
                fn = fields[0]
            img_data = {fn: data}

        for fd in fields:
            if isinstance(fd, tuple):
                self.fields.append(fd[1])
            elif isinstance(fd, DerivedField):
                self.fields.append(fd.name[1])
            else:
                self.fields.append(fd)

        first = True
        for name, field in zip(self.fields, fields):
            if name not in exclude_fields:
                if hasattr(img_data[field], "units"):
                    self.field_units[name] = str(img_data[field].units)
                else:
                    self.field_units[name] = "dimensionless"
                mylog.info("Making a FITS image of field %s" % name)
                if first:
                    hdu = _astropy.pyfits.PrimaryHDU(np.array(img_data[field]))
                    first = False
                else:
                    hdu = _astropy.pyfits.ImageHDU(np.array(img_data[field]))
                hdu.name = name
                hdu.header["btype"] = name
                hdu.header["bunit"] = re.sub('()', '', self.field_units[name])
                self.hdulist.append(hdu)

        self.shape = self.hdulist[0].shape
        self.dimensionality = len(self.shape)

        if wcs is None:
            w = _astropy.pywcs.WCS(header=self.hdulist[0].header,
                                   naxis=self.dimensionality)
            if isinstance(img_data, FixedResolutionBuffer):
                # FRBs are a special case where we have coordinate
                # information, so we take advantage of this and
                # construct the WCS object
                dx = (img_data.bounds[1]-img_data.bounds[0]).to(units).v
                dy = (img_data.bounds[3]-img_data.bounds[2]).to(units).v
                dx /= self.shape[0]
                dy /= self.shape[1]
                xctr = 0.5*(img_data.bounds[1]+img_data.bounds[0]).to(units).v
                yctr = 0.5*(img_data.bounds[3]+img_data.bounds[2]).to(units).v
                center = [xctr, yctr]
                cdelt = [dx, dy]
            elif isinstance(img_data, YTCoveringGrid):
                cdelt = img_data.dds.to(units).v
                center = 0.5*(img_data.left_edge+img_data.right_edge).to(units).v
            else:
                # If img_data is just an array, we assume the center is the 
                # origin and use the image width to determine the cell widths
                if not iterable(width):
                    width = [width]*self.dimensionality
                if isinstance(width[0], YTQuantity):
                    cdelt = [wh.to(units).v/n for wh, n in zip(width, self.shape)]
                else:
                    cdelt = [float(wh)/n for wh, n in zip(width, self.shape)]
                center = [0.0]*self.dimensionality
            w.wcs.crpix = 0.5*(np.array(self.shape)+1)
            w.wcs.crval = center
            w.wcs.cdelt = cdelt
            w.wcs.ctype = ["linear"]*self.dimensionality
            w.wcs.cunit = [units]*self.dimensionality
            self.set_wcs(w)
        else:
            self.set_wcs(wcs)

    def set_wcs(self, wcs, wcsname=None, suffix=None):
        """
        Set the WCS coordinate information for all images
        with a WCS object *wcs*.
        """
        if wcsname is None:
            wcs.wcs.name = 'yt'
        else:
            wcs.wcs.name = wcsname
        h = wcs.to_header()
        if suffix is None:
            self.wcs = wcs
        else:
            setattr(self, "wcs"+suffix, wcs)
        for img in self.hdulist:
            for k, v in h.items():
                kk = k
                if suffix is not None:
                    kk += suffix
                img.header[kk] = v

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
                raise KeyError("%s not an image!" % field)
            idx = self.fields.index(field)
            self.hdulist[idx].header[key] = value

    def update_all_headers(self, key, value):
        mylog.warning("update_all_headers is deprecated. "
                      "Use update_header('all', key, value) instead.")
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
            header = 'No.    Name      Ver    Type      Cards   Dimensions   Format     Units'
            format = '{:3d}  {:10}  {:3} {:11}  {:5d}   {}   {}   {}'
        else:
            header = 'No.    Name         Type      Cards   Dimensions   Format     Units'
            format = '{:3d}  {:10}  {:11}  {:5d}   {}   {}   {}'
        if self.hdulist._file is None:
            name = '(No file associated with this FITSImageData)'
        else:
            name = self.hdulist._file.name
        results = ['Filename: {}'.format(name), header]
        for line in hinfo:
            units = self.field_units[self.hdulist[line[0]].header['btype']]
            summary = tuple(list(line[:-1]) + [units])
            if output:
                results.append(format.format(*summary))
            else:
                results.append(summary)

        if output:
            output.write('\n'.join(results))
            output.write('\n')
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
        clobber : overwrite, optional
            Whether or not to overwrite a previously existing file.
            Default: False

        Additional keyword arguments are passed to
        :meth:`~astropy.io.fits.HDUList.writeto`.
        """
        if "clobber" in kwargs:
            issue_deprecation_warning("The \"clobber\" keyword argument "
                                      "is deprecated. Use the \"overwrite\" "
                                      "argument, which has the same effect, "
                                      "instead.")
            overwrite = kwargs.pop("clobber")
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
        Glue (http://www.glueviz.org) for interactive analysis. Optionally 
        add a *label*. If you are already within the Glue environment, you 
        can pass a *data_collection* object, otherwise Glue will be started.
        """
        from glue.core import DataCollection, Data
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
            raise KeyError("%s not an image!" % field)
        idx = self.fields.index(field)
        new_data = YTArray(self.hdulist[idx].data, 
                           self.field_units[field]).to(units)
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
            raise KeyError("%s not an image!" % key)
        idx = self.fields.index(key)
        im = self.hdulist.pop(idx)
        self.field_units.pop(key)
        self.fields.remove(key)
        data = _astropy.pyfits.HDUList([im])
        return FITSImageData(data)

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
        return cls(f)

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
        w = image_list[0].wcs
        img_shape = image_list[0].shape
        data = []
        first = True
        for fid in image_list:
            assert_same_wcs(w, fid.wcs)
            if img_shape != fid.shape:
                raise RuntimeError("Images do not have the same shape!")
            for hdu in fid.hdulist:
                if first:
                    data.append(hdu)
                    first = False
                else:
                    data.append(_astropy.pyfits.ImageHDU(hdu.data, header=hdu.header))
        data = _astropy.pyfits.HDUList(data)
        return cls(data)

    def create_sky_wcs(self, sky_center, sky_scale,
                       ctype=["RA---TAN","DEC--TAN"],
                       crota=None, cd=None, pc=None,
                       wcsname="celestial",
                       replace_old_wcs=True):
        """
        Takes a Cartesian WCS and converts it to one in a
        celestial coordinate system.

        Parameters
        ----------
        sky_center : iterable of floats
            Reference coordinates of the WCS in degrees.
        sky_scale : tuple or YTQuantity
            Conversion between an angle unit and a length unit,
            e.g. (3.0, "arcsec/kpc")
        ctype : list of strings, optional
            The type of the coordinate system to create.
        crota : 2-element ndarray, optional
            Rotation angles between cartesian coordinates and
            the celestial coordinates.
        cd : 2x2-element ndarray, optional
            Dimensioned coordinate transformation matrix.
        pc : 2x2-element ndarray, optional
            Coordinate transformation matrix.
        replace_old_wcs : boolean, optional
            Whether or not to overwrite the default WCS of the 
            FITSImageData instance. If false, a second WCS will 
            be added to the header. Default: True.
        """
        old_wcs = self.wcs
        naxis = old_wcs.naxis
        crval = [sky_center[0], sky_center[1]]
        if isinstance(sky_scale, YTQuantity):
            scaleq = sky_scale
        else:
            scaleq = YTQuantity(sky_scale[0],sky_scale[1])
        if scaleq.units.dimensions != dimensions.angle/dimensions.length:
            raise RuntimeError("sky_scale %s not in correct " % sky_scale +
                               "dimensions of angle/length!")
        deltas = old_wcs.wcs.cdelt
        units = [str(unit) for unit in old_wcs.wcs.cunit]
        new_dx = (YTQuantity(-deltas[0], units[0])*scaleq).in_units("deg")
        new_dy = (YTQuantity(deltas[1], units[1])*scaleq).in_units("deg")
        new_wcs = _astropy.pywcs.WCS(naxis=naxis)
        cdelt = [new_dx.v, new_dy.v]
        cunit = ["deg"]*2
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
        mylog.info("Changing FITS file unit to kpc.")
        unit = "kpc"
    elif unit == "au":
        unit = "AU"
    return unit

axis_wcs = [[1,2],[0,2],[0,1]]

def construct_image(ds, axis, data_source, center, width=None, image_res=None):
    if width is None:
        width = ds.domain_width[axis_wcs[axis]]
        unit = ds.get_smallest_appropriate_unit(width[0])
        mylog.info("Making an image of the entire domain, "+
                   "so setting the center to the domain center.")
    else:
        width = ds.coordinates.sanitize_width(axis, width, None)
        unit = str(width[0].units)
    if image_res is None:
        ddims = ds.domain_dimensions*ds.refine_by**ds.index.max_level
        if iterable(axis):
            nx = ddims.max()
            ny = ddims.max()
        else:
            nx, ny = [ddims[idx] for idx in axis_wcs[axis]]
    else:
        if iterable(image_res):
            nx, ny = image_res
        else:
            nx, ny = image_res, image_res
    dx = width[0]/nx
    dy = width[1]/ny
    crpix = [0.5*(nx+1), 0.5*(ny+1)]
    if unit == "unitary":
        unit = ds.get_smallest_appropriate_unit(ds.domain_width.max())
    elif unit == "code_length":
        unit = ds.get_smallest_appropriate_unit(ds.quan(1.0, "code_length"))
    unit = sanitize_fits_unit(unit)
    cunit = [unit]*2
    ctype = ["LINEAR"]*2
    cdelt = [dx.in_units(unit), dy.in_units(unit)]
    if iterable(axis):
        crval = center.in_units(unit)
    else:
        crval = [center[idx].in_units(unit) for idx in axis_wcs[axis]]
    if hasattr(data_source, 'to_frb'):
        if iterable(axis):
            frb = data_source.to_frb(width[0], (nx, ny), height=width[1])
        else:
            frb = data_source.to_frb(width[0], (nx, ny), center=center,
                                     height=width[1])
    else:
        frb = None
    w = _astropy.pywcs.WCS(naxis=2)
    w.wcs.crpix = crpix
    w.wcs.cdelt = cdelt
    w.wcs.crval = crval
    w.wcs.cunit = cunit
    w.wcs.ctype = ctype
    return w, frb

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
         Width can have four different formats to support windows with variable
         x and y widths.  They are:

         ==================================     =======================
         format                                 example
         ==================================     =======================
         (float, string)                        (10,'kpc')
         ((float, string), (float, string))     ((10,'kpc'),(15,'kpc'))
         float                                  0.2
         (float, float)                         (0.2, 0.3)
         ==================================     =======================

         For example, (10, 'kpc') requests a plot window that is 10 kiloparsecs
         wide in the x and y directions, ((10,'kpc'),(15,'kpc')) requests a
         window that is 10 kiloparsecs wide along the x axis and 15
         kiloparsecs wide along the y axis.  In the other two examples, code
         units are assumed, for example (0.2, 0.3) requests a plot that has an
         x width of 0.2 and a y width of 0.3 in code units.  If units are
         provided the resulting plot axis labels will use the supplied units.
    image_res : an int or 2-tuple of ints
        Specify the resolution of the resulting image. If not provided, it will
        be determined based on the minimum cell size of the dataset.
    """
    def __init__(self, ds, axis, fields, center="c", width=None, 
                 image_res=None, **kwargs):
        fields = ensure_list(fields)
        axis = fix_axis(axis, ds)
        center, dcenter = ds.coordinates.sanitize_center(center, axis)
        slc = ds.slice(axis, center[axis], **kwargs)
        w, frb = construct_image(ds, axis, slc, dcenter, width=width, 
                                 image_res=image_res)
        super(FITSSlice, self).__init__(frb, fields=fields, wcs=w)


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
    weight_field : string
        The field used to weight the projection.
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
         Width can have four different formats to support windows with variable
         x and y widths.  They are:

         ==================================     =======================
         format                                 example
         ==================================     =======================
         (float, string)                        (10,'kpc')
         ((float, string), (float, string))     ((10,'kpc'),(15,'kpc'))
         float                                  0.2
         (float, float)                         (0.2, 0.3)
         ==================================     =======================

         For example, (10, 'kpc') requests a plot window that is 10 kiloparsecs
         wide in the x and y directions, ((10,'kpc'),(15,'kpc')) requests a
         window that is 10 kiloparsecs wide along the x axis and 15
         kiloparsecs wide along the y axis.  In the other two examples, code
         units are assumed, for example (0.2, 0.3) requests a plot that has an
         x width of 0.2 and a y width of 0.3 in code units.  If units are
         provided the resulting plot axis labels will use the supplied units.
    image_res : an int or 2-tuple of ints
        Specify the resolution of the resulting image. If not provided, it will
        be determined based on the minimum cell size of the dataset.
    """
    def __init__(self, ds, axis, fields, center="c", width=None,
                 weight_field=None, image_res=None, **kwargs):
        fields = ensure_list(fields)
        axis = fix_axis(axis, ds)
        center, dcenter = ds.coordinates.sanitize_center(center, axis)
        prj = ds.proj(fields[0], axis, weight_field=weight_field, **kwargs)
        w, frb = construct_image(ds, axis, prj, dcenter, width=width, 
                                 image_res=image_res)
        super(FITSProjection, self).__init__(frb, fields=fields, wcs=w)

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
    center : A sequence of floats, a string, or a tuple.
        The coordinate of the center of the image. If set to 'c', 'center' or
        left blank, the plot is centered on the middle of the domain. If set to
        'max' or 'm', the center will be located at the maximum of the
        ('gas', 'density') field. Centering on the max or min of a specific
        field is supported by providing a tuple such as ("min","temperature") 
        or ("max","dark_matter_density"). Units can be specified by passing in 
        *center* as a tuple containing a coordinate and string unit name or by 
        passing in a YTArray. If a list or unitless array is supplied, code 
        units are assumed.
    width : tuple or a float.
        Width can have four different formats to support windows with variable
        x and y widths.  They are:

        ==================================     =======================
        format                                 example
        ==================================     =======================
        (float, string)                        (10,'kpc')
        ((float, string), (float, string))     ((10,'kpc'),(15,'kpc'))
        float                                  0.2
        (float, float)                         (0.2, 0.3)
        ==================================     =======================

        For example, (10, 'kpc') requests a plot window that is 10 kiloparsecs
        wide in the x and y directions, ((10,'kpc'),(15,'kpc')) requests a
        window that is 10 kiloparsecs wide along the x axis and 15
        kiloparsecs wide along the y axis.  In the other two examples, code
        units are assumed, for example (0.2, 0.3) requests a plot that has an
        x width of 0.2 and a y width of 0.3 in code units.  If units are
        provided the resulting plot axis labels will use the supplied units.
    image_res : an int or 2-tuple of ints
        Specify the resolution of the resulting image.
    north_vector : a sequence of floats
        A vector defining the 'up' direction in the plot.  This
        option sets the orientation of the slicing plane.  If not
        set, an arbitrary grid-aligned north-vector is chosen.
    """
    def __init__(self, ds, normal, fields, center='c', width=None,
                 image_res=512, north_vector=None):
        fields = ensure_list(fields)
        center, dcenter = ds.coordinates.sanitize_center(center, 4)
        cut = ds.cutting(normal, center, north_vector=north_vector)
        center = ds.arr([0.0] * 2, 'code_length')
        w, frb = construct_image(ds, normal, cut, center, width=width,
                                 image_res=image_res)
        super(FITSOffAxisSlice, self).__init__(frb, fields=fields, wcs=w)


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
         Width can have four different formats to support windows with variable
         x and y widths.  They are:

         ==================================     =======================
         format                                 example
         ==================================     =======================
         (float, string)                        (10,'kpc')
         ((float, string), (float, string))     ((10,'kpc'),(15,'kpc'))
         float                                  0.2
         (float, float)                         (0.2, 0.3)
         ==================================     =======================

         For example, (10, 'kpc') requests a plot window that is 10 kiloparsecs
         wide in the x and y directions, ((10,'kpc'),(15,'kpc')) requests a
         window that is 10 kiloparsecs wide along the x axis and 15
         kiloparsecs wide along the y axis.  In the other two examples, code
         units are assumed, for example (0.2, 0.3) requests a plot that has an
         x width of 0.2 and a y width of 0.3 in code units.  If units are
         provided the resulting plot axis labels will use the supplied units.
    depth : A tuple or a float
         A tuple containing the depth to project through and the string
         key of the unit: (width, 'unit').  If set to a float, code units
         are assumed
    weight_field : string
         The name of the weighting field.  Set to None for no weight.
    image_res : an int or 2-tuple of ints
        Specify the resolution of the resulting image. 
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
    """
    def __init__(self, ds, normal, fields, center='c', width=(1.0, 'unitary'),
                 weight_field=None, image_res=512, data_source=None,
                 north_vector=None, depth=(1.0, "unitary"), method='integrate'):
        fields = ensure_list(fields)
        center, dcenter = ds.coordinates.sanitize_center(center, 4)
        buf = {}
        width = ds.coordinates.sanitize_width(normal, width, depth)
        wd = tuple(el.in_units('code_length').v for el in width)
        if not iterable(image_res):
            image_res = (image_res, image_res)
        res = (image_res[0], image_res[1])
        if data_source is None:
            source = ds
        else:
            source = data_source
        for field in fields:
            buf[field] = off_axis_projection(source, center, normal, wd,
                                             res, field, north_vector=north_vector,
                                             method=method, weight=weight_field).swapaxes(0,1)
        center = ds.arr([0.0] * 2, 'code_length')
        w, not_an_frb = construct_image(ds, normal, buf, center, width=width, image_res=image_res)
        super(FITSOffAxisProjection, self).__init__(buf, fields=fields, wcs=w)
