"""
FITSImageBuffer Class
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from yt.funcs import mylog, iterable, fix_axis, ensure_list
from yt.visualization.fixed_resolution import FixedResolutionBuffer
from yt.data_objects.construction_data_containers import YTCoveringGridBase
from yt.utilities.on_demand_imports import _astropy
from yt.units.yt_array import YTQuantity
import re

pyfits = _astropy.pyfits
pywcs = _astropy.pywcs

if pyfits is None:
    HDUList = object
else:
    HDUList = pyfits.HDUList

class FITSImageBuffer(HDUList):

    def __init__(self, data, fields=None, units="cm", wcs=None):
        r""" Initialize a FITSImageBuffer object.

        FITSImageBuffer contains a list of FITS ImageHDU instances, and
        optionally includes WCS information. It inherits from HDUList, so
        operations such as `writeto` are enabled. Images can be constructed
        from ImageArrays, NumPy arrays, dicts of such arrays,
        FixedResolutionBuffers, and YTCoveringGrids. The latter two are the
        most powerful because WCS information can be constructed from their
        coordinates.

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
            The units of the WCS coordinates. Only applies
            to FixedResolutionBuffers or YTCoveringGrids. Defaults to "cm".
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
        >>> f_kpc = FITSImageBuffer(frb, fields="kT", units="kpc")
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
        >>> f_deg = FITSImageBuffer(frb, fields="kT", wcs=w)
        >>> f_deg.writeto("temp.fits")
        """
        
        super(FITSImageBuffer, self).__init__()

        if isinstance(fields, basestring): fields = [fields]
            
        exclude_fields = ['x','y','z','px','py','pz',
                          'pdx','pdy','pdz','weight_field']
        
        if hasattr(data, 'keys'):
            img_data = data
        else:
            img_data = {}
            if fields is None:
                mylog.error("Please specify a field name for this array.")
                raise KeyError
            img_data[fields[0]] = data

        if fields is None: fields = img_data.keys()
        if len(fields) == 0:
            mylog.error("Please specify one or more fields to write.")
            raise KeyError

        first = True

        for key in fields:
            if key not in exclude_fields:
                mylog.info("Making a FITS image of field %s" % (key))
                if first:
                    hdu = pyfits.PrimaryHDU(np.array(img_data[key]))
                    first = False
                else:
                    hdu = pyfits.ImageHDU(np.array(img_data[key]))
                hdu.name = key
                hdu.header["btype"] = key
                if hasattr(img_data[key], "units"):
                    hdu.header["bunit"] = re.sub('()', '', str(img_data[key].units))
                self.append(hdu)

        self.dimensionality = len(self[0].data.shape)
        
        if self.dimensionality == 2:
            self.nx, self.ny = self[0].data.shape
        elif self.dimensionality == 3:
            self.nx, self.ny, self.nz = self[0].data.shape

        has_coords = (isinstance(img_data, FixedResolutionBuffer) or
                      isinstance(img_data, YTCoveringGridBase))

        if wcs is None:
            w = pywcs.WCS(header=self[0].header, naxis=self.dimensionality)
            w.wcs.crpix = 0.5*(np.array(self.shape)+1)
            if isinstance(img_data, FixedResolutionBuffer):
                # FRBs are a special case where we have coordinate
                # information, so we take advantage of this and
                # construct the WCS object
                dx = (img_data.bounds[1]-img_data.bounds[0]).in_units(units)/self.nx
                dy = (img_data.bounds[3]-img_data.bounds[2]).in_units(units)/self.ny
                xctr = 0.5*(img_data.bounds[1]+img_data.bounds[0]).in_units(units)
                yctr = 0.5*(img_data.bounds[3]+img_data.bounds[2]).in_units(units)
                center = [xctr, yctr]
            elif isinstance(img_data, YTCoveringGridBase):
                dx, dy, dz = img_data.dds.in_units(units)
                center = 0.5*(img_data.left_edge+img_data.right_edge).in_units(units)
            else:
                # We default to pixel coordinates if nothing is provided
                dx, dy, dz = 1.0, 1.0, 1.0
                center = 0.5*(np.array(self.shape)+1)
            w.wcs.crval = center
            if has_coords:
                w.wcs.cunit = [units]*self.dimensionality
            if self.dimensionality == 2:
                w.wcs.cdelt = [dx,dy]
            elif self.dimensionality == 3:
                w.wcs.cdelt = [dx,dy,dz]
            w.wcs.ctype = ["linear"]*self.dimensionality
            self._set_wcs(w)
        else:
            self._set_wcs(wcs)

    def _set_wcs(self, wcs):
        """
        Set the WCS coordinate information for all images
        with a WCS object *wcs*.
        """
        self.wcs = wcs
        h = self.wcs.to_header()
        for img in self:
            for k, v in h.items():
                img.header[k] = v

    def update_all_headers(self, key, value):
        """
        Update the FITS headers for all images with the
        same *key*, *value* pair.
        """
        for img in self: img.header[key] = value
            
    def keys(self):
        return [f.name for f in self]

    def has_key(self, key):
        return key in self.keys()

    def values(self):
        return [self[k] for k in self.keys()]

    def items(self):
        return [(k, self[k]) for k in self.keys()]

    def writeto(self, fileobj, **kwargs):
        pyfits.HDUList(self).writeto(fileobj, **kwargs)

    @property
    def shape(self):
        if self.dimensionality == 2:
            return self.nx, self.ny
        elif self.dimensionality == 3:
            return self.nx, self.ny, self.nz

    def to_glue(self, label="yt", data_collection=None):
        """
        Takes the data in the FITSImageBuffer and exports it to
        Glue (http://www.glueviz.org) for interactive
        analysis. Optionally add a *label*. If you are already within
        the Glue environment, you can pass a *data_collection* object,
        otherwise Glue will be started.
        """
        from glue.core import DataCollection, Data
        from glue.core.coordinates import coordinates_from_header
        from glue.qt.glue_application import GlueApplication

        field_dict = dict((key,self[key].data) for key in self.keys())
        
        image = Data(label=label)
        image.coords = coordinates_from_header(self.wcs.to_header())
        for k,v in field_dict.items():
            image.add_component(v, k)
        if data_collection is None:
            dc = DataCollection([image])
            app = GlueApplication(dc)
            app.start()
        else:
            data_collection.append(image)

    def to_aplpy(self, **kwargs):
        """
        Use APLpy (http://aplpy.github.io) for plotting. Returns an
        `aplpy.FITSFigure` instance. All keyword arguments are passed to the
        `aplpy.FITSFigure` constructor.
        """
        import aplpy
        return aplpy.FITSFigure(self, **kwargs)

axis_wcs = [[1,2],[0,2],[0,1]]

def sanitize_fits_unit(unit, dx):
    if unit == "Mpc":
        mylog.info("Changing FITS file unit to kpc.")
        unit = "kpc"
        dx *= 1000.
    elif unit == "au":
        unit = "AU"
    return unit, dx

def construct_image(data_source, center=None, width=None):
    ds = data_source.ds
    axis = data_source.axis
    if center is None or width is None:
        center = ds.domain_center[axis_wcs[axis]]
    if width is None:
        width = ds.domain_width[axis_wcs[axis]]
        mylog.info("Making an image of the entire domain, "+
                   "so setting the center to the domain center.")
    else:
        width = ds.coordinates.sanitize_width(axis, width, None)
    dx = ds.index.get_smallest_dx()
    nx, ny = [int((w/dx).in_units("dimensionless")) for w in width]
    crpix = [0.5*(nx+1), 0.5*(ny+1)]
    if hasattr(ds, "wcs"):
        # This is a FITS dataset, so we use it to construct the WCS
        cunit = [str(ds.wcs.wcs.cunit[idx]) for idx in axis_wcs[axis]]
        ctype = [ds.wcs.wcs.ctype[idx] for idx in axis_wcs[axis]]
        crval = [ds.wcs.wcs.crval[idx] for idx in axis_wcs[axis]]
        cdelt = [ds.wcs.wcs.cdelt[idx] for idx in axis_wcs[axis]]
    else:
        # This is some other kind of dataset                                                                                    
        unit = str(width[0].units)
        if unit == "unitary":
            unit = ds.get_smallest_appropriate_unit(ds.domain_width.max())
        unit = sanitize_fits_unit(unit)
        cunit = [unit]*2
        ctype = ["LINEAR"]*2
        crval = [center[idx].in_units(unit) for idx in axis_wcs[axis]]
        cdelt = [dx.in_units(unit)]*2
    frb = data_source.to_frb(width[0], (nx,ny), center=center, height=width[1])
    w = pywcs.WCS(naxis=2)
    w.wcs.crpix = crpix
    w.wcs.cdelt = cdelt
    w.wcs.crval = crval
    w.wcs.cunit = cunit
    w.wcs.ctype = ctype
    return w, frb

class FITSSlice(FITSImageBuffer):
    r"""
    Generate a FITSImageBuffer of an on-axis slice.

    Parameters
    ----------
    ds : FITSDataset
        The FITS dataset object.
    axis : character or integer
        The axis of the slice. One of "x","y","z", or 0,1,2.
    fields : string or list of strings
        The fields to slice
    center : A sequence floats, a string, or a tuple.
         The coordinate of the origin of the image. If set to 'c', 'center' or
         left blank, the plot is centered on the middle of the domain. If set to
         'max' or 'm', the center will be located at the maximum of the
         ('gas', 'density') field. Units can be specified by passing in center
         as a tuple containing a coordinate and string unit name or by passing
         in a YTArray.  If a list or unitless array is supplied, code units are
         assumed.
    width : 
    """
    def __init__(self, ds, axis, fields, center="c", width=None, **kwargs):
        fields = ensure_list(fields)
        axis = fix_axis(axis, ds)
        center = ds.coordinates.sanitize_center(center, axis)
        slc = ds.slice(axis, center[axis], **kwargs)
        w, frb = construct_image(slc, center=dcenter, width=width)
        super(FITSSlice, self).__init__(frb, fields=fields, wcs=w)
        for i, field in enumerate(fields):
            self[i].header["bunit"] = str(frb[field].units)

class FITSProjection(FITSImageBuffer):
    r"""
    Generate a FITSImageBuffer of an on-axis projection.

    Parameters
    ----------
    ds : FITSDataset
        The FITS dataset object.
    axis : character or integer
        The axis along which to project. One of "x","y","z", or 0,1,2.
    fields : string or list of strings
        The fields to project
    weight_field : string
        The field used to weight the projection.
    center : A sequence floats, a string, or a tuple.
        The coordinate of the origin of the image. If set to 'c', 'center' or
        left blank, the plot is centered on the middle of the domain. If set to
        'max' or 'm', the center will be located at the maximum of the
        ('gas', 'density') field. Units can be specified by passing in center
        as a tuple containing a coordinate and string unit name or by passing
        in a YTArray.  If a list or unitless array is supplied, code units are
        assumed.
    width :
    """
    def __init__(self, ds, axis, fields, center="c", width=None, weight_field=None, **kwargs):
        fields = ensure_list(fields)
        axis = fix_axis(axis, ds)
        center, dcenter = ds.coordinates.sanitize_center(center, axis)
        prj = ds.proj(fields[0], axis, weight_field=weight_field, **kwargs)
        w, frb = construct_image(prj, center=dcenter, width=width)
        super(FITSProjection, self).__init__(frb, fields=fields, wcs=w)
        for i, field in enumerate(fields):
            self[i].header["bunit"] = str(frb[field].units)



        

    
