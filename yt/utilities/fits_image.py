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
from yt.visualization.plot_window import get_sanitized_center
from yt.data_objects.construction_data_containers import YTCoveringGridBase
from yt.utilities.on_demand_imports import _astropy
from yt.units.yt_array import YTQuantity

pyfits = _astropy.pyfits
pywcs = _astropy.pywcs

if pyfits is None:
    HDUList = object
else:
    HDUList = pyfits.HDUList

class FITSImageBuffer(HDUList):

    def __init__(self, data, fields=None, units="cm",
                 center=None, scale=None, wcs=None):
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
            The units of the WCS coordinates, default "cm". 
        center : array_like, optional
            The coordinates [xctr,yctr] of the images in units
            *units*. If *units* is not specified, defaults to the origin. 
        scale : tuple of floats, optional
            Pixel scale in unit *units*. Will be ignored if *data* is
            a FixedResolutionBuffer or a YTCoveringGrid. Must be
            specified otherwise, or if *units* is "deg".
        wcs : `astropy.wcs.WCS` instance, optional
            Supply an AstroPy WCS instance to override automatic WCS creation.

        Examples
        --------

        >>> ds = load("sloshing_nomag2_hdf5_plt_cnt_0150")
        >>> prj = ds.proj(2, "kT", weight_field="density")
        >>> frb = prj.to_frb((0.5, "Mpc"), 800)
        >>> # This example just uses the FRB and puts the coords in kpc.
        >>> f_kpc = FITSImageBuffer(frb, fields="kT", units="kpc")
        >>> # This example specifies sky coordinates.
        >>> scale = [1./3600.]*2 # One arcsec per pixel
        >>> f_deg = FITSImageBuffer(frb, fields="kT", units="deg",
                                    scale=scale, center=(30., 45.))
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
                    hdu.header["bunit"] = str(img_data[key].units)
                self.append(hdu)

        self.dimensionality = len(self[0].data.shape)
        
        if self.dimensionality == 2:
            self.nx, self.ny = self[0].data.shape
        elif self.dimensionality == 3:
            self.nx, self.ny, self.nz = self[0].data.shape

        has_coords = (isinstance(img_data, FixedResolutionBuffer) or
                      isinstance(img_data, YTCoveringGridBase))
        
        if center is None:
            if units == "deg":
                mylog.error("Please specify center=(RA, Dec) in degrees.")
                raise ValueError
            elif not has_coords:
                mylog.warning("Setting center to the origin.")
                center = [0.0]*self.dimensionality

        if scale is None:
            if units == "deg" or not has_coords and wcs is None:
                mylog.error("Please specify scale=(dx,dy[,dz]) in %s." % (units))
                raise ValueError

        if wcs is None:
            w = pywcs.WCS(header=self[0].header, naxis=self.dimensionality)
            w.wcs.crpix = 0.5*(np.array(self.shape)+1)
            proj_type = ["linear"]*self.dimensionality
            if isinstance(img_data, FixedResolutionBuffer) and units != "deg":
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
            elif units == "deg" and self.dimensionality == 2:
                dx = -scale[0]
                dy = scale[1]
                proj_type = ["RA---TAN","DEC--TAN"]
            else:
                dx = scale[0]
                dy = scale[1]
                if self.dimensionality == 3: dz = scale[2]
            
            w.wcs.crval = center
            w.wcs.cunit = [units]*self.dimensionality
            w.wcs.ctype = proj_type
        
            if self.dimensionality == 2:
                w.wcs.cdelt = [dx,dy]
            elif self.dimensionality == 3:
                w.wcs.cdelt = [dx,dy,dz]

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

def construct_image(data_source, center=None):
    ds = data_source.ds
    axis = data_source.axis
    if hasattr(ds, "wcs"):
        # This is a FITS dataset
        nx, ny = ds.domain_dimensions[axis_wcs[axis]]
        crpix = [ds.wcs.wcs.crpix[idx] for idx in axis_wcs[axis]]
        cdelt = [ds.wcs.wcs.cdelt[idx] for idx in axis_wcs[axis]]
        crval = [ds.wcs.wcs.crval[idx] for idx in axis_wcs[axis]]
        cunit = [str(ds.wcs.wcs.cunit[idx]) for idx in axis_wcs[axis]]
        ctype = [ds.wcs.wcs.ctype[idx] for idx in axis_wcs[axis]]
    else:
        # This is some other kind of dataset
        unit = ds.get_smallest_appropriate_unit(ds.domain_width.max())
        if center is None:
            crval = [0.0,0.0]
        else:
            crval = [(ds.domain_center-center)[idx].in_units(unit) for idx in axis_wcs[axis]]
        dx = ds.index.get_smallest_dx()
        nx, ny = (ds.domain_width[axis_wcs[axis]]/dx).ndarray_view().astype("int")
        crpix = [0.5*(nx+1), 0.5*(ny+1)]
        cdelt = [dx.in_units(unit)]*2
        cunit = [unit]*2
        ctype = ["LINEAR"]*2
    frb = data_source.to_frb((1.0,"unitary"), (nx,ny))
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
    """
    def __init__(self, ds, axis, fields, center="c", **kwargs):
        fields = ensure_list(fields)
        axis = fix_axis(axis, ds)
        center = get_sanitized_center(center, ds)
        slc = ds.slice(axis, center[axis], **kwargs)
        w, frb = construct_image(slc, center=center)
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
    """
    def __init__(self, ds, axis, fields, center="c", weight_field=None, **kwargs):
        fields = ensure_list(fields)
        axis = fix_axis(axis, ds)
        center = get_sanitized_center(center, ds)
        prj = ds.proj(fields[0], axis, weight_field=weight_field, **kwargs)
        w, frb = construct_image(prj, center=center)
        super(FITSProjection, self).__init__(frb, fields=fields, wcs=w)
        for i, field in enumerate(fields):
            self[i].header["bunit"] = str(frb[field].units)



        

    
