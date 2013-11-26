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
from yt.funcs import mylog, iterable
from yt.visualization.fixed_resolution import FixedResolutionBuffer
from yt.data_objects.construction_data_containers import YTCoveringGridBase

try:
    from astropy.io.fits import HDUList, ImageHDU
    from astropy import wcs as pywcs
except ImportError:
    ImportError("You need to have AstroPy installed.")

class FITSImageBuffer(HDUList):

    def __init__(self, data, fields=None, center=None, units="cm", wcs=None, D_A=None):
        r""" Initialize a FITSImageBuffer object.

        FITSImageBuffer contains a list of FITS ImageHDU instances, and optionally includes
        WCS information. It inherits from HDUList, so operations such as `writeto` are
        enabled. Images can be constructed from ImageArrays, NumPy arrays, dicts of such
        arrays, or FixedResolutionBuffers. The latter is the most powerful because WCS
        coordinate information can be constructed from the buffer's coordinates. 

        Parameters
        ----------
        data : FixedResolutionBuffer, or ImageArray, numpy.ndarray, or dict of such arrays
            The data to be made into a FITS image or images.
        fields : single string or list of strings, optional
            The field names for the data. If *fields* is none and *data* has keys,
            it will use these for the fields. If *data* is just a single array one field name
            must be specified.
        center : array_like, optional
            The coordinates [xctr,yctr] of the images in units *units*. If *wcs* is set
            this is ignored.
        units : string, optional
            The units of the WCS coordinates. If *wcs* is set this is ignored.
        wcs : astropy.wcs object, optional
            A WCS object. If set, it overrides the other coordinate-related keyword arguments.
        D_A : float or (value, unit) pair, optional
            Angular diameter distance. Only used when constructing WCS information
            for a FixedResolutionBuffer and *units* == "deg". 

        Examples
        --------

        >>> ds = load("sloshing_nomag2_hdf5_plt_cnt_0150")
        >>> prj = ds.h.proj(2, "TempkeV", weight_field="Density")
        >>> frb = prj.to_frb((0.5, "mpc"), 800)
        >>> f_kpc = FITSImageBuffer(frb, fields="TempkeV", units="kpc")
        >>> D_A = 300.0 # in Mpc
        >>> f_deg = FITSImageBuffer(frb, fields="TempkeV", units="deg",
                                    D_A=(D_A, "mpc"), center=(30., 45.))
        >>> f_deg.writeto("temp.fits")
        """
        
        HDUList.__init__(self)

        if isinstance(fields, basestring): fields = [fields]
            
        exclude_fields = ['x','y','z','px','py','pz',
                          'pdx','pdy','pdz','weight_field']

        if center is not None: center = np.array(center)
        
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
        
        for key in fields:
            if key not in exclude_fields:
                mylog.info("Making a FITS image of field %s" % (key))
                hdu = ImageHDU(np.array(img_data[key]).transpose(), name=key)
                self.append(hdu)

        self.dimensionality = len(self[0].data.shape)
        
        if self.dimensionality == 2:
            self.nx, self.ny = self[0].data.shape
        elif self.dimensionality == 3:
            self.nx, self.ny, self.nz = self[0].data.shape

        self.axes = [None]*self.dimensionality

        if wcs is not None:
            w = wcs
        elif isinstance(img_data, FixedResolutionBuffer):
            # FRBs are a special case where we have coordinate
            # information, so we take advantage of this and
            # construct the WCS object
            dx = (img_data.bounds[1]-img_data.bounds[0])/self.nx
            dy = (img_data.bounds[3]-img_data.bounds[2])/self.ny
            if units == "deg":
                if center == None:
                    mylog.error("Please specify center=(RA, Dec)")
                    raise ValueError
                if D_A == None:
                    mylog.error("Please specify D_A.")
                    raise ValueError
                if iterable(D_A):
                    D_A = D_A[0]/img_data.pf.units[D_A[1]]
                dx /= -D_A
                dy /= D_A
                dx = np.rad2deg(dx)
                dy = np.rad2deg(dy)
                xctr = center[0]
                yctr = center[1]
                proj_type = ["RA---TAN","DEC--TAN"]
            else:
                dx *= img_data.pf.units[units]
                dy *= img_data.pf.units[units]
                xctr = 0.5*(img_data.bounds[1]+img_data.bounds[0])
                yctr = 0.5*(img_data.bounds[3]+img_data.bounds[2])
                xctr *= img_data.pf.units[units]
                yctr *= img_data.pf.units[units]
                proj_type = ["linear"]*2
            w = pywcs.WCS(header=self[0].header, naxis=2)
            w.wcs.crpix = [0.5*(self.nx+1), 0.5*(self.ny+1)]
            w.wcs.cdelt = [dx,dy]
            w.wcs.crval = [xctr,yctr]
            w.wcs.ctype = proj_type
            w.wcs.cunit = [units]*2
        elif isinstance(img_data, YTCoveringGridBase):
            dx, dy, dz = img_data.dds
            dx *= img_data.pf.units[units]
            dy *= img_data.pf.units[units]
            dz *= img_data.pf.units[units]
            center = 0.5*(img_data.left_edge+img_data.right_edge)
            center *= img_data.pf.units[units]
            w = pywcs.WCS(header=self[0].header, naxis=3)
            w.wcs.crpix = [0.5*(self.nx+1), 0.5*(self.ny+1), 0.5*(self.ny+1)]
            w.wcs.cdelt = [dx,dy,dz]
            w.wcs.crval = center
            w.wcs.ctype = ["linear"]*3
            w.wcs.cunit = [units]*3
        else:
            w = None

        if w is None:
            self.wcs = None
        else:
            self.set_wcs(w)
            
    def set_wcs(self, wcs):
        """
        Set the WCS coordinate information for all images
        with a WCS object *wcs*.
        """
        self.wcs = wcs
        h = self.wcs.to_header()
        for img in self:
            for k, v in h.items():
                img.header.update(k,v)

    def update_all_headers(self, key, value):
        """
        Update the FITS headers for all images with the
        same *key*, *value* pair.
        """
        for img in self: img.header.update(key,value)
            
    def keys(self):
        return [f.name for f in self]

    def has_key(self, key):
        return key in self.keys()

    def values(self):
        return [self[k] for k in self.keys()]

    def items(self):
        return [(k, self[k]) for k in self.keys()]

    def __add__(self, other):
        if len(set(self.keys()).intersection(set(other.keys()))) > 0:
            mylog.error("There are duplicate extension names! Don't know which ones you want to keep!")
            raise KeyError
        new_buffer = {}
        for im1 in self:
            new_buffer[im1.name] = im1.data
        for im2 in other:
            new_buffer[im2.name] = im2.data
        new_wcs = self.wcs
        return FITSImageBuffer(new_buffer, wcs=new_wcs)
    
    @property
    def shape(self):
        if self.dimensionality == 2:
            return self.nx, self.ny
        elif self.dimensionality == 3:
            return self.nx, self.ny, self.nz

    
