"""
Fixed resolution buffer support, along with a primitive image analysis tool.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008-2011 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from yt.funcs import *
from yt.utilities.definitions import \
    x_dict, \
    y_dict, \
    axis_names
from .volume_rendering.api import off_axis_projection
from yt.data_objects.image_array import ImageArray
import _MPL
import numpy as np
import weakref

class FixedResolutionBuffer(object):
    _exclude_fields = ('pz','pdz','dx','x','y','z')
    def __init__(self, data_source, bounds, buff_size, antialias = True,
                 periodic = False):
        r"""
        FixedResolutionBuffer(data_source, bounds, buff_size, antialias = True)

        This accepts a 2D data object, such as a Projection or Slice, and
        implements a protocol for generating a pixelized, fixed-resolution
        image buffer.

        yt stores 2D AMR data internally as a set of 2D coordinates and the
        half-width of individual pixels.  Converting this to an image buffer
        requires a deposition step, where individual variable-resolution pixels
        are deposited into a buffer of some resolution, to create an image.
        This object is an interface to that pixelization step: it can deposit
        multiple fields.  It acts as a standard AMRData object, such that
        dict-style access returns an image of a given field.

        Parameters
        ----------
        data_source : :class:`yt.data_objects.data_containers.AMRProjBase` or :class:`yt.data_objects.data_containers.AMRSliceBase`
            This is the source to be pixelized, which can be a projection or a
            slice.  (For cutting planes, see
            `yt.visualization.fixed_resolution.ObliqueFixedResolutionBuffer`.)
        bounds : sequence of floats
            Bounds are the min and max in the image plane that we want our
            image to cover.  It's in the order of (xmin, xmax, ymin, ymax),
            where the coordinates are all in the appropriate code units.
        buff_size : sequence of ints
            The size of the image to generate.
        antialias : boolean
            This can be true or false.  It determines whether or not sub-pixel
            rendering is used during data deposition.
        periodic : boolean
            This can be true or false, and governs whether the pixelization
            will span the domain boundaries.

        See Also
        --------
        :class:`yt.visualization.fixed_resolution.ObliqueFixedResolutionBuffer` : A similar object,
                                                         used for cutting
                                                         planes.

        Examples
        --------
        To make a projection and then several images, you can generate a
        single FRB and then access multiple fields:

        >>> proj = pf.h.proj(0, "Density")
        >>> frb1 = FixedResolutionBuffer(proj, (0.2, 0.3, 0.4, 0.5),
                        (1024, 1024))
        >>> print frb1["Density"].max()
        1.0914e-9
        >>> print frb1["Temperature"].max()
        104923.1
        """
        self.data_source = data_source
        self.pf = data_source.pf
        self.bounds = bounds
        self.buff_size = buff_size
        self.antialias = antialias
        self.data = {}
        self.axis = data_source.axis
        self.periodic = periodic

        h = getattr(data_source, "hierarchy", None)
        if h is not None:
            h.plots.append(weakref.proxy(self))

        # Handle periodicity, just in case
        if self.data_source.axis < 3:
            DLE = self.pf.domain_left_edge
            DRE = self.pf.domain_right_edge
            DD = float(self.periodic)*(DRE - DLE)
            axis = self.data_source.axis
            xax = x_dict[axis]
            yax = y_dict[axis]
            self._period = (DD[xax], DD[yax])
            self._edges = ( (DLE[xax], DRE[xax]), (DLE[yax], DRE[yax]) )
        
    def keys(self):
        return self.data.keys()
    
    def __delitem__(self, item):
        del self.data[item]
    
    def __getitem__(self, item):
        if item in self.data: return self.data[item]
        mylog.info("Making a fixed resolution buffer of (%s) %d by %d" % \
            (item, self.buff_size[0], self.buff_size[1]))
        buff = _MPL.Pixelize(self.data_source['px'],
                             self.data_source['py'],
                             self.data_source['pdx'],
                             self.data_source['pdy'],
                             self.data_source[item],
                             self.buff_size[0], self.buff_size[1],
                             self.bounds, int(self.antialias),
                             self._period, int(self.periodic),
                             ).transpose()
        ia = ImageArray(buff, info=self._get_info(item))
        self[item] = ia
        return ia 

    def __setitem__(self, item, val):
        self.data[item] = val

    def _get_data_source_fields(self):
        exclude = self.data_source._key_fields + list(self._exclude_fields)
        for f in self.data_source.fields:
            if f not in exclude:
                self[f]

    def _get_info(self, item):
        info = {}
        info['data_source'] = self.data_source.__str__()  
        info['axis'] = self.data_source.axis
        info['field'] = str(item)
        info['units'] = self.data_source.pf.field_info[item].get_units()
        info['xlim'] = self.bounds[:2]
        info['ylim'] = self.bounds[2:]
        info['length_to_cm'] = self.data_source.pf['cm']
        info['projected_units'] = \
                self.data_source.pf.field_info[item].get_projected_units()
        info['center'] = self.data_source.center
        
        try:
            info['coord'] = self.data_source.coord
        except AttributeError:
            pass
        
        try:
            info['weight_field'] = self.data_source.weight_field
        except AttributeError:
            pass
        
        info['label'] = self.data_source.pf.field_info[item].display_name
        if info['label'] is None:
            info['label'] = r'$\rm{'+item+r'}$'
        elif info['label'].find('$') == -1:
            info['label'] = r'$\rm{'+info['label']+r'}$'
        if info['units'] is None or info['units'] == '':
            pass
        else:
            info['label'] += r'$\/\/('+info['units']+r')$'
        
        return info

    def convert_to_pixel(self, coords):
        r"""This function converts coordinates in code-space to pixel-space.

        Parameters
        ----------
        coords : sequence of array_like
            This is (x_coord, y_coord).  Because of the way the math is done,
            these can both be arrays.

        Returns
        -------
        output : sequence of array_like
            This returns px_coord, py_coord

        """
        dpx = (self.bounds[1]-self.bounds[0])/self.buff_size[0]
        dpy = (self.bounds[3]-self.bounds[2])/self.buff_size[1]
        px = (coords[0] - self.bounds[0])/dpx
        py = (coords[1] - self.bounds[2])/dpy
        return (px, py)

    def convert_distance_x(self, distance):
        r"""This function converts code-space distance into pixel-space
        distance in the x-coordiante.

        Parameters
        ----------
        distance : array_like
            This is x-distance in code-space you would like to convert.

        Returns
        -------
        output : array_like
            The return value is the distance in the y-pixel coordinates.

        """
        dpx = (self.bounds[1]-self.bounds[0])/self.buff_size[0]
        return distance/dpx
        
    def convert_distance_y(self, distance):
        r"""This function converts code-space distance into pixel-space
        distance in the y-coordiante.

        Parameters
        ----------
        distance : array_like
            This is y-distance in code-space you would like to convert.

        Returns
        -------
        output : array_like
            The return value is the distance in the x-pixel coordinates.

        """
        dpy = (self.bounds[3]-self.bounds[2])/self.buff_size[1]
        return distance/dpy

    def export_hdf5(self, filename, fields = None):
        r"""Export a set of fields to a set of HDF5 datasets.

        This function will export any number of fields into datasets in a new
        HDF5 file.
        
        Parameters
        ----------
        filename : string
            This file will be opened in "append" mode.
        fields : list of strings
            These fields will be pixelized and output.
        """
        import h5py
        if fields is None: fields = self.data.keys()
        output = h5py.File(filename, "a")
        for field in fields:
            output.create_dataset(field,data=self[field])
        output.close()

    def export_fits(self, filename_prefix, fields = None, clobber=False,
                    other_keys=None, gzip_file=False, units="1"):

        """
        This will export a set of FITS images of either the fields specified
        or all the fields already in the object.  The output filename is
        *filename_prefix*. If clobber is set to True, this will overwrite any
        existing FITS file.

        This requires the *pyfits* module, which is a standalone module
        provided by STSci to interface with FITS-format files.
        """
        r"""Export a set of pixelized fields to a FITS file.

        This will export a set of FITS images of either the fields specified
        or all the fields already in the object.  The output filename is the
        the specified prefix.

        Parameters
        ----------
        filename_prefix : string
            This prefix will be prepended to every FITS file name.
        fields : list of strings
            These fields will be pixelized and output.
        clobber : boolean
            If the file exists, this governs whether we will overwrite.
        other_keys : dictionary, optional
            A set of header keys and values to write into the FITS header.
        gzip_file : boolean, optional
            gzip the file after writing, default False
        units : string, optional
            the length units that the coordinates are written in, default '1'
        """
        
        import pyfits
        from os import system
        
        extra_fields = ['x','y','z','px','py','pz','pdx','pdy','pdz','weight_field']
        if filename_prefix.endswith('.fits'): filename_prefix=filename_prefix[:-5]
        if fields is None: 
            fields = [field for field in self.data_source.fields 
                      if field not in extra_fields]

        nx, ny = self.buff_size
        dx = (self.bounds[1]-self.bounds[0])/nx*self.pf[units]
        dy = (self.bounds[3]-self.bounds[2])/ny*self.pf[units]
        xmin = self.bounds[0]*self.pf[units]
        ymin = self.bounds[2]*self.pf[units]
        simtime = self.pf.current_time

        hdus = []

        first = True
        
        for field in fields:

            if (first) :
                hdu = pyfits.PrimaryHDU(self[field])
                first = False
            else :
                hdu = pyfits.ImageHDU(self[field])
                
            if self.data_source.has_key('weight_field'):
                weightname = self.data_source._weight
                if weightname is None: weightname = 'None'
                field = field +'_'+weightname

            hdu.header.update("Field", field)
            hdu.header.update("Time", simtime)

            hdu.header.update('WCSNAMEP', "PHYSICAL")            
            hdu.header.update('CTYPE1P', "LINEAR")
            hdu.header.update('CTYPE2P', "LINEAR")
            hdu.header.update('CRPIX1P', 0.5)
            hdu.header.update('CRPIX2P', 0.5)
            hdu.header.update('CRVAL1P', xmin)
            hdu.header.update('CRVAL2P', ymin)
            hdu.header.update('CDELT1P', dx)
            hdu.header.update('CDELT2P', dy)
                    
            hdu.header.update('CTYPE1', "LINEAR")
            hdu.header.update('CTYPE2', "LINEAR")                                
            hdu.header.update('CUNIT1', units)
            hdu.header.update('CUNIT2', units)
            hdu.header.update('CRPIX1', 0.5)
            hdu.header.update('CRPIX2', 0.5)
            hdu.header.update('CRVAL1', xmin)
            hdu.header.update('CRVAL2', ymin)
            hdu.header.update('CDELT1', dx)
            hdu.header.update('CDELT2', dy)

            if (other_keys is not None) :

                for k,v in other_keys.items() :

                    hdu.header.update(k,v)

            hdus.append(hdu)

            del hdu
            
        hdulist = pyfits.HDUList(hdus)

        hdulist.writeto("%s.fits" % (filename_prefix), clobber=clobber)
        
        if (gzip_file) :
            clob = ""
            if (clobber) : clob = "-f"
            system("gzip "+clob+" %s.fits" % (filename_prefix))
        
    def open_in_ds9(self, field, take_log=True):
        """
        This will open a given field in the DS9 viewer.

        Displaying fields can often be much easier in an interactive viewer,
        particularly one as versatile as DS9.  This function will pixelize a
        field and export it to an interactive DS9 package.  This requires the
        *numdisplay* package, which is a simple download from STSci.
        Furthermore, it presupposed that it can connect to DS9 -- that is, that
        DS9 is already open.

        Parameters
        ----------
        field : strings
            This field will be pixelized and displayed.
        take_log : boolean
            DS9 seems to have issues with logging fields in-memory.  This will
            pre-log the field before sending it to DS9.
        """
        import numdisplay
        numdisplay.open()
        if take_log: data=np.log10(self[field])
        else: data=self[field]
        numdisplay.display(data)    

    @property
    def limits(self):
        rv = dict(x = None, y = None, z = None)
        xax = x_dict[self.axis]
        yax = y_dict[self.axis]
        xn = axis_names[xax]
        yn = axis_names[yax]
        rv[xn] = (self.bounds[0], self.bounds[1])
        rv[yn] = (self.bounds[2], self.bounds[3])
        return rv

class ObliqueFixedResolutionBuffer(FixedResolutionBuffer):
    """
    This object is a subclass of :class:`yt.visualization.fixed_resolution.FixedResolutionBuffer`
    that supports non-aligned input data objects, primarily cutting planes.
    """
    def __getitem__(self, item):
        if item in self.data: return self.data[item]
        indices = np.argsort(self.data_source['dx'])[::-1]
        buff = _MPL.CPixelize( self.data_source['x'],   self.data_source['y'],   self.data_source['z'],
                               self.data_source['px'],  self.data_source['py'],
                               self.data_source['pdx'], self.data_source['pdy'], self.data_source['pdz'],
                               self.data_source.center, self.data_source._inv_mat, indices,
                               self.data_source[item],
                               self.buff_size[0], self.buff_size[1],
                               self.bounds).transpose()
        ia = ImageArray(buff, info=self._get_info(item))
        self[item] = ia
        return ia 


class OffAxisProjectionFixedResolutionBuffer(FixedResolutionBuffer):
    def __init__(self, data_source, bounds, buff_size, antialias = True,                                                         
                 periodic = False):
        self.data = {}
        FixedResolutionBuffer.__init__(self, data_source, bounds, buff_size, antialias, periodic)

    def __getitem__(self, item):
        if item in self.data: return self.data[item]
        mylog.info("Making a fixed resolutuion buffer of (%s) %d by %d" % \
            (item, self.buff_size[0], self.buff_size[1]))
        ds = self.data_source
        width = (self.bounds[1] - self.bounds[0],
                 self.bounds[3] - self.bounds[2],
                 self.bounds[5] - self.bounds[4])
        buff = off_axis_projection(ds.pf, ds.center, ds.normal_vector,
                                   width, ds.resolution, item,
                                   weight=ds.weight_field, volume=ds.volume,
                                   no_ghost=ds.no_ghost, interpolated=ds.interpolated,
                                   north_vector=ds.north_vector)
        ia = ImageArray(buff.swapaxes(0,1), info=self._get_info(item))
        self[item] = ia
        return ia 


