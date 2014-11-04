"""
Fixed resolution buffer support, along with a primitive image analysis tool.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.funcs import *
from yt.units.unit_object import Unit
from .volume_rendering.api import off_axis_projection
from yt.data_objects.image_array import ImageArray
from yt.utilities.lib.misc_utilities import \
    pixelize_cylinder
from . import _MPL
import numpy as np
import weakref
import re
import string

class FixedResolutionBuffer(object):
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

    >>> proj = ds.proj(0, "Density")
    >>> frb1 = FixedResolutionBuffer(proj, (0.2, 0.3, 0.4, 0.5),
                    (1024, 1024))
    >>> print frb1["Density"].max()
    1.0914e-9
    >>> print frb1["Temperature"].max()
    104923.1
    """
    _exclude_fields = ('pz','pdz','dx','x','y','z',
                       ('index','dx'),('index','x'),('index','y'),('index','z'))
    def __init__(self, data_source, bounds, buff_size, antialias = True,
                 periodic = False):
        self.data_source = data_source
        self.ds = data_source.ds
        self.bounds = bounds
        self.buff_size = buff_size
        self.antialias = antialias
        self.data = {}
        self.axis = data_source.axis
        self.periodic = periodic

        ds = getattr(data_source, "ds", None)
        if ds is not None:
            ds.plots.append(weakref.proxy(self))

        # Handle periodicity, just in case
        if self.data_source.axis < 3:
            DLE = self.ds.domain_left_edge
            DRE = self.ds.domain_right_edge
            DD = float(self.periodic)*(DRE - DLE)
            axis = self.data_source.axis
            xax = self.ds.coordinates.x_axis[axis]
            yax = self.ds.coordinates.y_axis[axis]
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
        bounds = []
        for b in self.bounds:
            if hasattr(b, "in_units"):
                b = float(b.in_units("code_length"))
            bounds.append(b)
        buff = self.ds.coordinates.pixelize(self.data_source.axis,
            self.data_source, item, bounds, self.buff_size,
            int(self.antialias))
        # Need to add _period and self.periodic
        # self._period, int(self.periodic)
        ia = ImageArray(buff, input_units=self.data_source[item].units,
                        info=self._get_info(item))
        self.data[item] = ia
        return self.data[item]

    def __setitem__(self, item, val):
        self.data[item] = val

    def _get_data_source_fields(self):
        exclude = self.data_source._key_fields + list(self._exclude_fields)
        fields = getattr(self.data_source, "fields", [])
        fields += getattr(self.data_source, "field_data", {}).keys()
        for f in fields:
            if f not in exclude and f[0] not in self.data_source.ds.particle_types:
                self[f]


    def _is_ion( self, fname ):
        p = re.compile("_p[0-9]+_")
        result = False
        if p.search( fname ) != None:
            result = True
        return result

    def _ion_to_label( self, fname ):
        pnum2rom = {
            "0":"I", "1":"II", "2":"III", "3":"IV", "4":"V",
            "5":"VI", "6":"VII", "7":"VIII", "8":"IX", "9":"X",
            "10":"XI", "11":"XII", "12":"XIII", "13":"XIV", "14":"XV",
            "15":"XVI", "16":"XVII", "17":"XVIII", "18":"XIX", "19":"XX"}

        p = re.compile("_p[0-9]+_")
        m = p.search( fname )
        if m != None:
            pstr = m.string[m.start()+1:m.end()-1]
            segments = fname.split("_")
            for i,s in enumerate(segments):
                segments[i] = string.capitalize(s)
                if s == pstr:
                    ipstr = i
            element = segments[ipstr-1]
            roman = pnum2rom[pstr[1:]] 
            label = element + '\/' + roman + '\/' + \
                string.join( segments[ipstr+1:], '\/' ) 
        else:
            label = fname
        return label




    def _get_info(self, item):
        info = {}
        ftype, fname = field = self.data_source._determine_fields(item)[0]
        finfo = self.data_source.ds._get_field_info(*field)
        info['data_source'] = self.data_source.__str__()  
        info['axis'] = self.data_source.axis
        info['field'] = str(item)
        info['xlim'] = self.bounds[:2]
        info['ylim'] = self.bounds[2:]
        info['length_unit'] = self.data_source.ds.length_unit
        info['length_to_cm'] = info['length_unit'].in_cgs().to_ndarray()
        info['center'] = self.data_source.center
        
        try:
            info['coord'] = self.data_source.coord
        except AttributeError:
            pass
        
        try:
            info['weight_field'] = self.data_source.weight_field
        except AttributeError:
            pass
        
        info['label'] = finfo.display_name
        if info['label'] is None:
            if self._is_ion( fname ):
                fname = self._ion_to_label( fname )
                info['label'] = r'$\rm{'+fname+r'}$'
                info['label'] = r'$\rm{'+fname.replace('_','\/')+r'}$'
            else:    
                info['label'] = r'$\rm{'+fname+r'}$'
                info['label'] = r'$\rm{'+fname.replace('_','\/').title()+r'}$'
        elif info['label'].find('$') == -1:
            info['label'] = info['label'].replace(' ','\/')
            info['label'] = r'$\rm{'+info['label']+r'}$'
        
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

    def export_fits(self, filename, fields=None, clobber=False,
                    other_keys=None, units="cm"):
        r"""Export a set of pixelized fields to a FITS file.

        This will export a set of FITS images of either the fields specified
        or all the fields already in the object.

        Parameters
        ----------
        filename : string
            The name of the FITS file to be written.
        fields : list of strings
            These fields will be pixelized and output.
        clobber : boolean
            If the file exists, this governs whether we will overwrite.
        other_keys : dictionary, optional
            A set of header keys and values to write into the FITS header.
        units : string, optional
            the length units that the coordinates are written in, default 'cm'
            If units are set to "deg" then assume that sky coordinates are
            requested.
        """

        from yt.utilities.fits_image import FITSImageBuffer

        extra_fields = ['x','y','z','px','py','pz','pdx','pdy','pdz','weight_field']
        if fields is None: 
            fields = [field[-1] for field in self.data_source.field_data
                      if field not in extra_fields]

        fib = FITSImageBuffer(self, fields=fields, units=units)
        if other_keys is not None:
            for k,v in other_keys.items():
                fib.update_all_headers(k,v)
        fib.writeto(filename, clobber=clobber)
        
    @property
    def limits(self):
        rv = dict(x = None, y = None, z = None)
        xax = self.ds.coordinates.x_axis[self.axis]
        yax = self.ds.coordinates.y_axis[self.axis]
        xn = self.ds.coordinates.axis_name[xax]
        yn = self.ds.coordinates.axis_name[yax]
        rv[xn] = (self.bounds[0], self.bounds[1])
        rv[yn] = (self.bounds[2], self.bounds[3])
        return rv

class CylindricalFixedResolutionBuffer(FixedResolutionBuffer):
    """
    This object is a subclass of
    :class:`yt.visualization.fixed_resolution.FixedResolutionBuffer`
    that supports non-aligned input data objects, primarily cutting planes.
    """
    def __init__(self, data_source, radius, buff_size, antialias = True) :

        self.data_source = data_source
        self.ds = data_source.ds
        self.radius = radius
        self.buff_size = buff_size
        self.antialias = antialias
        self.data = {}
        
        ds = getattr(data_source, "ds", None)
        if ds is not None:
            ds.plots.append(weakref.proxy(self))

    def __getitem__(self, item) :
        if item in self.data: return self.data[item]
        buff = pixelize_cylinder(self.data_source["r"], self.data_source["dr"],
                                 self.data_source["theta"], self.data_source["dtheta"],
                                 self.buff_size, self.data_source[item].astype("float64"),
                                 self.radius)
        self[item] = buff
        return buff
        
class ObliqueFixedResolutionBuffer(FixedResolutionBuffer):
    """
    This object is a subclass of
    :class:`yt.visualization.fixed_resolution.FixedResolutionBuffer`
    that supports non-aligned input data objects, primarily cutting planes.
    """
    def __getitem__(self, item):
        if item in self.data: return self.data[item]
        indices = np.argsort(self.data_source['dx'])[::-1]
        bounds = []
        for b in self.bounds:
            if hasattr(b, "in_units"):
                b = float(b.in_units("code_length"))
            bounds.append(b)
        buff = _MPL.CPixelize( self.data_source['x'],   self.data_source['y'],   self.data_source['z'],
                               self.data_source['px'],  self.data_source['py'],
                               self.data_source['pdx'], self.data_source['pdy'], self.data_source['pdz'],
                               self.data_source.center, self.data_source._inv_mat, indices,
                               self.data_source[item],
                               self.buff_size[0], self.buff_size[1],
                               bounds).transpose()
        ia = ImageArray(buff, input_units=self.data_source[item].units,
                        info=self._get_info(item))
        self[item] = ia
        return ia 


class OffAxisProjectionFixedResolutionBuffer(FixedResolutionBuffer):
    """
    This object is a subclass of
    :class:`yt.visualization.fixed_resolution.FixedResolutionBuffer`
    that supports off axis projections.  This calls the volume renderer.
    """
    def __init__(self, data_source, bounds, buff_size, antialias = True,
                 periodic = False):
        self.data = {}
        FixedResolutionBuffer.__init__(self, data_source, bounds, buff_size, antialias, periodic)

    def __getitem__(self, item):
        if item in self.data: return self.data[item]
        mylog.info("Making a fixed resolutuion buffer of (%s) %d by %d" % \
            (item, self.buff_size[0], self.buff_size[1]))
        dd = self.data_source
        width = self.ds.arr((self.bounds[1] - self.bounds[0],
                             self.bounds[3] - self.bounds[2],
                             self.bounds[5] - self.bounds[4]))
        buff = off_axis_projection(dd.ds, dd.center, dd.normal_vector,
                                   width, dd.resolution, item,
                                   weight=dd.weight_field, volume=dd.volume,
                                   no_ghost=dd.no_ghost, interpolated=dd.interpolated,
                                   north_vector=dd.north_vector, method=dd.method)
        units = Unit(dd.ds.field_info[item].units, registry=dd.ds.unit_registry)
        if dd.weight_field is None and dd.method == "integrate":
            units *= Unit('cm', registry=dd.ds.unit_registry)
        ia = ImageArray(buff.swapaxes(0,1), input_units=units, info=self._get_info(item))
        self[item] = ia
        return ia 


