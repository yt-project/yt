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

from yt.frontends.ytdata.utilities import \
    save_as_dataset
from yt.funcs import \
    get_output_filename, \
    mylog, \
    ensure_list, \
    deprecate
from .volume_rendering.api import off_axis_projection
from .fixed_resolution_filters import apply_filter, filter_registry
from yt.data_objects.image_array import ImageArray
from yt.utilities.lib.pixelization_routines import \
    pixelize_cylinder
from yt.utilities.lib.api import add_points_to_greyscale_image
from yt.frontends.stream.api import load_uniform_grid

import numpy as np
import weakref
import types

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
    multiple fields.  It acts as a standard YTDataContainer object, such that
    dict-style access returns an image of a given field.

    Parameters
    ----------
    data_source : :class:`yt.data_objects.construction_data_containers.YTQuadTreeProj` or :class:`yt.data_objects.selection_data_containers.YTSlice`
        This is the source to be pixelized, which can be a projection, slice or
        cutting plane.
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

    Examples
    --------
    To make a projection and then several images, you can generate a
    single FRB and then access multiple fields:

    >>> proj = ds.proj(0, "density")
    >>> frb1 = FixedResolutionBuffer(proj, (0.2, 0.3, 0.4, 0.5),
    ...                              (1024, 1024))
    >>> print frb1["density"].max()
    1.0914e-9 g/cm**3
    >>> print frb1["temperature"].max()
    104923.1 K
    """
    _exclude_fields = ('pz','pdz','dx','x','y','z',
        'r', 'dr', 'phi', 'dphi', 'theta', 'dtheta',
                       ('index','dx'),('index','x'),('index','y'),('index','z'),
                       ('index', 'r'), ('index', 'dr'),
                       ('index', 'phi'), ('index', 'dphi'),
                       ('index', 'theta'), ('index', 'dtheta'))
    def __init__(self, data_source, bounds, buff_size, antialias = True,
                 periodic = False):
        self.data_source = data_source
        self.ds = data_source.ds
        self.bounds = bounds
        self.buff_size = (int(buff_size[0]), int(buff_size[1]))
        self.antialias = antialias
        self.data = {}
        self._filters = []
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

        self.setup_filters()

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

        for name, (args, kwargs) in self._filters:
            buff = filter_registry[name](*args[1:], **kwargs).apply(buff)

        # FIXME FIXME FIXME we shouldn't need to do this for projections
        # but that will require fixing data object access for particle
        # projections
        try:
            if hasattr(item, 'name'):
                it = item.name
            else:
                it = item
            units = self.data_source._projected_units[it]
        except (KeyError, AttributeError):
            units = self.data_source[item].units

        ia = ImageArray(buff, input_units=units, info=self._get_info(item))
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

        info['label'] = finfo.get_latex_display_name()

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
        distance in the x-coordinate.

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
        distance in the y-coordinate.

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

    def set_unit(self, field, unit, equivalency=None, equivalency_kwargs=None):
        """Sets a new unit for the requested field

        parameters
        ----------
        field : string or field tuple
           The name of the field that is to be changed.

        unit : string or Unit object
           The name of the new unit.

        equivalency : string, optional
           If set, the equivalency to use to convert the current units to
           the new requested unit. If None, the unit conversion will be done
           without an equivalency

        equivalency_kwargs : string, optional
           Keyword arguments to be passed to the equivalency. Only used if
           ``equivalency`` is set.
        """
        if equivalency_kwargs is None:
            equivalency_kwargs = {}
        field = self.data_source._determine_fields(field)[0]
        if equivalency is None:
            self[field].convert_to_units(unit)
        else:
            equiv_array = self[field].to_equivalent(
                unit, equivalency, **equivalency_kwargs)
            # equiv_array isn't necessarily an ImageArray. This is an issue
            # inherent to the way the unit system handles YTArray
            # subclasses and I don't see how to modify the unit system to
            # fix this. Instead, we paper over this issue and hard code
            # that equiv_array is an ImageArray
            self[field] = ImageArray(
                equiv_array, equiv_array.units, equiv_array.units.registry,
                self[field].info)


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
        if fields is None: fields = list(self.data.keys())
        output = h5py.File(filename, "a")
        for field in fields:
            output.create_dataset(field,data=self[field])
        output.close()

    def export_fits(self, filename, fields=None, overwrite=False,
                    other_keys=None, units="cm", **kwargs):
        r"""Export a set of pixelized fields to a FITS file.

        This will export a set of FITS images of either the fields specified
        or all the fields already in the object.

        Parameters
        ----------
        filename : string
            The name of the FITS file to be written.
        fields : list of strings
            These fields will be pixelized and output. If "None", the keys of the
            FRB will be used.
        overwrite : boolean
            If the file exists, this governs whether we will overwrite.
        other_keys : dictionary, optional
            A set of header keys and values to write into the FITS header.
        units : string, optional
            the length units that the coordinates are written in, default 'cm'.
        """

        from yt.visualization.fits_image import FITSImageData

        if fields is None:
            fields = list(self.data.keys())
        else:
            fields = ensure_list(fields)

        if len(fields) == 0:
            raise RuntimeError(
                "No fields to export. Either pass a field or list of fields to "
                "export_fits or access a field from the fixed resolution buffer "
                "object."
            )

        fib = FITSImageData(self, fields=fields, units=units)
        if other_keys is not None:
            for k,v in other_keys.items():
                fib.update_all_headers(k,v)
        fib.writeto(filename, overwrite=overwrite, **kwargs)

    def export_dataset(self, fields=None, nprocs=1):
        r"""Export a set of pixelized fields to an in-memory dataset that can be
        analyzed as any other in yt. Unit information and other parameters (e.g.,
        geometry, current_time, etc.) will be taken from the parent dataset.

        Parameters
        ----------
        fields : list of strings, optional
            These fields will be pixelized and output. If "None", the keys of the
            FRB will be used.
        nprocs: integer, optional
            If greater than 1, will create this number of subarrays out of data

        Examples
        --------
        >>> import yt
        >>> ds = yt.load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150")
        >>> slc = ds.slice(2, 0.0)
        >>> frb = slc.to_frb((500.,"kpc"), 500)
        >>> ds2 = frb.export_dataset(fields=["density","temperature"], nprocs=32)
        """
        nx, ny = self.buff_size
        data = {}
        if fields is None:
            fields = list(self.keys())
        for field in fields:
            arr = self[field]
            data[field] = (arr.d.T.reshape(nx,ny,1), str(arr.units))
        bounds = [b.in_units("code_length").v for b in self.bounds]
        bbox = np.array([[bounds[0],bounds[1]],[bounds[2],bounds[3]],[0.,1.]])
        return load_uniform_grid(data, [nx,ny,1],
                                 length_unit=self.ds.length_unit,
                                 bbox=bbox,
                                 sim_time=self.ds.current_time.in_units("s").v,
                                 mass_unit=self.ds.mass_unit,
                                 time_unit=self.ds.time_unit,
                                 velocity_unit=self.ds.velocity_unit,
                                 magnetic_unit=self.ds.magnetic_unit,
                                 periodicity=(False,False,False),
                                 geometry=self.ds.geometry,
                                 nprocs=nprocs)

    def save_as_dataset(self, filename=None, fields=None):
        r"""Export a fixed resolution buffer to a reloadable yt dataset.

        This function will take a fixed resolution buffer and output a 
        dataset containing either the fields presently existing or fields 
        given in the ``fields`` list.  The resulting dataset can be
        reloaded as a yt dataset.

        Parameters
        ----------
        filename : str, optional
            The name of the file to be written.  If None, the name 
            will be a combination of the original dataset and the type 
            of data container.
        fields : list of strings or tuples, optional
            If this is supplied, it is the list of fields to be saved to
            disk.  If not supplied, all the fields that have been queried
            will be saved.

        Returns
        -------
        filename : str
            The name of the file that has been created.

        Examples
        --------

        >>> import yt
        >>> ds = yt.load("enzo_tiny_cosmology/DD0046/DD0046")
        >>> proj = ds.proj("density", "x", weight_field="density")
        >>> frb = proj.to_frb(1.0, (800, 800))
        >>> fn = frb.save_as_dataset(fields=["density"])
        >>> ds2 = yt.load(fn)
        >>> print (ds2.data["density"])
        [[  1.25025353e-30   1.25025353e-30   1.25025353e-30 ...,   7.90820691e-31
            7.90820691e-31   7.90820691e-31]
         [  1.25025353e-30   1.25025353e-30   1.25025353e-30 ...,   7.90820691e-31
            7.90820691e-31   7.90820691e-31]
         [  1.25025353e-30   1.25025353e-30   1.25025353e-30 ...,   7.90820691e-31
            7.90820691e-31   7.90820691e-31]
         ...,
         [  1.55834239e-30   1.55834239e-30   1.55834239e-30 ...,   8.51353199e-31
            8.51353199e-31   8.51353199e-31]
         [  1.55834239e-30   1.55834239e-30   1.55834239e-30 ...,   8.51353199e-31
            8.51353199e-31   8.51353199e-31]
         [  1.55834239e-30   1.55834239e-30   1.55834239e-30 ...,   8.51353199e-31
            8.51353199e-31   8.51353199e-31]] g/cm**3

        """

        keyword = "%s_%s_frb" % (str(self.ds), self.data_source._type_name)
        filename = get_output_filename(filename, keyword, ".h5")

        data = {}
        if fields is not None:
            for f in self.data_source._determine_fields(fields):
                data[f] = self[f]
        else:
            data.update(self.data)

        ftypes = dict([(field, "grid") for field in data])
        extra_attrs = dict([(arg, getattr(self.data_source, arg, None))
                            for arg in self.data_source._con_args +
                            self.data_source._tds_attrs])
        extra_attrs["con_args"] = self.data_source._con_args
        extra_attrs["left_edge"] = self.ds.arr([self.bounds[0],
                                                self.bounds[2]])
        extra_attrs["right_edge"] = self.ds.arr([self.bounds[1],
                                                 self.bounds[3]])
        extra_attrs["ActiveDimensions"] = self.buff_size
        extra_attrs["level"] = 0
        extra_attrs["data_type"] = "yt_frb"
        extra_attrs["container_type"] = self.data_source._type_name
        extra_attrs["dimensionality"] = self.data_source._dimensionality
        save_as_dataset(self.ds, filename, data, field_types=ftypes,
                        extra_attrs=extra_attrs)

        return filename

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

    def setup_filters(self):
        ignored = ['FixedResolutionBufferFilter']
        for key in filter_registry:
            if key in ignored:
                continue
            filtername = filter_registry[key]._filter_name
            FilterMaker = filter_registry[key]
            filt = apply_filter(FilterMaker)
            filt.__doc__ = FilterMaker.__doc__
            self.__dict__['apply_' + filtername] = \
                types.MethodType(filt, self)

class ObliqueFixedResolutionBuffer(FixedResolutionBuffer):
    @deprecate("FixedResolutionBuffer")
    def __init__(self, *args, **kwargs):
        super(ObliqueFixedResolutionBuffer, self).__init__(*args, **kwargs)

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
        buff = np.zeros(self.buff_size, dtype="f8")
        pixelize_cylinder(buff, self.data_source["r"], self.data_source["dr"],
                          self.data_source["theta"], self.data_source["dtheta"],
                          self.data_source[item].astype("float64"), self.radius)
        self[item] = buff
        return buff

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
        mylog.info("Making a fixed resolution buffer of (%s) %d by %d" % \
            (item, self.buff_size[0], self.buff_size[1]))
        dd = self.data_source
        width = self.ds.arr((self.bounds[1] - self.bounds[0],
                             self.bounds[3] - self.bounds[2],
                             self.bounds[5] - self.bounds[4]))
        buff = off_axis_projection(dd.dd, dd.center, dd.normal_vector,
                                   width, dd.resolution, item,
                                   weight=dd.weight_field, volume=dd.volume,
                                   no_ghost=dd.no_ghost,
                                   interpolated=dd.interpolated,
                                   north_vector=dd.north_vector,
                                   method=dd.method)
        ia = ImageArray(buff.swapaxes(0,1), info=self._get_info(item))
        self[item] = ia
        return ia


class ParticleImageBuffer(FixedResolutionBuffer):
    """

    This object is a subclass of
    :class:`yt.visualization.fixed_resolution.FixedResolutionBuffer`
    that supports particle plots. It splats points onto an image
    buffer.

    """
    def __init__(self, data_source, bounds, buff_size, antialias=True,
                 periodic=False):
        self.data = {}
        FixedResolutionBuffer.__init__(self, data_source, bounds, buff_size,
                                       antialias, periodic)

        # set up the axis field names
        axis = self.axis
        xax = self.ds.coordinates.x_axis[axis]
        yax = self.ds.coordinates.y_axis[axis]
        ax_field_template = 'particle_position_%s'
        self.x_field = ax_field_template % self.ds.coordinates.axis_name[xax]
        self.y_field = ax_field_template % self.ds.coordinates.axis_name[yax]

    def __getitem__(self, item):
        if item in self.data:
            return self.data[item]

        mylog.info("Splatting (%s) onto a %d by %d mesh" %
                (item, self.buff_size[0], self.buff_size[1]))

        bounds = []
        for b in self.bounds:
            if hasattr(b, "in_units"):
                b = float(b.in_units("code_length"))
            bounds.append(b)

        ftype = item[0]
        x_data = self.data_source.dd[ftype, self.x_field]
        y_data = self.data_source.dd[ftype, self.y_field]
        data = self.data_source.dd[item]

        # handle periodicity
        dx = x_data.in_units("code_length").d - bounds[0]
        dy = y_data.in_units("code_length").d - bounds[2]
        if self.periodic:
            dx %= float(self._period[0].in_units("code_length"))
            dy %= float(self._period[1].in_units("code_length"))

        # convert to pixels
        px = dx / (bounds[1] - bounds[0])
        py = dy / (bounds[3] - bounds[2])

        # select only the particles that will actually show up in the image
        mask = np.logical_and(np.logical_and(px >= 0.0, px <= 1.0),
                              np.logical_and(py >= 0.0, py <= 1.0))

        weight_field = self.data_source.weight_field
        if weight_field is None:
            weight_data = np.ones_like(data.v)
        else:
            weight_data = self.data_source.dd[weight_field]
        splat_vals = weight_data[mask]*data[mask]

        # splat particles
        buff = np.zeros(self.buff_size)
        buff_mask = np.zeros(self.buff_size).astype('int')
        add_points_to_greyscale_image(buff,
                                      buff_mask,
                                      px[mask],
                                      py[mask],
                                      splat_vals)
        # remove values in no-particle region
        buff[buff_mask==0] = np.nan
        ia = ImageArray(buff, input_units=data.units,
                        info=self._get_info(item))

        # divide by the weight_field, if needed
        if weight_field is not None:
            weight_buff = np.zeros(self.buff_size)
            weight_buff_mask = np.zeros(self.buff_size).astype('int')
            add_points_to_greyscale_image(weight_buff,
                                          weight_buff_mask,
                                          px[mask],
                                          py[mask],
                                          weight_data[mask])
            weight_array = ImageArray(weight_buff,
                                      input_units=weight_data.units,
                                      info=self._get_info(item))
            # remove values in no-particle region
            weight_buff[weight_buff_mask==0] = np.nan
            locs = np.where(weight_array > 0)
            ia[locs] /= weight_array[locs]

        self.data[item] = ia
        return self.data[item]

    # over-ride the base class version, since we don't want to exclude
    # particle fields
    def _get_data_source_fields(self):
        exclude = self.data_source._key_fields + list(self._exclude_fields)
        fields = getattr(self.data_source, "fields", [])
        fields += getattr(self.data_source, "field_data", {}).keys()
        for f in fields:
            if f not in exclude:
                self[f]
