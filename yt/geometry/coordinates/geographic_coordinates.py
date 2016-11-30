"""
Definitions for geographic coordinate systems




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from .coordinate_handler import \
    CoordinateHandler, \
    _unknown_coord, \
    _get_coord_fields
from yt.utilities.lib.pixelization_routines import \
    pixelize_cylinder, pixelize_aitoff

class GeographicCoordinateHandler(CoordinateHandler):
    radial_axis = "altitude"
    name = "geographic"

    def __init__(self, ds, ordering = None):
        if not ordering:
            ordering = ('latitude', 'longitude', self.radial_axis)
        super(GeographicCoordinateHandler, self).__init__(ds, ordering)
        self.image_units = {}
        self.image_units[self.axis_id['latitude']] = (None, None)
        self.image_units[self.axis_id['longitude']] = (None, None)
        self.image_units[self.axis_id[self.radial_axis]] = ('deg', 'deg')

    def setup_fields(self, registry):
        # return the fields for r, z, theta
        registry.add_field(("index", "dx"), function=_unknown_coord)
        registry.add_field(("index", "dy"), function=_unknown_coord)
        registry.add_field(("index", "dz"), function=_unknown_coord)
        registry.add_field(("index", "x"), function=_unknown_coord)
        registry.add_field(("index", "y"), function=_unknown_coord)
        registry.add_field(("index", "z"), function=_unknown_coord)
        f1, f2 = _get_coord_fields(self.axis_id['latitude'], "")
        registry.add_field(("index", "dlatitude"), function = f1,
                           display_field = False,
                           units = "")
        registry.add_field(("index", "latitude"), function = f2,
                           display_field = False,
                           units = "")

        f1, f2 = _get_coord_fields(self.axis_id['longitude'], "")
        registry.add_field(("index", "dlongitude"), function = f1,
                           display_field = False,
                           units = "")
        registry.add_field(("index", "longitude"), function = f2,
                           display_field = False,
                           units = "")

        f1, f2 = _get_coord_fields(self.axis_id[self.radial_axis])
        registry.add_field(("index", "d%s" % (self.radial_axis,)), function = f1,
                           display_field = False,
                           units = "code_length")
        registry.add_field(("index", self.radial_axis), function = f2,
                           display_field = False,
                           units = "code_length")

        def _SphericalVolume(field, data):
            # r**2 sin theta dr dtheta dphi
            # We can use the transformed coordinates here.
            vol = data["index", "r"]**2.0
            vol *= data["index", "dr"]
            vol *= np.sin(data["index", "theta"])
            vol *= data["index", "dtheta"]
            vol *= data["index", "dphi"]
            return vol
        registry.add_field(("index", "cell_volume"),
                 function=_SphericalVolume,
                 units = "code_length**3")

        def _path_radial_axis(field, data):
            return data["index", "d%s" % self.radial_axis]
        registry.add_field(("index", "path_element_%s" % self.radial_axis),
                 function = _path_radial_axis,
                 units = "code_length")
        def _path_latitude(field, data):
            # We use r here explicitly
            return data["index", "r"] * \
                data["index", "dlatitude"] * np.pi/180.0
        registry.add_field(("index", "path_element_latitude"),
                 function = _path_latitude,
                 units = "code_length")
        def _path_longitude(field, data):
            # We use r here explicitly
            return data["index", "r"] \
                    * data["index", "dlongitude"] * np.pi/180.0 \
                    * np.sin((data["index", "latitude"] + 90.0) * np.pi/180.0)
        registry.add_field(("index", "path_element_longitude"),
                 function = _path_longitude,
                 units = "code_length")

        def _latitude_to_theta(field, data):
            # latitude runs from -90 to 90
            return (data["latitude"] + 90) * np.pi/180.0
        registry.add_field(("index", "theta"),
                 function = _latitude_to_theta,
                 units = "")
        def _dlatitude_to_dtheta(field, data):
            return data["dlatitude"] * np.pi/180.0
        registry.add_field(("index", "dtheta"),
                 function = _dlatitude_to_dtheta,
                 units = "")

        def _longitude_to_phi(field, data):
            # longitude runs from -180 to 180
            return (data["longitude"] + 180) * np.pi/180.0
        registry.add_field(("index", "phi"),
                 function = _longitude_to_phi,
                 units = "")
        def _dlongitude_to_dphi(field, data):
            return data["dlongitude"] * np.pi/180.0
        registry.add_field(("index", "dphi"),
                 function = _dlongitude_to_dphi,
                 units = "")

        self._setup_radial_fields(registry)

    def _setup_radial_fields(self, registry):
        # This stays here because we don't want to risk the field detector not
        # properly getting the data_source, etc.
        def _altitude_to_radius(field, data):
            surface_height = data.get_field_parameter("surface_height")
            if surface_height is None:
                if hasattr(data.ds, "surface_height"):
                    surface_height = data.ds.surface_height
                else:
                    surface_height = data.ds.quan(0.0, "code_length")
            return data["altitude"] + surface_height
        registry.add_field(("index", "r"),
                 function=_altitude_to_radius,
                 units = "code_length")
        registry.alias(("index", "dr"), ("index", "daltitude"))

    def _retrieve_radial_offset(self, data_source = None):
        # This returns the factor by which the radial field is multiplied and
        # the scalar its offset by.  Typically the "factor" will *only* be
        # either 1.0 or -1.0.  The order will be factor * r + offset.
        # Altitude is the radius from the central zone minus the radius of the
        # surface.  Depth to radius is negative value of depth plus the
        # outermost radius.
        surface_height = None
        if data_source is not None:
            surface_height = data_source.get_field_parameter("surface_height")
        if surface_height is None:
            if hasattr(self.ds, "surface_height"):
                surface_height = self.ds.surface_height
            else:
                surface_height = self.ds.quan(0.0, "code_length")
        return surface_height, 1.0

    def pixelize(self, dimension, data_source, field, bounds, size,
                 antialias = True, periodic = True):
        if self.axis_name[dimension] in ('latitude', 'longitude'):
            return self._cyl_pixelize(data_source, field, bounds, size,
                                          antialias, dimension)
        elif self.axis_name[dimension] == self.radial_axis:
            return self._ortho_pixelize(data_source, field, bounds, size,
                                        antialias, dimension, periodic)
        else:
            raise NotImplementedError

    def _ortho_pixelize(self, data_source, field, bounds, size, antialias,
                        dim, periodic):
        # For axis=2, x axis will be latitude, y axis will be longitude
        px = (data_source["px"].d + 90) * np.pi/180
        pdx = data_source["pdx"].d * np.pi/180
        py = (data_source["py"].d + 180) * np.pi/180
        pdy = data_source["pdy"].d * np.pi/180
        # First one in needs to be the equivalent of "theta", which is
        # longitude
        buff = pixelize_aitoff(py, pdy, px, pdx,
                               size, data_source[field], None,
                               None).transpose()
        return buff

    def _cyl_pixelize(self, data_source, field, bounds, size, antialias,
                      dimension):
        offset, factor = self._retrieve_radial_offset(data_source)
        r = factor * data_source['py'] + offset
        # Because of the axis-ordering, dimensions 0 and 1 both have r as py
        # and the angular coordinate as px.  But we need to figure out how to
        # convert our coordinate back to an actual angle, based on which
        # dimension we're in.
        pdx = data_source['pdx'].d * np.pi/180
        if self.axis_name[self.x_axis[dimension]] == 'latitude':
            px = (data_source['px'].d + 90) * np.pi/180
            do_transpose = True
        elif self.axis_name[self.x_axis[dimension]] == 'longitude':
            px = (data_source['px'].d + 180) * np.pi/180
            do_transpose = False
        else:
            # We should never get here!
            raise NotImplementedError
        buff = pixelize_cylinder(r, data_source['pdy'],
                                 px, pdx,
                                 size, data_source[field], bounds)
        if do_transpose:
            buff = buff.transpose()
        return buff

    def convert_from_cartesian(self, coord):
        raise NotImplementedError

    def convert_to_cartesian(self, coord):
        offset, factor = self._retrieve_radial_offset()
        if isinstance(coord, np.ndarray) and len(coord.shape) > 1:
            rad = self.axis_id[self.radial_axis]
            lon = self.axis_id['longitude']
            lat = self.axis_id['latitude']
            r = factor * coord[:,rad] + offset
            theta = coord[:,lon] * np.pi/180
            phi = coord[:,lat] * np.pi/180
            nc = np.zeros_like(coord)
            # r, theta, phi
            nc[:,lat] = np.cos(phi) * np.sin(theta)*r
            nc[:,lon] = np.sin(phi) * np.sin(theta)*r
            nc[:,rad] = np.cos(theta) * r
        else:
            a, b, c = coord
            theta = b * np.pi/180
            phi = a * np.pi/180
            r = factor * c + offset
            nc = (np.cos(phi) * np.sin(theta)*r,
                  np.sin(phi) * np.sin(theta)*r,
                  np.cos(theta) * r)
        return nc

    def convert_to_cylindrical(self, coord):
        raise NotImplementedError

    def convert_from_cylindrical(self, coord):
        raise NotImplementedError

    def convert_to_spherical(self, coord):
        raise NotImplementedError

    def convert_from_spherical(self, coord):
        raise NotImplementedError

    _image_axis_name = None
    @property
    def image_axis_name(self):    
        if self._image_axis_name is not None:
            return self._image_axis_name
        # This is the x and y axes labels that get displayed.  For
        # non-Cartesian coordinates, we usually want to override these for
        # Cartesian coordinates, since we transform them.
        rv = {self.axis_id['latitude']:
                 ('x / \\sin(\mathrm{latitude})',
                  'y / \\sin(\mathrm{latitude})'),
              self.axis_id['longitude']:
                 ('R', 'z'),
              self.axis_id[self.radial_axis]:
                 ('longitude', 'latitude')}
        for i in list(rv.keys()):
            rv[self.axis_name[i]] = rv[i]
            rv[self.axis_name[i].capitalize()] = rv[i]
        self._image_axis_name = rv
        return rv

    _x_pairs = (('latitude', 'longitude'),
                ('longitude', 'latitude'),
                ('altitude', 'latitude'))

    _y_pairs = (('latitude', 'altitude'),
                ('longitude', 'altitude'),
                ('altitude', 'longitude'))

    @property
    def period(self):
        return self.ds.domain_width

    def sanitize_center(self, center, axis):
        center, display_center = super(
            GeographicCoordinateHandler, self).sanitize_center(center, axis)
        name = self.axis_name[axis]
        if name == self.radial_axis:
            display_center = center
        elif name == 'latitude':
            display_center = (0.0 * display_center[0],
                              0.0 * display_center[1],
                              0.0 * display_center[2])
        elif name == 'longitude':
            ri = self.axis_id[self.radial_axis]
            c = (self.ds.domain_right_edge[ri] +
                 self.ds.domain_left_edge[ri])/2.0
            display_center = [0.0 * display_center[0], 
                              0.0 * display_center[1],
                              0.0 * display_center[2]]
            display_center[self.axis_id['latitude']] = c
        return center, display_center

    def sanitize_width(self, axis, width, depth):
        name = self.axis_name[axis]
        if width is not None:
            width = super(GeographicCoordinateHandler, self).sanitize_width(
              axis, width, depth)
        elif name == self.radial_axis:
            rax = self.radial_axis
            width = [self.ds.domain_width[self.x_axis[rax]],
                     self.ds.domain_width[self.y_axis[rax]]]
        elif name == 'latitude':
            ri = self.axis_id[self.radial_axis]
            # Remember, in spherical coordinates when we cut in theta,
            # we create a conic section
            width = [2.0*self.ds.domain_width[ri],
                     2.0*self.ds.domain_width[ri]]
        elif name == 'longitude':
            ri = self.axis_id[self.radial_axis]
            width = [self.ds.domain_width[ri],
                     2.0*self.ds.domain_width[ri]]
        return width

class InternalGeographicCoordinateHandler(GeographicCoordinateHandler):
    radial_axis = "depth"
    name = "internal_geographic"
    
    def _setup_radial_fields(self, registry):
        # Altitude is the radius from the central zone minus the radius of the
        # surface.
        def _depth_to_radius(field, data):
            outer_radius = data.get_field_parameter("outer_radius")
            if outer_radius is None:
                if hasattr(data.ds, "outer_radius"):
                    outer_radius = data.ds.outer_radius
                else:
                    # Otherwise, we assume that the depth goes to full depth,
                    # so we can look at the domain right edge in depth.
                    rax = self.axis_id[self.radial_axis]
                    outer_radius = data.ds.domain_right_edge[rax]
            return -1.0 * data["depth"] + outer_radius
        registry.add_field(("index", "r"),
                 function=_depth_to_radius,
                 units = "code_length")
        registry.alias(("index", "dr"), ("index", "ddepth"))

    def _retrieve_radial_offset(self, data_source = None):
        # Depth means switching sign and adding to full radius
        outer_radius = None
        if data_source is not None:
            outer_radius = data_source.get_field_parameter("outer_radius")
        if outer_radius is None:
            if hasattr(self.ds, "outer_radius"):
                outer_radius = self.ds.outer_radius
            else:
                # Otherwise, we assume that the depth goes to full depth,
                # so we can look at the domain right edge in depth.
                rax = self.axis_id[self.radial_axis]
                outer_radius = self.ds.domain_right_edge[rax]
        return outer_radius, -1.0

    _x_pairs = (('latitude', 'longitude'),
                ('longitude', 'latitude'),
                ('depth', 'latitude'))

    _y_pairs = (('latitude', 'depth'),
                ('longitude', 'depth'),
                ('depth', 'longitude'))

    def sanitize_center(self, center, axis):
        center, display_center = super(
            GeographicCoordinateHandler, self).sanitize_center(center, axis)
        name = self.axis_name[axis]
        if name == self.radial_axis:
            display_center = center
        elif name == 'latitude':
            display_center = (0.0 * display_center[0],
                              0.0 * display_center[1],
                              0.0 * display_center[2])
        elif name == 'longitude':
            ri = self.axis_id[self.radial_axis]
            offset, factor = self._retrieve_radial_offset()
            outermost = factor * self.ds.domain_left_edge[ri] + offset
            display_center = [0.0 * display_center[0], 
                              0.0 * display_center[1],
                              0.0 * display_center[2]]
            display_center[self.axis_id['latitude']] = outermost/2.0
        return center, display_center

    def sanitize_width(self, axis, width, depth):
        name = self.axis_name[axis]
        if width is not None:
            width = super(GeographicCoordinateHandler, self).sanitize_width(
              axis, width, depth)
        elif name == self.radial_axis:
            rax = self.radial_axis
            width = [self.ds.domain_width[self.y_axis[rax]],
                     self.ds.domain_width[self.x_axis[rax]]]
        elif name == 'latitude':
            ri = self.axis_id[self.radial_axis]
            # Remember, in spherical coordinates when we cut in theta,
            # we create a conic section
            offset, factor = self._retrieve_radial_offset()
            outermost = factor * self.ds.domain_left_edge[ri] + offset
            width = [2.0*outermost, 2.0*outermost]
        elif name == 'longitude':
            ri = self.axis_id[self.radial_axis]
            offset, factor = self._retrieve_radial_offset()
            outermost = factor * self.ds.domain_left_edge[ri] + offset
            width = [outermost, 2.0*outermost]
        return width

