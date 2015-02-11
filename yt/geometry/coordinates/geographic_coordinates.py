"""
Geographic fields




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
import yt.visualization._MPL as _MPL
from yt.utilities.lib.misc_utilities import \
    pixelize_cylinder, pixelize_aitoff

class GeographicCoordinateHandler(CoordinateHandler):

    def __init__(self, ds, ordering = 'latlonalt'):
        if ordering != 'latlonalt': raise NotImplementedError
        super(GeographicCoordinateHandler, self).__init__(ds)

    def setup_fields(self, registry):
        # return the fields for r, z, theta
        registry.add_field(("index", "dx"), function=_unknown_coord)
        registry.add_field(("index", "dy"), function=_unknown_coord)
        registry.add_field(("index", "dz"), function=_unknown_coord)
        registry.add_field(("index", "x"), function=_unknown_coord)
        registry.add_field(("index", "y"), function=_unknown_coord)
        registry.add_field(("index", "z"), function=_unknown_coord)
        f1, f2 = _get_coord_fields(0, "")
        registry.add_field(("index", "dlatitude"), function = f1,
                           display_field = False,
                           units = "")
        registry.add_field(("index", "latitude"), function = f2,
                           display_field = False,
                           units = "")

        f1, f2 = _get_coord_fields(1, "")
        registry.add_field(("index", "dlongitude"), function = f1,
                           display_field = False,
                           units = "")
        registry.add_field(("index", "longitude"), function = f2,
                           display_field = False,
                           units = "")

        f1, f2 = _get_coord_fields(2)
        registry.add_field(("index", "daltitude"), function = f1,
                           display_field = False,
                           units = "code_length")
        registry.add_field(("index", "altitude"), function = f2,
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

        def _path_altitude(field, data):
            return data["index", "daltitude"]
        registry.add_field(("index", "path_element_altitude"),
                 function = _path_altitude,
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

        # Altitude is the radius from the central zone minus the radius of the
        # surface.
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

    def pixelize(self, dimension, data_source, field, bounds, size,
                 antialias = True, periodic = True):
        if dimension in (0, 1):
            return self._cyl_pixelize(data_source, field, bounds, size,
                                          antialias, dimension)
        elif dimension == 2:
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
        surface_height = data_source.get_field_parameter("surface_height")
        if surface_height is None:
            if hasattr(data_source.ds, "surface_height"):
                surface_height = data_source.ds.surface_height
            else:
                surface_height = data_source.ds.quan(0.0, "code_length")
        r = data_source['py'] + surface_height
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
        raise NotImplementedError

    def convert_to_cylindrical(self, coord):
        raise NotImplementedError

    def convert_from_cylindrical(self, coord):
        raise NotImplementedError

    def convert_to_spherical(self, coord):
        raise NotImplementedError

    def convert_from_spherical(self, coord):
        raise NotImplementedError

    # Despite being mutables, we uses these here to be clear about how these
    # are generated and to ensure that they are not re-generated unnecessarily
    axis_name = { 0  : 'latitude',  1  : 'longitude',  2  : 'altitude',
                 'latitude' : 'latitude',
                 'longitude' : 'longitude', 
                 'altitude' : 'altitude',
                 'Latitude' : 'latitude',
                 'Longitude' : 'longitude', 
                 'Altitude' : 'altitude',
                 'lat' : 'latitude',
                 'lon' : 'longitude', 
                 'alt' : 'altitude' }

    _image_axis_name = None
    @property
    def image_axis_name(self):    
        if self._image_axis_name is not None:
            return self._image_axis_name
        # This is the x and y axes labels that get displayed.  For
        # non-Cartesian coordinates, we usually want to override these for
        # Cartesian coordinates, since we transform them.
        rv = {0: ('x / \\sin(\mathrm{longitude})',
                  'y / \\sin(\mathrm{longitude})'),
              1: ('R', 'z'),
              2: ('longitude', 'latitude')}
        for i in rv.keys():
            rv[self.axis_name[i]] = rv[i]
            rv[self.axis_name[i].upper()] = rv[i]
        self._image_axis_name = rv
        return rv

    axis_id = { 'latitude' : 0, 'longitude' : 1, 'altitude' : 2,
                 0  : 0,  1  : 1,  2  : 2}

    x_axis = { 'latitude' : 1, 'longitude' : 0, 'altitude' : 0,
                0  : 1,  1  : 0,  2  : 0}

    y_axis = { 'latitude' : 2, 'longitude' : 2, 'altitude' : 1,
                0  : 2,  1  : 2,  2  : 1}

    @property
    def period(self):
        return self.ds.domain_width

    def sanitize_center(self, center, axis):
        center, display_center = super(
            GeographicCoordinateHandler, self).sanitize_center(center, axis)
        if axis == 2:
            display_center = center
        elif axis == 1:
            display_center = (0.0 * display_center[0],
                              0.0 * display_center[1],
                              0.0 * display_center[2])
        elif axis == 0:
            display_center = (self.ds.domain_width[0]/2.0 +
                              self.ds.domain_left_edge[0],
                              0.0 * display_center[1],
                              0.0 * display_center[2])
        return center, display_center

    def convert_to_cartesian(self, coord):
        if isinstance(coord, np.ndarray) and len(coord.shape) > 1:
            r = coord[:,2] + self.ds.surface_height
            theta = coord[:,1] * np.pi/180
            phi = coord[:,0] * np.pi/180
            nc = np.zeros_like(coord)
            # r, theta, phi
            nc[:,0] = np.cos(phi) * np.sin(theta)*r
            nc[:,1] = np.sin(phi) * np.sin(theta)*r
            nc[:,2] = np.cos(theta) * r
        else:
            a, b, c = coord
            theta = b * np.pi/180
            phi = a * np.pi/180
            r = self.ds.surface_height + c
            nc = (np.cos(phi) * np.sin(theta)*r,
                  np.sin(phi) * np.sin(theta)*r,
                  np.cos(theta) * r)
        return nc

