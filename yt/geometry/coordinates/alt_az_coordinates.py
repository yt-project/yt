"""
Definitions for geographic coordinate systems




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2018, yt Development Team.
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

class AltAzCoordinateHandler(CoordinateHandler):
    name = "alt_az"

    def __init__(self, ds, ordering = None):
        if not ordering:
            ordering = ('range', 'azimuth', 'elevation')
        super(AltAzCoordinateHandler, self).__init__(ds, ordering)
        self.image_units = {}
        self.image_units[self.axis_id['azimuth']] = (None, None)
        self.image_units[self.axis_id['elevation']] = (None, None)
        self.image_units[self.axis_id["range"]] = ('deg', 'deg')

    def setup_fields(self, registry):
        registry.add_field(("index", "dx"), sampling_type="cell",  function=_unknown_coord)
        registry.add_field(("index", "dy"), sampling_type="cell",  function=_unknown_coord)
        registry.add_field(("index", "dz"), sampling_type="cell",  function=_unknown_coord)
        registry.add_field(("index", "x"), sampling_type="cell",  function=_unknown_coord)
        registry.add_field(("index", "y"), sampling_type="cell",  function=_unknown_coord)
        registry.add_field(("index", "z"), sampling_type="cell",  function=_unknown_coord)
        f1, f2 = _get_coord_fields(self.axis_id['azimuth'], "")
        registry.add_field(("index", "dazimuth"), sampling_type="cell",  function = f1,
                           display_field = False,
                           units = "")
        registry.add_field(("index", "azimuth"), sampling_type="cell",  function = f2,
                           display_field = False,
                           units = "")

        f1, f2 = _get_coord_fields(self.axis_id['elevation'], "")
        registry.add_field(("index", "delevation"), sampling_type="cell",  function = f1,
                           display_field = False,
                           units = "")
        registry.add_field(("index", "elevation"), sampling_type="cell",  function = f2,
                           display_field = False,
                           units = "")

        f1, f2 = _get_coord_fields(self.axis_id["range"])
        registry.add_field(("index", "drange"), sampling_type="cell",  function = f1,
                           display_field = False,
                           units = "code_length")
        registry.add_field(("index", "range"), sampling_type="cell",  function = f2,
                           display_field = False,
                           units = "code_length")

        def _SphericalVolume(field, data):
            # We can use the transformed coordinates here.
            # Here we compute the spherical volume element exactly
            r = data["index", "r"]
            dr = data["index", "dr"]
            theta = data["index", "theta"]
            dtheta = data["index", "dtheta"]
            vol = ((r+0.5*dr)**3-(r-0.5*dr)**3) / 3.0
            vol *= np.cos(theta-0.5*dtheta)-np.cos(theta+0.5*dtheta)
            vol *= data["index", "dphi"]
            return vol
        registry.add_field(("index", "cell_volume"), sampling_type="cell",
                 function=_SphericalVolume,
                 units = "code_length**3")

        def _path_range(field, data):
            return data["index", "drange"]
        registry.add_field(("index", "path_element_range"), sampling_type="cell",
                 function = _path_range,
                 units = "code_length")
        def _path_elevation(field, data):
            # We use r here explicitly
            return data["index", "r"] * \
                data["index", "delevation"] * np.pi/180.0
        registry.add_field(("index", "path_element_elevation"), sampling_type="cell",
                 function = _path_elevation,
                 units = "code_length")
        def _path_azimuth(field, data):
            # We use r here explicitly
            return data["index", "r"] \
                * data["index", "dazimuth"] * np.pi/180.0 \
                * np.sin((data["index", "elevation"] + 90.0) * np.pi/180.0)
        registry.add_field(("index", "path_element_azimuth"), sampling_type="cell",
                 function = _path_azimuth,
                 units = "code_length")

        def _elevation_to_theta(field, data):
            # elevation runs from -90 to 90
            return (data["elevation"] + 90) * np.pi/180.0
        registry.add_field(("index", "theta"), sampling_type="cell",
                 function = _elevation_to_theta,
                 units = "")
        def _delevation_to_dtheta(field, data):
            return data["delevation"] * np.pi/180.0
        registry.add_field(("index", "dtheta"), sampling_type="cell",
                 function = _delevation_to_dtheta,
                 units = "")

        def _azimuth_to_phi(field, data):
            # azimuth runs from -180 to 180
            return (data["azimuth"] + 180) * np.pi/180.0
        registry.add_field(("index", "phi"), sampling_type="cell",
                 function = _azimuth_to_phi,
                 units = "")
        def _dazimuth_to_dphi(field, data):
            return data["dazimuth"] * np.pi/180.0
        registry.add_field(("index", "dphi"), sampling_type="cell",
                 function = _dazimuth_to_dphi,
                 units = "")

        self._setup_radial_fields(registry)

    def _setup_radial_fields(self, registry):
        # This stays here because we don't want to risk the field detector not
        # properly getting the data_source, etc.
        def _range_to_radius(field, data):
            return data["range"]
        registry.add_field(("index", "r"), sampling_type="cell",
                 function=_range_to_radius,
                 units = "code_length")
        registry.alias(("index", "dr"), ("index", "drange"))

    def pixelize(self, dimension, data_source, field, bounds, size,
                 antialias = True, periodic = True):
        if self.axis_name[dimension] in ('elevation', 'azimuth'):
            return self._cyl_pixelize(data_source, field, bounds, size,
                                          antialias, dimension)
        elif self.axis_name[dimension] == "range":
            return self._ortho_pixelize(data_source, field, bounds, size,
                                        antialias, dimension, periodic)
        else:
            raise NotImplementedError

    def _ortho_pixelize(self, data_source, field, bounds, size, antialias,
                        dim, periodic):
        # For a radial axis, px will correspond to azimuth and py will
        # correspond to elevation.
        px = (data_source["px"].d + 180) * np.pi/180
        pdx = data_source["pdx"].d * np.pi/180
        py = (data_source["py"].d + 90) * np.pi/180
        pdy = data_source["pdy"].d * np.pi/180
        # First one in needs to be the equivalent of "theta", which is
        # azimuth
        buff = pixelize_aitoff(px, pdx, py, pdy,
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
        if self.axis_name[self.x_axis[dimension]] == 'elevation':
            px = (data_source['px'].d + 90) * np.pi/180
            do_transpose = True
        elif self.axis_name[self.x_axis[dimension]] == 'azimuth':
            px = (data_source['px'].d + 180) * np.pi/180
            do_transpose = False
        else:
            # We should never get here!
            raise NotImplementedError
        buff = np.zeros((size[1], size[0]), dtype="f8")
        pixelize_cylinder(buff, r, data_source['pdy'],
                          px, pdx, data_source[field], bounds)
        if do_transpose:
            buff = buff.transpose()
        return buff

    def convert_from_cartesian(self, coord):
        raise NotImplementedError

    def convert_to_cartesian(self, coord):
        offset, factor = self._retrieve_radial_offset()
        if isinstance(coord, np.ndarray) and len(coord.shape) > 1:
            rad = self.axis_id["range"]
            az = self.axis_id['azimuth']
            el = self.axis_id['elevation']
            r = factor * coord[:,rad] + offset
            theta = coord[:,az] * np.pi/180
            phi = coord[:,el] * np.pi/180
            nc = np.zeros_like(coord)
            # r, theta, phi
            nc[:,el] = np.cos(phi) * np.sin(theta)*r
            nc[:,az] = np.sin(phi) * np.sin(theta)*r
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
        rv = {self.axis_id['elevation']:
                 ('x / \\sin(\mathrm{elevation})',
                  'y / \\sin(\mathrm{elevation})'),
              self.axis_id['azimuth']:
                 ('R', 'z'),
              self.axis_id["range"]:
                 ('azimuth', 'elevation')}
        for i in list(rv.keys()):
            rv[self.axis_name[i]] = rv[i]
            rv[self.axis_name[i].capitalize()] = rv[i]
        self._image_axis_name = rv
        return rv

    _x_pairs = (('elevation', 'azimuth'),
                ('azimuth', 'elevation'),
                ('altitude', 'azimuth'))

    _y_pairs = (('elevation', 'altitude'),
                ('azimuth', 'altitude'),
                ('altitude', 'elevation'))

    @property
    def period(self):
        return self.ds.domain_width

    def sanitize_center(self, center, axis):
        center, display_center = super(
            AltAzCoordinateHandler, self).sanitize_center(center, axis)
        name = self.axis_name[axis]
        if name == "range":
            display_center = center
        elif name == 'elevation':
            display_center = (0.0 * display_center[0],
                              0.0 * display_center[1],
                              0.0 * display_center[2])
        elif name == 'azimuth':
            ri = self.axis_id["range"]
            c = (self.ds.domain_right_edge[ri] +
                 self.ds.domain_left_edge[ri])/2.0
            display_center = [0.0 * display_center[0],
                              0.0 * display_center[1],
                              0.0 * display_center[2]]
            display_center[self.axis_id['elevation']] = c
        return center, display_center

    def sanitize_width(self, axis, width, depth):
        name = self.axis_name[axis]
        if width is not None:
            width = super(AltAzCoordinateHandler, self).sanitize_width(
              axis, width, depth)
        elif name == "range":
            width = [self.ds.domain_width[self.x_axis[name]],
                     self.ds.domain_width[self.y_axis[name]]]
        elif name == 'elevation':
            ri = self.axis_id["range"]
            # Remember, in spherical coordinates when we cut in theta,
            # we create a conic section
            width = [2.0*self.ds.domain_width[ri],
                     2.0*self.ds.domain_width[ri]]
        elif name == 'azimuth':
            ri = self.axis_id["range"]
            width = [self.ds.domain_width[ri],
                     2.0*self.ds.domain_width[ri]]
        return width
