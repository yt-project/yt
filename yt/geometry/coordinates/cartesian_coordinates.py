"""
Definitions for cartesian coordinate systems




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
    _get_coord_fields, \
    _get_vert_fields, \
    cartesian_to_cylindrical, \
    cylindrical_to_cartesian
from yt.funcs import mylog
from yt.units.yt_array import uvstack, YTArray
from yt.utilities.lib.pixelization_routines import \
    pixelize_element_mesh, pixelize_off_axis_cartesian, \
    pixelize_cartesian, pixelize_cartesian_nodal, \
    pixelize_element_mesh_line
from yt.data_objects.unstructured_mesh import SemiStructuredMesh
from yt.utilities.nodal_data_utils import get_nodal_data

def _sample_ray(ray, npoints, field):
    """
    Private function that uses a ray object for calculating the field values
    that will be the y-axis values in a LinePlot object.

    Parameters
    ----------
    ray : YTOrthoRay, YTRay, or LightRay
        Ray object from which to sample field values
    npoints : int
        The number of points to sample
    field : str or field tuple
        The name of the field to sample
    """
    start_point = ray.start_point
    end_point = ray.end_point
    sample_dr = (end_point - start_point)/(npoints-1)
    sample_points = [np.arange(npoints)*sample_dr[i] for i in range(3)]
    sample_points = uvstack(sample_points).T + start_point
    ray_coordinates = uvstack([ray[d] for d in 'xyz']).T
    ray_dds = uvstack([ray['d'+d] for d in 'xyz']).T
    ray_field = ray[field]
    field_values = ray.ds.arr(np.zeros(npoints), ray_field.units)
    for i, sample_point in enumerate(sample_points):
        ray_contains = ((sample_point >= (ray_coordinates - ray_dds/2)) &
                        (sample_point <= (ray_coordinates + ray_dds/2)))
        ray_contains = ray_contains.all(axis=-1)
        # use argmax to find the first nonzero index, sometimes there
        # are two indices if the sampling point happens to fall exactly at
        # a cell boundary
        field_values[i] = ray_field[np.argmax(ray_contains)]
    dr = np.sqrt((sample_dr**2).sum())
    x = np.arange(npoints)/(npoints-1)*(dr*npoints)
    return x, field_values

class CartesianCoordinateHandler(CoordinateHandler):
    name = "cartesian"

    def __init__(self, ds, ordering = ('x','y','z')):
        super(CartesianCoordinateHandler, self).__init__(ds, ordering)

    def setup_fields(self, registry):
        for axi, ax in enumerate(self.axis_order):
            f1, f2 = _get_coord_fields(axi)
            registry.add_field(("index", "d%s" % ax), sampling_type="cell",  function = f1,
                               display_field = False,
                               units = "code_length")
            registry.add_field(("index", "path_element_%s" % ax), sampling_type="cell",  function = f1,
                               display_field = False,
                               units = "code_length")
            registry.add_field(("index", "%s" % ax), sampling_type="cell",  function = f2,
                               display_field = False,
                               units = "code_length")
            f3 = _get_vert_fields(axi)
            registry.add_field(("index", "vertex_%s" % ax), sampling_type="cell",  function = f3,
                               display_field = False,
                               units = "code_length")
        def _cell_volume(field, data):
            rv  = data["index", "dx"].copy(order='K')
            rv *= data["index", "dy"]
            rv *= data["index", "dz"]
            return rv
        registry.add_field(("index", "cell_volume"), sampling_type="cell",  function=_cell_volume,
                           display_field=False, units = "code_length**3")
        registry.check_derived_fields(
            [("index", "dx"), ("index", "dy"), ("index", "dz"),
             ("index", "x"), ("index", "y"), ("index", "z"),
             ("index", "cell_volume")])

    def pixelize(self, dimension, data_source, field, bounds, size,
                 antialias = True, periodic = True):
        """
        Method for pixelizing datasets in preparation for
        two-dimensional image plots. Relies on several sampling
        routines written in cython
        """
        index = data_source.ds.index
        if (hasattr(index, 'meshes') and
           not isinstance(index.meshes[0], SemiStructuredMesh)):
            ftype, fname = field
            if ftype == "all":
                mesh_id = 0
                indices = np.concatenate([mesh.connectivity_indices for mesh in index.mesh_union])
            else:
                mesh_id = int(ftype[-1]) - 1
                indices = index.meshes[mesh_id].connectivity_indices

            coords = index.meshes[mesh_id].connectivity_coords
            offset = index.meshes[mesh_id]._index_offset
            ad = data_source.ds.all_data()
            field_data = ad[field]
            buff_size = size[0:dimension] + (1,) + size[dimension:]

            ax = data_source.axis
            xax = self.x_axis[ax]
            yax = self.y_axis[ax]
            c = np.float64(data_source.center[dimension].d)

            extents = np.zeros((3, 2))
            extents[ax] = np.array([c, c])
            extents[xax] = bounds[0:2]
            extents[yax] = bounds[2:4]

            # if this is an element field, promote to 2D here
            if len(field_data.shape) == 1:
                field_data = np.expand_dims(field_data, 1)
            # if this is a higher-order element, we demote to 1st order
            # here, for now.
            elif field_data.shape[1] == 27:
                # hexahedral
                mylog.warning("High order elements not yet supported, " +
                              "dropping to 1st order.")
                field_data = field_data[:, 0:8]
                indices = indices[:, 0:8]

            img = pixelize_element_mesh(coords,
                                        indices,
                                        buff_size, field_data, extents,
                                        index_offset=offset)

            # re-order the array and squeeze out the dummy dim
            return np.squeeze(np.transpose(img, (yax, xax, ax)))

        elif self.axis_id.get(dimension, dimension) < 3:
            return self._ortho_pixelize(data_source, field, bounds, size,
                                        antialias, dimension, periodic)
        else:
            return self._oblique_pixelize(data_source, field, bounds, size,
                                          antialias)


    def pixelize_line(self, field, start_point, end_point, npoints):
        """
        Method for sampling datasets along a line in preparation for
        one-dimensional line plots. For UnstructuredMesh, relies on a
        sampling routine written in cython
        """
        if npoints < 2:
            raise ValueError("Must have at least two sample points in order "
                             "to draw a line plot.")
        index = self.ds.index
        if (hasattr(index, 'meshes') and
           not isinstance(index.meshes[0], SemiStructuredMesh)):
            ftype, fname = field
            if ftype == "all":
                mesh_id = 0
                indices = np.concatenate([mesh.connectivity_indices for mesh in index.mesh_union])
            else:
                mesh_id = int(ftype[-1]) - 1
                indices = index.meshes[mesh_id].connectivity_indices

            coords = index.meshes[mesh_id].connectivity_coords
            if coords.shape[1] != end_point.size != start_point.size:
                raise ValueError("The coordinate dimension doesn't match the "
                                 "start and end point dimensions.")


            offset = index.meshes[mesh_id]._index_offset
            ad = self.ds.all_data()
            field_data = ad[field]
            if field_data.shape[1] == 27:
                # hexahedral
                mylog.warning("High order elements not yet supported, " +
                              "dropping to 1st order.")
                field_data = field_data[:, 0:8]
                indices = indices[:, 0:8]

            arc_length, plot_values = pixelize_element_mesh_line(coords, indices,
                                                                 start_point,
                                                                 end_point,
                                                                 npoints, field_data,
                                                                 index_offset=offset)
            arc_length = YTArray(arc_length, start_point.units)
            plot_values = YTArray(plot_values, field_data.units)
        else:
            ray = self.ds.ray(start_point, end_point)
            arc_length, plot_values = _sample_ray(ray, npoints, field)
        return arc_length, plot_values



    def _ortho_pixelize(self, data_source, field, bounds, size, antialias,
                        dim, periodic):
        # We should be using fcoords
        period = self.period[:2].copy() # dummy here
        period[0] = self.period[self.x_axis[dim]]
        period[1] = self.period[self.y_axis[dim]]
        if hasattr(period, 'in_units'):
            period = period.in_units("code_length").d

        buff = np.zeros((size[1], size[0]), dtype="f8")

        finfo = self.ds._get_field_info(field)
        nodal_flag = finfo.nodal_flag
        if np.any(nodal_flag):
            nodal_data = get_nodal_data(data_source, field)
            coord = data_source.coord.d
            pixelize_cartesian_nodal(buff,
                                     data_source['px'], data_source['py'], data_source['pz'],
                                     data_source['pdx'], data_source['pdy'], data_source['pdz'],
                                     nodal_data, coord, bounds, int(antialias),
                                     period, int(periodic))
        else:
            pixelize_cartesian(buff, data_source['px'], data_source['py'],
                               data_source['pdx'], data_source['pdy'],
                               data_source[field],
                               bounds, int(antialias),
                               period, int(periodic))
        return buff

    def _oblique_pixelize(self, data_source, field, bounds, size, antialias):
        indices = np.argsort(data_source['pdx'])[::-1].astype(np.int_)
        buff = np.zeros((size[1], size[0]), dtype="f8")
        pixelize_off_axis_cartesian(buff,
                              data_source['index', 'x'], data_source['index', 'y'],
                              data_source['index', 'z'], data_source['px'],
                              data_source['py'], data_source['pdx'],
                              data_source['pdy'], data_source['pdz'],
                              data_source.center, data_source._inv_mat, indices,
                              data_source[field], bounds)
        return buff

    def convert_from_cartesian(self, coord):
        return coord

    def convert_to_cartesian(self, coord):
        return coord

    def convert_to_cylindrical(self, coord):
        center = self.ds.domain_center
        return cartesian_to_cylindrical(coord, center)

    def convert_from_cylindrical(self, coord):
        center = self.ds.domain_center
        return cylindrical_to_cartesian(coord, center)

    def convert_to_spherical(self, coord):
        raise NotImplementedError

    def convert_from_spherical(self, coord):
        raise NotImplementedError

    _x_pairs = (('x', 'y'), ('y', 'z'), ('z', 'x'))
    _y_pairs = (('x', 'z'), ('y', 'x'), ('z', 'y'))

    @property
    def period(self):
        return self.ds.domain_width
