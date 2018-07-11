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
    pixelize_cartesian, \
    pixelize_cartesian_nodal, \
    pixelize_sph_kernel_slice, \
    pixelize_sph_kernel_projection, \
    pixelize_element_mesh_line, \
    pixelize_sph_kernel_gather_arbitrary_grid
from yt.data_objects.unstructured_mesh import SemiStructuredMesh
from yt.utilities.nodal_data_utils import get_nodal_data
from yt.utilities.lib.particle_kdtree_tools import \
    generate_nn_list, generate_nn_list_guess
from yt.extern.tqdm import tqdm

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
            registry.add_field(("index", "d%s" % ax),
                               sampling_type="cell",
                               function = f1,
                               display_field = False,
                               units = "code_length")

            registry.add_field(("index", "path_element_%s" % ax),
                               sampling_type="cell",
                               function = f1,
                               display_field = False,
                               units = "code_length")

            registry.add_field(("index", "%s" % ax),
                               sampling_type="cell",
                               function = f2,
                               display_field = False,
                               units = "code_length")

            f3 = _get_vert_fields(axi)
            registry.add_field(("index", "vertex_%s" % ax),
                               sampling_type="cell",
                               function = f3,
                               display_field = False,
                               units = "code_length")
        def _cell_volume(field, data):
            rv  = data["index", "dx"].copy(order='K')
            rv *= data["index", "dy"]
            rv *= data["index", "dz"]
            return rv

        registry.add_field(("index", "cell_volume"),
                           sampling_type="cell",
                           function=_cell_volume,
                           display_field=False,
                           units = "code_length**3")
        registry.alias(('index', 'volume'), ('index', 'cell_volume'))

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
        from yt.frontends.sph.data_structures import ParticleDataset
        from yt.frontends.stream.data_structures import StreamParticlesDataset
        from yt.data_objects.selection_data_containers import \
            YTSlice
        from yt.data_objects.construction_data_containers import \
            YTKDTreeProj
        # We should be using fcoords
        field = data_source._determine_fields(field)[0]
        period = self.period[:2].copy() # dummy here
        period[0] = self.period[self.x_axis[dim]]
        period[1] = self.period[self.y_axis[dim]]
        if hasattr(period, 'in_units'):
            period = period.in_units("code_length").d

        buff = np.zeros((size[1], size[0]), dtype="f8")
        particle_datasets = (ParticleDataset, StreamParticlesDataset)
        is_sph_field = field[0] in 'gas'
        if hasattr(data_source.ds, '_sph_ptype'):
            is_sph_field |= field[0] in data_source.ds._sph_ptype

        finfo = self.ds._get_field_info(field)
        if np.any(finfo.nodal_flag):
            nodal_data = get_nodal_data(data_source, field)
            coord = data_source.coord.d
            pixelize_cartesian_nodal(buff,
                                     data_source['px'], data_source['py'], data_source['pz'],
                                     data_source['pdx'], data_source['pdy'], data_source['pdz'],
                                     nodal_data, coord, bounds, int(antialias),
                                     period, int(periodic))
        elif isinstance(data_source.ds, particle_datasets) and is_sph_field:
            ptype = data_source.ds._sph_ptype
            ounits = data_source.ds.field_info[field].output_units
            px_name = 'particle_position_%s' % self.axis_name[self.x_axis[dim]]
            py_name = 'particle_position_%s' % self.axis_name[self.y_axis[dim]]
            if isinstance(data_source, YTKDTreeProj):
                weight = data_source.weight_field
                le = data_source.data_source.left_edge.in_units('code_length')
                re = data_source.data_source.right_edge.in_units('code_length')
                le[self.x_axis[dim]] = bounds[0]
                le[self.y_axis[dim]] = bounds[2]
                re[self.x_axis[dim]] = bounds[1]
                re[self.y_axis[dim]] = bounds[3]
                proj_reg = data_source.ds.region(
                    left_edge=le, right_edge=re, center=data_source.center,
                    data_source=data_source.data_source
                )
                bnds = data_source.ds.arr(
                    bounds, 'code_length').in_units('cm').tolist()
                buff = np.zeros(size, dtype='float64')
                if weight is None:
                    for chunk in proj_reg.chunks([], 'io'):
                        data_source._initialize_projected_units([field], chunk)
                        pixelize_sph_kernel_projection(
                            buff,
                            chunk[ptype, px_name].in_units('cm'),
                            chunk[ptype, py_name].in_units('cm'),
                            chunk[ptype, 'smoothing_length'].in_units('cm'),
                            chunk[ptype, 'particle_mass'].in_units('g'),
                            chunk[ptype, 'density'].in_units('g/cm**3'),
                            chunk[field].in_units(ounits),
                            bnds)
                # if there is a weight field, take two projections:
                # one of field*weight, the other of just weight, and divide them
                else:
                    weight_buff = np.zeros(size, dtype='float64')
                    wounits = data_source.ds.field_info[weight].output_units
                    for chunk in proj_reg.chunks([], 'io'):
                        data_source._initialize_projected_units([field], chunk)
                        data_source._initialize_projected_units([weight], chunk)
                        pixelize_sph_kernel_projection(
                            buff,
                            chunk[ptype, px_name].in_units('cm'),
                            chunk[ptype, py_name].in_units('cm'),
                            chunk[ptype, 'smoothing_length'].in_units('cm'),
                            chunk[ptype, 'particle_mass'].in_units('g'),
                            chunk[ptype, 'density'].in_units('g/cm**3'),
                            chunk[field].in_units(ounits),
                            bnds, kernel_name="cubic",
                            weight_field=chunk[weight].in_units(wounits))
                    mylog.info("Making a fixed resolution buffer of (%s) %d by %d" % \
                                (weight, size[0], size[1]))
                    for chunk in proj_reg.chunks([], 'io'):
                        data_source._initialize_projected_units([weight], chunk)
                        pixelize_sph_kernel_projection(
                            weight_buff,
                            chunk[ptype, px_name].in_units('cm'),
                            chunk[ptype, py_name].in_units('cm'),
                            chunk[ptype, 'smoothing_length'].in_units('cm'),
                            chunk[ptype, 'particle_mass'].in_units('g'),
                            chunk[ptype, 'density'].in_units('g/cm**3'),
                            chunk[weight].in_units(wounits),
                            bnds)
                    buff /= weight_buff
            elif isinstance(data_source, YTSlice):
                smoothing_style = "scatter"
                if hasattr(data_source.ds, 'sph_smoothing_style'):
                    smoothing_style = data_source.ds.sph_smoothing_style

                buff = np.zeros(size, dtype='float64')
                if smoothing_style == "scatter":
                    for chunk in data_source.chunks([], 'io'):
                        pixelize_sph_kernel_slice(
                            buff,
                            chunk[ptype, px_name],
                            chunk[ptype, py_name],
                            chunk[ptype, 'smoothing_length'],
                            chunk[ptype, 'particle_mass'],
                            chunk[ptype, 'density'],
                            chunk[field].in_units(ounits),
                            bounds)
                if smoothing_style == "gather":
                    tree = data_source.ds.index.kdtree
                    grid_size = np.array([size[0], size[1], 1], dtype="int64")
                    grid_bounds = np.array([0, 10000, 0, 10000, 0, 10000],
                                           dtype="float64")
                    # set up the nearest neighbour arrays to store the distances
                    # and the particle ids
                    pids = np.zeros((size[0], size[1], 1,
                                 self.ds._num_neighbors), dtype="int64") - 1
                    dists = np.zeros((size[0], size[1], 1,
                                 self.ds._num_neighbors), dtype="float64") - 1
                    queue_sizes = np.zeros(grid_size, dtype="int64")
                    dest_num = np.zeros(grid_size, dtype="float64")

                    # find all the nearest neighbors, we have to loop through
                    # EVERY SINGLE PARTICLE and fill our nearest neighbour lists
                    # before we can do the deposition
                    offsets = np.array([0], dtype="int64")

                    # first we loop through and fill the particles with our best
                    # guess
                    pbar = tqdm(desc="Generating nearest neighbor lists")
                    for i, chunk in enumerate(
                            self.ds.all_data().chunks([field], 'io')):
                        offsets = np.append(offsets,
                                offsets[-1]+chunk[(ptype,'density')].shape[0])

                        generate_nn_list_guess(pids, dists, queue_sizes,
                            tree, grid_bounds, grid_size,
                            chunk[(ptype,'particle_position')].in_base("code").d,
                            offset=offsets[i], gather_type=1)
                        pbar.update(i)

                    # now loop through and traverse the tree. It is much quicker
                    # to do it this way, as the better the initial guess the faster
                    # the main process
                    for i, chunk in enumerate(
                            self.ds.all_data().chunks([field], 'io')):
                        generate_nn_list(pids, dists, tree, grid_bounds, grid_size,
                            chunk[(ptype,'particle_position')].in_base("code").d,
                            offset=offsets[i], gather_type=1)
                        pbar.update(i+offsets.shape[0])
                    pbar.close()

                    # perform the deposition onto the pixels -> do it twice to
                    # allow normalization
                    pbar = tqdm(desc="Interpolating SPH field {}".format(field))
                    for i, chunk in enumerate(
                            data_source.ds.all_data().chunks([field], 'io')):
                        pixelize_sph_kernel_gather_arbitrary_grid(dest_num, pids,
                            dists,
                            chunk[(ptype,'smoothing_length')].in_base("code").d,
                            chunk[(ptype,'particle_mass')].in_base("code").d,
                            chunk[(ptype,'density')].in_base("code").d,
                            chunk[field].d, tree=tree, offset=offsets[i],
                            use_normalization=False)
                        pbar.update(i)
                    pbar.close()
                    buff[:, :] = dest_num[:, :, 0]
            else:
                raise NotImplementedError(
                    "A pixelization routine has not been implemented for %s "
                    "data objects" % str(type(data_source)))
            buff = buff.transpose()
        else:
            pixelize_cartesian(buff,
                               data_source['px'], data_source['py'],
                               data_source['pdx'], data_source['pdy'],
                               data_source[field],
                               bounds, int(antialias),
                               period, int(periodic))

        return buff

    def _oblique_pixelize(self, data_source, field, bounds, size, antialias):
        indices = np.argsort(data_source['pdx'])[::-1].astype(np.int_)
        buff = np.zeros((size[1], size[0]), dtype="f8")
        pixelize_off_axis_cartesian(buff,
                              data_source['x'], data_source['y'],
                              data_source['z'], data_source['px'],
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
