import numpy as np

from yt.data_objects.index_subobjects.unstructured_mesh import SemiStructuredMesh
from yt.funcs import mylog
from yt.units._numpy_wrapper_functions import uconcatenate, uvstack
from yt.units.yt_array import YTArray
from yt.utilities.lib.pixelization_routines import (
    interpolate_sph_grid_gather,
    normalization_2d_utility,
    pixelize_cartesian,
    pixelize_cartesian_nodal,
    pixelize_element_mesh,
    pixelize_element_mesh_line,
    pixelize_off_axis_cartesian,
    pixelize_sph_kernel_projection,
    pixelize_sph_kernel_slice,
)
from yt.utilities.math_utils import compute_stddev_image
from yt.utilities.nodal_data_utils import get_nodal_data

from .coordinate_handler import (
    CoordinateHandler,
    _get_coord_fields,
    _get_vert_fields,
    cartesian_to_cylindrical,
    cylindrical_to_cartesian,
)


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
    sample_dr = (end_point - start_point) / (npoints - 1)
    sample_points = [np.arange(npoints) * sample_dr[i] for i in range(3)]
    sample_points = uvstack(sample_points).T + start_point
    ray_coordinates = uvstack([ray[("index", d)] for d in "xyz"]).T
    ray_dds = uvstack([ray[("index", f"d{d}")] for d in "xyz"]).T
    ray_field = ray[field]
    field_values = ray.ds.arr(np.zeros(npoints), ray_field.units)
    for i, sample_point in enumerate(sample_points):
        ray_contains = (sample_point >= (ray_coordinates - ray_dds / 2)) & (
            sample_point <= (ray_coordinates + ray_dds / 2)
        )
        ray_contains = ray_contains.all(axis=-1)
        # use argmax to find the first nonzero index, sometimes there
        # are two indices if the sampling point happens to fall exactly at
        # a cell boundary
        field_values[i] = ray_field[np.argmax(ray_contains)]
    dr = np.sqrt((sample_dr**2).sum())
    x = np.arange(npoints) / (npoints - 1) * (dr * npoints)
    return x, field_values


def all_data(data, ptype, fields, kdtree=False):
    field_data = {}
    fields = set(fields)
    for field in fields:
        field_data[field] = []

    for chunk in data.all_data().chunks([], "io"):
        for field in fields:
            field_data[field].append(chunk[ptype, field].in_base("code"))

    for field in fields:
        field_data[field] = uconcatenate(field_data[field])

    if kdtree is True:
        kdtree = data.index.kdtree
        for field in fields:
            if len(field_data[field].shape) == 1:
                field_data[field] = field_data[field][kdtree.idx]
            else:
                field_data[field] = field_data[field][kdtree.idx, :]

    return field_data


class CartesianCoordinateHandler(CoordinateHandler):
    name = "cartesian"

    def __init__(self, ds, ordering=("x", "y", "z")):
        super().__init__(ds, ordering)

    def setup_fields(self, registry):
        for axi, ax in enumerate(self.axis_order):
            f1, f2 = _get_coord_fields(axi)
            registry.add_field(
                ("index", f"d{ax}"),
                sampling_type="cell",
                function=f1,
                display_field=False,
                units="code_length",
            )

            registry.add_field(
                ("index", f"path_element_{ax}"),
                sampling_type="cell",
                function=f1,
                display_field=False,
                units="code_length",
            )

            registry.add_field(
                ("index", f"{ax}"),
                sampling_type="cell",
                function=f2,
                display_field=False,
                units="code_length",
            )

            f3 = _get_vert_fields(axi)
            registry.add_field(
                ("index", f"vertex_{ax}"),
                sampling_type="cell",
                function=f3,
                display_field=False,
                units="code_length",
            )

        self._register_volume(registry)
        self._check_fields(registry)

    def _register_volume(self, registry):
        def _cell_volume(field, data):
            rv = data["index", "dx"].copy(order="K")
            rv *= data["index", "dy"]
            rv *= data["index", "dz"]
            return rv

        registry.add_field(
            ("index", "cell_volume"),
            sampling_type="cell",
            function=_cell_volume,
            display_field=False,
            units="code_length**3",
        )
        registry.alias(("index", "volume"), ("index", "cell_volume"))

    def _check_fields(self, registry):
        registry.check_derived_fields(
            [
                ("index", "dx"),
                ("index", "dy"),
                ("index", "dz"),
                ("index", "x"),
                ("index", "y"),
                ("index", "z"),
                ("index", "cell_volume"),
            ]
        )

    def pixelize(
        self, dimension, data_source, field, bounds, size, antialias=True, periodic=True
    ):
        """
        Method for pixelizing datasets in preparation for
        two-dimensional image plots. Relies on several sampling
        routines written in cython
        """
        index = data_source.ds.index
        if hasattr(index, "meshes") and not isinstance(
            index.meshes[0], SemiStructuredMesh
        ):
            ftype, fname = field
            if ftype == "all":
                mesh_id = 0
                indices = np.concatenate(
                    [mesh.connectivity_indices for mesh in index.mesh_union]
                )
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
                mylog.warning(
                    "High order elements not yet supported, dropping to 1st order."
                )
                field_data = field_data[:, 0:8]
                indices = indices[:, 0:8]

            img = pixelize_element_mesh(
                coords, indices, buff_size, field_data, extents, index_offset=offset
            )

            # re-order the array and squeeze out the dummy dim
            return np.squeeze(np.transpose(img, (yax, xax, ax)))

        elif self.axis_id.get(dimension, dimension) < 3:
            return self._ortho_pixelize(
                data_source, field, bounds, size, antialias, dimension, periodic
            )
        else:
            return self._oblique_pixelize(data_source, field, bounds, size, antialias)

    def pixelize_line(self, field, start_point, end_point, npoints):
        """
        Method for sampling datasets along a line in preparation for
        one-dimensional line plots. For UnstructuredMesh, relies on a
        sampling routine written in cython
        """
        if npoints < 2:
            raise ValueError(
                "Must have at least two sample points in order to draw a line plot."
            )
        index = self.ds.index
        if hasattr(index, "meshes") and not isinstance(
            index.meshes[0], SemiStructuredMesh
        ):
            ftype, fname = field
            if ftype == "all":
                mesh_id = 0
                indices = np.concatenate(
                    [mesh.connectivity_indices for mesh in index.mesh_union]
                )
            else:
                mesh_id = int(ftype[-1]) - 1
                indices = index.meshes[mesh_id].connectivity_indices

            coords = index.meshes[mesh_id].connectivity_coords
            if coords.shape[1] != end_point.size != start_point.size:
                raise ValueError(
                    "The coordinate dimension doesn't match the "
                    "start and end point dimensions."
                )

            offset = index.meshes[mesh_id]._index_offset
            ad = self.ds.all_data()
            field_data = ad[field]

            # if this is an element field, promote to 2D here
            if len(field_data.shape) == 1:
                field_data = np.expand_dims(field_data, 1)
            # if this is a higher-order element, we demote to 1st order
            # here, for now.
            elif field_data.shape[1] == 27:
                # hexahedral
                mylog.warning(
                    "High order elements not yet supported, dropping to 1st order."
                )
                field_data = field_data[:, 0:8]
                indices = indices[:, 0:8]

            arc_length, plot_values = pixelize_element_mesh_line(
                coords,
                indices,
                start_point,
                end_point,
                npoints,
                field_data,
                index_offset=offset,
            )
            arc_length = YTArray(arc_length, start_point.units)
            plot_values = YTArray(plot_values, field_data.units)
        else:
            ray = self.ds.ray(start_point, end_point)
            arc_length, plot_values = _sample_ray(ray, npoints, field)
        return arc_length, plot_values

    def _ortho_pixelize(
        self, data_source, field, bounds, size, antialias, dim, periodic
    ):
        from yt.data_objects.construction_data_containers import YTParticleProj
        from yt.data_objects.selection_objects.slices import YTSlice
        from yt.frontends.sph.data_structures import ParticleDataset
        from yt.frontends.stream.data_structures import StreamParticlesDataset

        # We should be using fcoords
        field = data_source._determine_fields(field)[0]
        finfo = data_source.ds.field_info[field]
        period = self.period[:2].copy()  # dummy here
        period[0] = self.period[self.x_axis[dim]]
        period[1] = self.period[self.y_axis[dim]]
        if hasattr(period, "in_units"):
            period = period.in_units("code_length").d

        buff = np.full((size[1], size[0]), np.nan, dtype="float64")
        particle_datasets = (ParticleDataset, StreamParticlesDataset)
        is_sph_field = finfo.is_sph_field

        finfo = self.ds._get_field_info(field)
        if np.any(finfo.nodal_flag):
            nodal_data = get_nodal_data(data_source, field)
            coord = data_source.coord.d
            pixelize_cartesian_nodal(
                buff,
                data_source["px"],
                data_source["py"],
                data_source["pz"],
                data_source["pdx"],
                data_source["pdy"],
                data_source["pdz"],
                nodal_data,
                coord,
                bounds,
                int(antialias),
                period,
                int(periodic),
            )
        elif isinstance(data_source.ds, particle_datasets) and is_sph_field:
            ptype = field[0]
            if ptype == "gas":
                ptype = data_source.ds._sph_ptypes[0]
            px_name = self.axis_name[self.x_axis[dim]]
            py_name = self.axis_name[self.y_axis[dim]]
            ounits = data_source.ds.field_info[field].output_units
            bnds = data_source.ds.arr(bounds, "code_length").tolist()
            if isinstance(data_source, YTParticleProj):
                weight = data_source.weight_field
                moment = data_source.moment
                le, re = data_source.data_source.get_bbox()
                xa = self.x_axis[dim]
                ya = self.y_axis[dim]
                # If we're not periodic, we need to clip to the boundary edges
                # or we get errors about extending off the edge of the region.
                if not self.ds.periodicity[xa]:
                    le[xa] = max(bounds[0], self.ds.domain_left_edge[xa])
                    re[xa] = min(bounds[1], self.ds.domain_right_edge[xa])
                else:
                    le[xa] = bounds[0]
                    re[xa] = bounds[1]
                if not self.ds.periodicity[ya]:
                    le[ya] = max(bounds[2], self.ds.domain_left_edge[ya])
                    re[ya] = min(bounds[3], self.ds.domain_right_edge[ya])
                else:
                    le[ya] = bounds[2]
                    re[ya] = bounds[3]
                # We actually need to clip these
                proj_reg = data_source.ds.region(
                    left_edge=le,
                    right_edge=re,
                    center=data_source.center,
                    data_source=data_source.data_source,
                )
                proj_reg.set_field_parameter("axis", data_source.axis)
                buff = np.zeros(size, dtype="float64")
                if weight is None:
                    for chunk in proj_reg.chunks([], "io"):
                        data_source._initialize_projected_units([field], chunk)
                        pixelize_sph_kernel_projection(
                            buff,
                            chunk[ptype, px_name].to("code_length"),
                            chunk[ptype, py_name].to("code_length"),
                            chunk[ptype, "smoothing_length"].to("code_length"),
                            chunk[ptype, "mass"].to("code_mass"),
                            chunk[ptype, "density"].to("code_density"),
                            chunk[field].in_units(ounits),
                            bnds,
                            check_period=int(periodic),
                            period=period,
                        )
                    # We use code length here, but to get the path length right
                    # we need to multiply by the conversion factor between
                    # code length and the unit system's length unit
                    default_path_length_unit = data_source.ds.unit_system["length"]
                    dl_conv = data_source.ds.quan(1.0, "code_length").to(
                        default_path_length_unit
                    )
                    buff *= dl_conv.v
                # if there is a weight field, take two projections:
                # one of field*weight, the other of just weight, and divide them
                else:
                    weight_buff = np.zeros(size, dtype="float64")
                    buff = np.zeros(size, dtype="float64")
                    wounits = data_source.ds.field_info[weight].output_units
                    for chunk in proj_reg.chunks([], "io"):
                        data_source._initialize_projected_units([field], chunk)
                        data_source._initialize_projected_units([weight], chunk)
                        pixelize_sph_kernel_projection(
                            buff,
                            chunk[ptype, px_name].to("code_length"),
                            chunk[ptype, py_name].to("code_length"),
                            chunk[ptype, "smoothing_length"].to("code_length"),
                            chunk[ptype, "mass"].to("code_mass"),
                            chunk[ptype, "density"].to("code_density"),
                            chunk[field].in_units(ounits),
                            bnds,
                            check_period=int(periodic),
                            period=period,
                            weight_field=chunk[weight].in_units(wounits),
                        )
                    mylog.info(
                        "Making a fixed resolution buffer of (%s) %d by %d",
                        weight,
                        size[0],
                        size[1],
                    )
                    for chunk in proj_reg.chunks([], "io"):
                        data_source._initialize_projected_units([weight], chunk)
                        pixelize_sph_kernel_projection(
                            weight_buff,
                            chunk[ptype, px_name].to("code_length"),
                            chunk[ptype, py_name].to("code_length"),
                            chunk[ptype, "smoothing_length"].to("code_length"),
                            chunk[ptype, "mass"].to("code_mass"),
                            chunk[ptype, "density"].to("code_density"),
                            chunk[weight].in_units(wounits),
                            bnds,
                            check_period=int(periodic),
                            period=period,
                        )
                    normalization_2d_utility(buff, weight_buff)
                    if moment == 2:
                        buff2 = np.zeros(size, dtype="float64")
                        for chunk in proj_reg.chunks([], "io"):
                            data_source._initialize_projected_units([field], chunk)
                            data_source._initialize_projected_units([weight], chunk)
                            pixelize_sph_kernel_projection(
                                buff2,
                                chunk[ptype, px_name].to("code_length"),
                                chunk[ptype, py_name].to("code_length"),
                                chunk[ptype, "smoothing_length"].to("code_length"),
                                chunk[ptype, "mass"].to("code_mass"),
                                chunk[ptype, "density"].to("code_density"),
                                chunk[field].in_units(ounits) ** 2,
                                bnds,
                                check_period=int(periodic),
                                period=period,
                                weight_field=chunk[weight].in_units(wounits),
                            )
                        normalization_2d_utility(buff2, weight_buff)
                        buff = compute_stddev_image(buff2, buff)
            elif isinstance(data_source, YTSlice):
                smoothing_style = getattr(self.ds, "sph_smoothing_style", "scatter")
                normalize = getattr(self.ds, "use_sph_normalization", True)

                if smoothing_style == "scatter":
                    buff = np.zeros(size, dtype="float64")
                    if normalize:
                        buff_den = np.zeros(size, dtype="float64")

                    for chunk in data_source.chunks([], "io"):
                        pixelize_sph_kernel_slice(
                            buff,
                            chunk[ptype, px_name].to("code_length"),
                            chunk[ptype, py_name].to("code_length"),
                            chunk[ptype, "smoothing_length"].to("code_length"),
                            chunk[ptype, "mass"].to("code_mass"),
                            chunk[ptype, "density"].to("code_density"),
                            chunk[field].in_units(ounits),
                            bnds,
                            check_period=int(periodic),
                            period=period,
                        )
                        if normalize:
                            pixelize_sph_kernel_slice(
                                buff_den,
                                chunk[ptype, px_name].to("code_length"),
                                chunk[ptype, py_name].to("code_length"),
                                chunk[ptype, "smoothing_length"].to("code_length"),
                                chunk[ptype, "mass"].to("code_mass"),
                                chunk[ptype, "density"].to("code_density"),
                                np.ones(chunk[ptype, "density"].shape[0]),
                                bnds,
                                check_period=int(periodic),
                                period=period,
                            )

                    if normalize:
                        normalization_2d_utility(buff, buff_den)

                if smoothing_style == "gather":
                    # Here we find out which axis are going to be the "x" and
                    # "y" axis for the actual visualisation and then we set the
                    # buffer size and bounds to match. The z axis of the plot
                    # is the axis we slice over and the buffer will be of size 1
                    # in that dimension
                    x, y, z = self.x_axis[dim], self.y_axis[dim], dim

                    buff_size = np.zeros(3, dtype="int64")
                    buff_size[x] = size[0]
                    buff_size[y] = size[1]
                    buff_size[z] = 1

                    buff_bounds = np.zeros(6, dtype="float64")
                    buff_bounds[2 * x : 2 * x + 2] = bounds[0:2]
                    buff_bounds[2 * y : 2 * y + 2] = bounds[2:4]
                    buff_bounds[2 * z] = data_source.coord
                    buff_bounds[2 * z + 1] = data_source.coord

                    # then we do the interpolation
                    buff_temp = np.zeros(buff_size, dtype="float64")

                    fields_to_get = [
                        "particle_position",
                        "density",
                        "mass",
                        "smoothing_length",
                        field[1],
                    ]
                    all_fields = all_data(self.ds, ptype, fields_to_get, kdtree=True)

                    num_neighbors = getattr(self.ds, "num_neighbors", 32)
                    interpolate_sph_grid_gather(
                        buff_temp,
                        all_fields["particle_position"].to("code_length"),
                        buff_bounds,
                        all_fields["smoothing_length"].to("code_length"),
                        all_fields["mass"].to("code_mass"),
                        all_fields["density"].to("code_density"),
                        all_fields[field[1]].in_units(ounits),
                        self.ds.index.kdtree,
                        num_neigh=num_neighbors,
                        use_normalization=normalize,
                    )

                    # We swap the axes back so the axis which was sliced over
                    # is the last axis, as this is the "z" axis of the plots.
                    if z != 2:
                        buff_temp = buff_temp.swapaxes(2, z)
                        if x == 2:
                            x = z
                        else:
                            y = z

                    buff = buff_temp[:, :, 0]

                    # Then we just transpose if the buffer x and y are
                    # different than the plot x and y
                    if y < x:
                        buff = buff.transpose()
            else:
                raise NotImplementedError(
                    "A pixelization routine has not been implemented for %s "
                    "data objects" % str(type(data_source))
                )
            buff = buff.transpose()
        else:
            pixelize_cartesian(
                buff,
                data_source["px"],
                data_source["py"],
                data_source["pdx"],
                data_source["pdy"],
                data_source[field],
                bounds,
                int(antialias),
                period,
                int(periodic),
            )
        return buff

    def _oblique_pixelize(self, data_source, field, bounds, size, antialias):
        from yt.frontends.ytdata.data_structures import YTSpatialPlotDataset

        indices = np.argsort(data_source["pdx"])[::-1].astype(np.int_)
        buff = np.full((size[1], size[0]), np.nan, dtype="float64")
        ftype = "index"
        if isinstance(data_source.ds, YTSpatialPlotDataset):
            ftype = "gas"
        pixelize_off_axis_cartesian(
            buff,
            data_source[ftype, "x"],
            data_source[ftype, "y"],
            data_source[ftype, "z"],
            data_source["px"],
            data_source["py"],
            data_source["pdx"],
            data_source["pdy"],
            data_source["pdz"],
            data_source.center,
            data_source._inv_mat,
            indices,
            data_source[field],
            bounds,
        )
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

    _x_pairs = (("x", "y"), ("y", "z"), ("z", "x"))
    _y_pairs = (("x", "z"), ("y", "x"), ("z", "y"))

    @property
    def period(self):
        return self.ds.domain_width
