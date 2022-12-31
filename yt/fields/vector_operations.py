import numpy as np

from yt.funcs import is_sequence, just_one
from yt.geometry.geometry_handler import is_curvilinear
from yt.utilities.lib.misc_utilities import obtain_relative_velocity_vector
from yt.utilities.math_utils import (
    get_cyl_r_component,
    get_cyl_theta_component,
    get_cyl_z_component,
    get_sph_phi_component,
    get_sph_r_component,
    get_sph_theta_component,
)

from .derived_field import NeedsParameter, ValidateParameter, ValidateSpatial


def get_bulk(data, basename, unit):
    if data.has_field_parameter(f"bulk_{basename}"):
        bulk = data.get_field_parameter(f"bulk_{basename}")
    else:
        bulk = [0, 0, 0] * unit
    return bulk


def create_magnitude_field(
    registry,
    basename,
    field_units,
    ftype="gas",
    slice_info=None,
    validators=None,
    sampling_type=None,
):

    axis_order = registry.ds.coordinates.axis_order

    field_components = [(ftype, f"{basename}_{ax}") for ax in axis_order]

    if sampling_type is None:
        sampling_type = "local"

    def _magnitude(field, data):
        fn = field_components[0]
        if data.has_field_parameter(f"bulk_{basename}"):
            fn = (fn[0], f"relative_{fn[1]}")
        d = data[fn]
        mag = (d) ** 2
        for idim in range(1, registry.ds.dimensionality):
            fn = field_components[idim]
            if data.has_field_parameter(f"bulk_{basename}"):
                fn = (fn[0], f"relative_{fn[1]}")
            mag += (data[fn]) ** 2
        return np.sqrt(mag)

    registry.add_field(
        (ftype, f"{basename}_magnitude"),
        sampling_type=sampling_type,
        function=_magnitude,
        units=field_units,
        validators=validators,
    )


def create_relative_field(
    registry, basename, field_units, ftype="gas", slice_info=None, validators=None
):

    axis_order = registry.ds.coordinates.axis_order

    field_components = [(ftype, f"{basename}_{ax}") for ax in axis_order]

    def relative_vector(ax):
        def _relative_vector(field, data):
            iax = axis_order.index(ax)
            d = data[field_components[iax]]
            bulk = get_bulk(data, basename, d.unit_quantity)
            return d - bulk[iax]

        return _relative_vector

    for d in axis_order:
        registry.add_field(
            (ftype, f"relative_{basename}_{d}"),
            sampling_type="local",
            function=relative_vector(d),
            units=field_units,
            validators=validators,
        )


def create_los_field(registry, basename, field_units, ftype="gas", slice_info=None):
    axis_order = registry.ds.coordinates.axis_order

    validators = [
        ValidateParameter(f"bulk_{basename}"),
        ValidateParameter("axis", {"axis": [0, 1, 2]}),
    ]

    field_comps = [(ftype, f"{basename}_{ax}") for ax in axis_order]

    def _los_field(field, data):
        if data.has_field_parameter(f"bulk_{basename}"):
            fns = [(fc[0], f"relative_{fc[1]}") for fc in field_comps]
        else:
            fns = field_comps
        ax = data.get_field_parameter("axis")
        if is_sequence(ax):
            # Make sure this is a unit vector
            ax /= np.sqrt(np.dot(ax, ax))
            ret = data[fns[0]] * ax[0] + data[fns[1]] * ax[1] + data[fns[2]] * ax[2]
        elif ax in [0, 1, 2]:
            ret = data[fns[ax]]
        else:
            raise NeedsParameter(["axis"])
        return ret

    registry.add_field(
        (ftype, f"{basename}_los"),
        sampling_type="local",
        function=_los_field,
        units=field_units,
        validators=validators,
        display_name=r"\mathrm{Line of Sight %s}" % basename.capitalize(),
    )


def create_squared_field(
    registry, basename, field_units, ftype="gas", slice_info=None, validators=None
):

    axis_order = registry.ds.coordinates.axis_order

    field_components = [(ftype, f"{basename}_{ax}") for ax in axis_order]

    def _squared(field, data):
        fn = field_components[0]
        if data.has_field_parameter(f"bulk_{basename}"):
            fn = (fn[0], f"relative_{fn[1]}")
        squared = data[fn] * data[fn]
        for idim in range(1, registry.ds.dimensionality):
            fn = field_components[idim]
            squared += data[fn] * data[fn]
        return squared

    registry.add_field(
        (ftype, f"{basename}_squared"),
        sampling_type="local",
        function=_squared,
        units=field_units,
        validators=validators,
    )


def create_vector_fields(registry, basename, field_units, ftype="gas", slice_info=None):
    from yt.units.unit_object import Unit

    # slice_info would be the left, the right, and the factor.
    # For example, with the old Enzo-ZEUS fields, this would be:
    # slice(None, -2, None)
    # slice(1, -1, None)
    # 1.0
    # Otherwise, we default to a centered difference.
    if slice_info is None:
        sl_left = slice(None, -2, None)
        sl_right = slice(2, None, None)
        div_fac = 2.0
    else:
        sl_left, sl_right, div_fac = slice_info

    axis_order = registry.ds.coordinates.axis_order

    xn, yn, zn = ((ftype, f"{basename}_{ax}") for ax in axis_order)

    # Is this safe?
    if registry.ds.dimensionality < 3:
        zn = ("index", "zeros")
    if registry.ds.dimensionality < 2:
        yn = ("index", "zeros")

    create_relative_field(
        registry,
        basename,
        field_units,
        ftype=ftype,
        slice_info=slice_info,
        validators=[ValidateParameter(f"bulk_{basename}")],
    )

    create_magnitude_field(
        registry,
        basename,
        field_units,
        ftype=ftype,
        slice_info=slice_info,
        validators=[ValidateParameter(f"bulk_{basename}")],
    )

    axis_names = registry.ds.coordinates.axis_order

    if not is_curvilinear(registry.ds.geometry):

        # The following fields are invalid for curvilinear geometries
        def _spherical_radius_component(field, data):
            """The spherical radius component of the vector field

            Relative to the coordinate system defined by the *normal* vector,
            *center*, and *bulk_* field parameters.
            """
            normal = data.get_field_parameter("normal")
            vectors = obtain_relative_velocity_vector(
                data, (xn, yn, zn), f"bulk_{basename}"
            )
            theta = data["index", "spherical_theta"]
            phi = data["index", "spherical_phi"]
            rv = get_sph_r_component(vectors, theta, phi, normal)
            # Now, anywhere that radius is in fact zero, we want to zero out our
            # return values.
            rv[np.isnan(theta)] = 0.0
            return rv

        registry.add_field(
            (ftype, f"{basename}_spherical_radius"),
            sampling_type="local",
            function=_spherical_radius_component,
            units=field_units,
            validators=[
                ValidateParameter("normal"),
                ValidateParameter("center"),
                ValidateParameter(f"bulk_{basename}"),
            ],
        )
        create_los_field(
            registry, basename, field_units, ftype=ftype, slice_info=slice_info
        )

        def _radial(field, data):
            return data[ftype, f"{basename}_spherical_radius"]

        def _radial_absolute(field, data):
            return np.abs(data[ftype, f"{basename}_spherical_radius"])

        def _tangential(field, data):
            return np.sqrt(
                data[ftype, f"{basename}_spherical_theta"] ** 2.0
                + data[ftype, f"{basename}_spherical_phi"] ** 2.0
            )

        registry.add_field(
            (ftype, f"radial_{basename}"),
            sampling_type="local",
            function=_radial,
            units=field_units,
            validators=[ValidateParameter("normal"), ValidateParameter("center")],
        )

        registry.add_field(
            (ftype, f"radial_{basename}_absolute"),
            sampling_type="local",
            function=_radial_absolute,
            units=field_units,
        )

        registry.add_field(
            (ftype, f"tangential_{basename}"),
            sampling_type="local",
            function=_tangential,
            units=field_units,
        )

        def _spherical_theta_component(field, data):
            """The spherical theta component of the vector field

            Relative to the coordinate system defined by the *normal* vector,
            *center*, and *bulk_* field parameters.
            """
            normal = data.get_field_parameter("normal")
            vectors = obtain_relative_velocity_vector(
                data, (xn, yn, zn), f"bulk_{basename}"
            )
            theta = data["index", "spherical_theta"]
            phi = data["index", "spherical_phi"]
            return get_sph_theta_component(vectors, theta, phi, normal)

        registry.add_field(
            (ftype, f"{basename}_spherical_theta"),
            sampling_type="local",
            function=_spherical_theta_component,
            units=field_units,
            validators=[
                ValidateParameter("normal"),
                ValidateParameter("center"),
                ValidateParameter(f"bulk_{basename}"),
            ],
        )

        def _spherical_phi_component(field, data):
            """The spherical phi component of the vector field

            Relative to the coordinate system defined by the *normal* vector,
            *center*, and *bulk_* field parameters.
            """
            normal = data.get_field_parameter("normal")
            vectors = obtain_relative_velocity_vector(
                data, (xn, yn, zn), f"bulk_{basename}"
            )
            phi = data["index", "spherical_phi"]
            return get_sph_phi_component(vectors, phi, normal)

        registry.add_field(
            (ftype, f"{basename}_spherical_phi"),
            sampling_type="local",
            function=_spherical_phi_component,
            units=field_units,
            validators=[
                ValidateParameter("normal"),
                ValidateParameter("center"),
                ValidateParameter(f"bulk_{basename}"),
            ],
        )

        def _cp_vectors(ax):
            def _cp_val(field, data):
                vec = data.get_field_parameter(f"cp_{ax}_vec")
                tr = data[xn[0], f"relative_{xn[1]}"] * vec.d[0]
                tr += data[yn[0], f"relative_{yn[1]}"] * vec.d[1]
                tr += data[zn[0], f"relative_{zn[1]}"] * vec.d[2]
                return tr

            return _cp_val

        for ax in "xyz":
            registry.add_field(
                (ftype, f"cutting_plane_{basename}_{ax}"),
                sampling_type="local",
                function=_cp_vectors(ax),
                units=field_units,
            )

        def _divergence(field, data):
            ds = div_fac * just_one(data["index", "dx"])
            f = data[xn[0], f"relative_{xn[1]}"][sl_right, 1:-1, 1:-1] / ds
            f -= data[xn[0], f"relative_{xn[1]}"][sl_left, 1:-1, 1:-1] / ds
            ds = div_fac * just_one(data["index", "dy"])
            f += data[yn[0], f"relative_{yn[1]}"][1:-1, sl_right, 1:-1] / ds
            f -= data[yn[0], f"relative_{yn[1]}"][1:-1, sl_left, 1:-1] / ds
            ds = div_fac * just_one(data["index", "dz"])
            f += data[zn[0], f"relative_{zn[1]}"][1:-1, 1:-1, sl_right] / ds
            f -= data[zn[0], f"relative_{zn[1]}"][1:-1, 1:-1, sl_left] / ds
            new_field = data.ds.arr(np.zeros(data[xn].shape, dtype=np.float64), f.units)
            new_field[1:-1, 1:-1, 1:-1] = f
            return new_field

        def _divergence_abs(field, data):
            return np.abs(data[ftype, f"{basename}_divergence"])

        field_units = Unit(field_units, registry=registry.ds.unit_registry)
        div_units = field_units / registry.ds.unit_system["length"]

        registry.add_field(
            (ftype, f"{basename}_divergence"),
            sampling_type="local",
            function=_divergence,
            units=div_units,
            validators=[ValidateSpatial(1), ValidateParameter(f"bulk_{basename}")],
        )

        registry.add_field(
            (ftype, f"{basename}_divergence_absolute"),
            sampling_type="local",
            function=_divergence_abs,
            units=div_units,
        )

        def _tangential_over_magnitude(field, data):
            tr = (
                data[ftype, f"tangential_{basename}"]
                / data[ftype, f"{basename}_magnitude"]
            )
            return np.abs(tr)

        registry.add_field(
            (ftype, f"tangential_over_{basename}_magnitude"),
            sampling_type="local",
            function=_tangential_over_magnitude,
            take_log=False,
        )

        def _cylindrical_radius_component(field, data):
            """The cylindrical radius component of the vector field

            Relative to the coordinate system defined by the *normal* vector,
            *center*, and *bulk_* field parameters.
            """
            normal = data.get_field_parameter("normal")
            vectors = obtain_relative_velocity_vector(
                data, (xn, yn, zn), f"bulk_{basename}"
            )
            theta = data["index", "cylindrical_theta"]
            return get_cyl_r_component(vectors, theta, normal)

        registry.add_field(
            (ftype, f"{basename}_cylindrical_radius"),
            sampling_type="local",
            function=_cylindrical_radius_component,
            units=field_units,
            validators=[ValidateParameter("normal")],
        )

        def _cylindrical_theta_component(field, data):
            """The cylindrical theta component of the vector field

            Relative to the coordinate system defined by the *normal* vector,
            *center*, and *bulk_* field parameters.
            """
            normal = data.get_field_parameter("normal")
            vectors = obtain_relative_velocity_vector(
                data, (xn, yn, zn), f"bulk_{basename}"
            )
            theta = data["index", "cylindrical_theta"].copy()
            theta = np.tile(theta, (3,) + (1,) * len(theta.shape))
            return get_cyl_theta_component(vectors, theta, normal)

        registry.add_field(
            (ftype, f"{basename}_cylindrical_theta"),
            sampling_type="local",
            function=_cylindrical_theta_component,
            units=field_units,
            validators=[
                ValidateParameter("normal"),
                ValidateParameter("center"),
                ValidateParameter(f"bulk_{basename}"),
            ],
        )

        def _cylindrical_z_component(field, data):
            """The cylindrical z component of the vector field

            Relative to the coordinate system defined by the *normal* vector,
            *center*, and *bulk_* field parameters.
            """
            normal = data.get_field_parameter("normal")
            vectors = obtain_relative_velocity_vector(
                data, (xn, yn, zn), f"bulk_{basename}"
            )
            return get_cyl_z_component(vectors, normal)

        registry.add_field(
            (ftype, f"{basename}_cylindrical_z"),
            sampling_type="local",
            function=_cylindrical_z_component,
            units=field_units,
            validators=[
                ValidateParameter("normal"),
                ValidateParameter("center"),
                ValidateParameter(f"bulk_{basename}"),
            ],
        )

    else:  # Create Cartesian fields for curvilinear coordinates

        def _cartesian_x(field, data):
            if registry.ds.geometry == "polar":

                return data[(ftype, f"{basename}_r")] * np.cos(data[(ftype, "theta")])

            elif registry.ds.geometry == "cylindrical":

                if data.ds.dimensionality == 2:
                    return data[(ftype, f"{basename}_r")]
                elif data.ds.dimensionality == 3:
                    return data[(ftype, f"{basename}_r")] * np.cos(
                        data[(ftype, "theta")]
                    ) - data[(ftype, f"{basename}_theta")] * np.sin(
                        data[(ftype, "theta")]
                    )

            elif registry.ds.geometry == "spherical":

                if data.ds.dimensionality == 2:
                    return data[(ftype, f"{basename}_r")] * np.sin(
                        data[(ftype, "theta")]
                    ) + data[(ftype, f"{basename}_theta")] * np.cos(
                        data[(ftype, "theta")]
                    )
                elif data.ds.dimensionality == 3:
                    return (
                        data[(ftype, f"{basename}_r")]
                        * np.sin(data[(ftype, "theta")])
                        * np.cos(data[(ftype, "phi")])
                        + data[(ftype, f"{basename}_theta")]
                        * np.cos(data[(ftype, "theta")])
                        * np.cos([(ftype, "phi")])
                        - data[(ftype, f"{basename}_phi")]
                        * np.sin(data[(ftype, "phi")])
                    )

        # it's redundant to define a cartesian x field for 1D data
        if registry.ds.dimensionality > 1:
            registry.add_field(
                (ftype, f"{basename}_cartesian_x"),
                sampling_type="local",
                function=_cartesian_x,
                units=field_units,
                display_field=True,
            )

        def _cartesian_y(field, data):
            if registry.ds.geometry == "polar":

                return data[(ftype, f"{basename}_r")] * np.sin(data[(ftype, "theta")])

            elif registry.ds.geometry == "cylindrical":

                if data.ds.dimensionality == 2:
                    return data[(ftype, f"{basename}_z")]
                elif data.ds.dimensionality == 3:
                    return data[(ftype, f"{basename}_r")] * np.sin(
                        data[(ftype, "theta")]
                    ) + data[(ftype, f"{basename}_theta")] * np.cos(
                        data[(ftype, "theta")]
                    )

            elif registry.ds.geometry == "spherical":

                if data.ds.dimensionality == 2:
                    return data[(ftype, f"{basename}_r")] * np.cos(
                        data[(ftype, "theta")]
                    ) - data[f"{basename}_theta"] * np.sin(data[(ftype, "theta")])
                elif data.ds.dimensionality == 3:
                    return (
                        data[(ftype, f"{basename}_r")]
                        * np.sin(data[(ftype, "theta")])
                        * np.sin(data[(ftype, "phi")])
                        + data[(ftype, f"{basename}_theta")]
                        * np.cos(data[(ftype, "theta")])
                        * np.sin([(ftype, "phi")])
                        + data[(ftype, f"{basename}_phi")]
                        * np.cos(data[(ftype, "phi")])
                    )

        if registry.ds.dimensionality >= 2:
            registry.add_field(
                (ftype, f"{basename}_cartesian_y"),
                sampling_type="local",
                function=_cartesian_y,
                units=field_units,
                display_field=True,
            )

        def _cartesian_z(field, data):
            if registry.ds.geometry == "cylindrical":
                return data[(ftype, f"{basename}_z")]
            elif registry.ds.geometry == "spherical":
                return data[(ftype, f"{basename}_r")] * np.cos(
                    data[(ftype, "theta")]
                ) - data[(ftype, f"{basename}_theta")] * np.sin(data[(ftype, "theta")])

        if registry.ds.dimensionality >= 2:
            registry.add_field(
                (ftype, f"{basename}_cartesian_z"),
                sampling_type="local",
                function=_cartesian_z,
                units=field_units,
                display_field=True,
            )

    if registry.ds.geometry == "spherical":

        def _cylindrical_radius_component(field, data):
            return (
                np.sin(data[(ftype, "theta")]) * data[(ftype, f"{basename}_r")]
                + np.cos(data[(ftype, "theta")]) * data[(ftype, f"{basename}_theta")]
            )

        registry.add_field(
            (ftype, f"{basename}_cylindrical_radius"),
            sampling_type="local",
            function=_cylindrical_radius_component,
            units=field_units,
            display_field=True,
        )

        registry.alias(
            (ftype, f"{basename}_cylindrical_z"),
            (ftype, f"{basename}_cartesian_z"),
        )

        # define vector components appropriate for 'theta'-normal plots.
        # The projection plane is called 'conic plane' in the code base as well as docs.
        # Contrary to 'poloidal' and 'toroidal', this isn't a widely spread
        # naming convention, but here it is exposed to users as part of dedicated
        # field names, so it needs to be stable.
        def _conic_x(field, data):
            rax = axis_names.index("r")
            pax = axis_names.index("phi")
            bc = data.get_field_parameter(f"bulk_{basename}")
            return np.cos(data[ftype, "phi"]) * (
                data[ftype, f"{basename}_r"] - bc[rax]
            ) - np.sin(data[ftype, "phi"]) * (data[ftype, f"{basename}_phi"] - bc[pax])

        def _conic_y(field, data):
            rax = axis_names.index("r")
            pax = axis_names.index("phi")
            bc = data.get_field_parameter(f"bulk_{basename}")
            return np.sin(data[(ftype, "phi")]) * (
                data[(ftype, f"{basename}_r")] - bc[rax]
            ) + np.cos(data[(ftype, "phi")]) * (
                data[(ftype, f"{basename}_phi")] - bc[pax]
            )

        if registry.ds.dimensionality == 3:
            registry.add_field(
                (ftype, f"{basename}_conic_x"),
                sampling_type="local",
                function=_conic_x,
                units=field_units,
                display_field=True,
            )
            registry.add_field(
                (ftype, f"{basename}_conic_y"),
                sampling_type="local",
                function=_conic_y,
                units=field_units,
                display_field=True,
            )


def create_averaged_field(
    registry,
    basename,
    field_units,
    ftype="gas",
    slice_info=None,
    validators=None,
    weight="mass",
):

    if validators is None:
        validators = []
    validators += [ValidateSpatial(1, [(ftype, basename)])]

    def _averaged_field(field, data):
        def atleast_4d(array):
            if array.ndim == 3:
                return array[..., None]
            else:
                return array

        nx, ny, nz, ngrids = atleast_4d(data[(ftype, basename)]).shape
        new_field = data.ds.arr(
            np.zeros((nx - 2, ny - 2, nz - 2, ngrids), dtype=np.float64),
            (just_one(data[(ftype, basename)]) * just_one(data[(ftype, weight)])).units,
        )
        weight_field = data.ds.arr(
            np.zeros((nx - 2, ny - 2, nz - 2, ngrids), dtype=np.float64),
            data[(ftype, weight)].units,
        )
        i_i, j_i, k_i = np.mgrid[0:3, 0:3, 0:3]

        for i, j, k in zip(i_i.ravel(), j_i.ravel(), k_i.ravel()):
            sl = (
                slice(i, nx - (2 - i)),
                slice(j, ny - (2 - j)),
                slice(k, nz - (2 - k)),
            )
            new_field += (
                atleast_4d(data[(ftype, basename)])[sl]
                * atleast_4d(data[(ftype, weight)])[sl]
            )
            weight_field += atleast_4d(data[(ftype, weight)])[sl]

        # Now some fancy footwork
        new_field2 = data.ds.arr(
            np.zeros((nx, ny, nz, ngrids)), data[(ftype, basename)].units
        )
        new_field2[1:-1, 1:-1, 1:-1] = new_field / weight_field

        if data[(ftype, basename)].ndim == 3:
            return new_field2[..., 0]
        else:
            return new_field2

    registry.add_field(
        (ftype, f"averaged_{basename}"),
        sampling_type="cell",
        function=_averaged_field,
        units=field_units,
        validators=validators,
    )
