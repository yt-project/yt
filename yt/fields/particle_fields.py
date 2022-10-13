import numpy as np

from yt.fields.derived_field import ValidateParameter, ValidateSpatial
from yt.units._numpy_wrapper_functions import uconcatenate, ucross
from yt.utilities.lib.misc_utilities import (
    obtain_position_vector,
    obtain_relative_velocity_vector,
)
from yt.utilities.math_utils import (
    get_cyl_r,
    get_cyl_r_component,
    get_cyl_theta,
    get_cyl_theta_component,
    get_cyl_z,
    get_cyl_z_component,
    get_sph_phi,
    get_sph_phi_component,
    get_sph_r_component,
    get_sph_theta,
    get_sph_theta_component,
    modify_reference_frame,
)

from .field_functions import get_radius
from .vector_operations import create_magnitude_field

sph_whitelist_fields = (
    "density",
    "temperature",
    "metallicity",
    "thermal_energy",
    "smoothing_length",
    "H_fraction",
    "He_fraction",
    "C_fraction",
    "Ca_fraction",
    "N_fraction",
    "O_fraction",
    "S_fraction",
    "Ne_fraction",
    "Mg_fraction",
    "Si_fraction",
    "Fe_fraction",
    "Na_fraction",
    "Al_fraction",
    "Ar_fraction",
    "Ni_fraction",
    "H_density",
    "He_density",
    "C_density",
    "Ca_density",
    "N_density",
    "O_density",
    "S_density",
    "Ne_density",
    "Mg_density",
    "Si_density",
    "Fe_density",
    "Na_density",
    "Al_density",
    "Ar_density",
    "Ni_density",
)


def _field_concat(fname):
    def _AllFields(field, data):
        v = []
        for ptype in data.ds.particle_types:
            data.ds._last_freq = (ptype, None)
            if ptype == "all" or ptype in data.ds.known_filters:
                continue
            v.append(data[ptype, fname].copy())
        rv = uconcatenate(v, axis=0)
        return rv

    return _AllFields


def _field_concat_slice(fname, axi):
    def _AllFields(field, data):
        v = []
        for ptype in data.ds.particle_types:
            data.ds._last_freq = (ptype, None)
            if ptype == "all" or ptype in data.ds.known_filters:
                continue
            v.append(data[ptype, fname][:, axi])
        rv = uconcatenate(v, axis=0)
        return rv

    return _AllFields


def particle_deposition_functions(ptype, coord_name, mass_name, registry):
    unit_system = registry.ds.unit_system
    orig = set(registry.keys())
    ptype_dn = ptype.replace("_", " ").title()

    def particle_count(field, data):
        pos = data[ptype, coord_name]
        d = data.deposit(pos, method="count")
        return data.apply_units(d, field.units)

    registry.add_field(
        ("deposit", f"{ptype}_count"),
        sampling_type="cell",
        function=particle_count,
        validators=[ValidateSpatial()],
        units="",
        display_name=r"\mathrm{%s Count}" % ptype_dn,
    )

    def particle_mass(field, data):
        pos = data[ptype, coord_name]
        pmass = data[ptype, mass_name]
        pmass.convert_to_units(field.units)
        d = data.deposit(pos, [pmass], method="sum")
        return data.apply_units(d, field.units)

    registry.add_field(
        ("deposit", f"{ptype}_mass"),
        sampling_type="cell",
        function=particle_mass,
        validators=[ValidateSpatial()],
        display_name=r"\mathrm{%s Mass}" % ptype_dn,
        units=unit_system["mass"],
    )

    def particle_density(field, data):
        pos = data[ptype, coord_name]
        pos.convert_to_units("code_length")
        mass = data[ptype, mass_name]
        mass.convert_to_units("code_mass")
        d = data.deposit(pos, [mass], method="sum")
        d = data.ds.arr(d, "code_mass")
        d /= data["index", "cell_volume"]
        return d

    registry.add_field(
        ("deposit", f"{ptype}_density"),
        sampling_type="cell",
        function=particle_density,
        validators=[ValidateSpatial()],
        display_name=r"\mathrm{%s Density}" % ptype_dn,
        units=unit_system["density"],
    )

    def particle_cic(field, data):
        pos = data[ptype, coord_name]
        d = data.deposit(pos, [data[ptype, mass_name]], method="cic")
        d = data.apply_units(d, data[ptype, mass_name].units)
        d /= data["index", "cell_volume"]
        return d

    registry.add_field(
        ("deposit", f"{ptype}_cic"),
        sampling_type="cell",
        function=particle_cic,
        validators=[ValidateSpatial()],
        display_name=r"\mathrm{%s CIC Density}" % ptype_dn,
        units=unit_system["density"],
    )

    def _get_density_weighted_deposit_field(fname, units, method):
        def _deposit_field(field, data):
            """
            Create a grid field for particle quantities weighted by particle
            mass, using cloud-in-cell deposit.
            """
            pos = data[ptype, "particle_position"]
            # Get back into density
            pden = data[ptype, "particle_mass"]
            top = data.deposit(pos, [pden * data[(ptype, fname)]], method=method)
            bottom = data.deposit(pos, [pden], method=method)
            top[bottom == 0] = 0.0
            bnz = bottom.nonzero()
            top[bnz] /= bottom[bnz]
            d = data.ds.arr(top, units=units)
            return d

        return _deposit_field

    for ax in "xyz":
        for method, name in zip(("cic", "sum"), ("cic", "nn")):
            function = _get_density_weighted_deposit_field(
                f"particle_velocity_{ax}", "code_velocity", method
            )
            registry.add_field(
                ("deposit", ("%s_" + name + "_velocity_%s") % (ptype, ax)),
                sampling_type="cell",
                function=function,
                units=unit_system["velocity"],
                take_log=False,
                validators=[ValidateSpatial(0)],
            )

    for method, name in zip(("cic", "sum"), ("cic", "nn")):
        function = _get_density_weighted_deposit_field("age", "code_time", method)
        registry.add_field(
            ("deposit", ("%s_" + name + "_age") % (ptype)),
            sampling_type="cell",
            function=function,
            units=unit_system["time"],
            take_log=False,
            validators=[ValidateSpatial(0)],
        )

    # Now some translation functions.

    def particle_ones(field, data):
        v = np.ones(data[ptype, coord_name].shape[0], dtype="float64")
        return data.apply_units(v, field.units)

    registry.add_field(
        (ptype, "particle_ones"),
        sampling_type="particle",
        function=particle_ones,
        units="",
        display_name=r"Particle Count",
    )

    def particle_mesh_ids(field, data):
        pos = data[ptype, coord_name]
        ids = np.zeros(pos.shape[0], dtype="float64") - 1
        # This is float64 in name only.  It will be properly cast inside the
        # deposit operation.
        # _ids = ids.view("float64")
        data.deposit(pos, [ids], method="mesh_id")
        return data.apply_units(ids, "")

    registry.add_field(
        (ptype, "mesh_id"),
        sampling_type="particle",
        function=particle_mesh_ids,
        validators=[ValidateSpatial()],
        units="",
    )

    return list(set(registry.keys()).difference(orig))


def particle_scalar_functions(ptype, coord_name, vel_name, registry):

    # Now we have to set up the various velocity and coordinate things.  In the
    # future, we'll actually invert this and use the 3-component items
    # elsewhere, and stop using these.

    # Note that we pass in _ptype here so that it's defined inside the closure.

    def _get_coord_funcs(axi, _ptype):
        def _particle_velocity(field, data):
            return data[_ptype, vel_name][:, axi]

        def _particle_position(field, data):
            return data[_ptype, coord_name][:, axi]

        return _particle_velocity, _particle_position

    for axi, ax in enumerate("xyz"):
        v, p = _get_coord_funcs(axi, ptype)
        registry.add_field(
            (ptype, f"particle_velocity_{ax}"),
            sampling_type="particle",
            function=v,
            units="code_velocity",
        )
        registry.add_field(
            (ptype, f"particle_position_{ax}"),
            sampling_type="particle",
            function=p,
            units="code_length",
        )


def particle_vector_functions(ptype, coord_names, vel_names, registry):

    unit_system = registry.ds.unit_system

    # This will column_stack a set of scalars to create vector fields.

    def _get_vec_func(_ptype, names):
        def particle_vectors(field, data):
            v = [data[_ptype, name].in_units(field.units) for name in names]
            return data.ds.arr(np.column_stack(v), v[0].units)

        return particle_vectors

    registry.add_field(
        (ptype, "particle_position"),
        sampling_type="particle",
        function=_get_vec_func(ptype, coord_names),
        units="code_length",
    )

    registry.add_field(
        (ptype, "particle_velocity"),
        sampling_type="particle",
        function=_get_vec_func(ptype, vel_names),
        units=unit_system["velocity"],
    )


def get_angular_momentum_components(ptype, data, spos, svel):
    if data.has_field_parameter("normal"):
        normal = data.get_field_parameter("normal")
    else:
        normal = data.ds.arr(
            [0.0, 0.0, 1.0], "code_length"
        )  # default to simulation axis
    pos = data.ds.arr([data[ptype, spos % ax] for ax in "xyz"]).T
    vel = data.ds.arr([data[ptype, f"relative_{svel % ax}"] for ax in "xyz"]).T
    return pos, vel, normal


def standard_particle_fields(
    registry, ptype, spos="particle_position_%s", svel="particle_velocity_%s"
):
    unit_system = registry.ds.unit_system

    def _particle_velocity_magnitude(field, data):
        """M{|v|}"""
        return np.sqrt(
            data[ptype, f"relative_{svel % 'x'}"] ** 2
            + data[ptype, f"relative_{svel % 'y'}"] ** 2
            + data[ptype, f"relative_{svel % 'z'}"] ** 2
        )

    registry.add_field(
        (ptype, "particle_velocity_magnitude"),
        sampling_type="particle",
        function=_particle_velocity_magnitude,
        take_log=False,
        units=unit_system["velocity"],
    )

    def _particle_specific_angular_momentum(field, data):
        """Calculate the angular of a particle velocity.

        Returns a vector for each particle.
        """
        center = data.get_field_parameter("center")
        pos, vel, normal = get_angular_momentum_components(ptype, data, spos, svel)
        L, r_vec, v_vec = modify_reference_frame(center, normal, P=pos, V=vel)
        # adding in the unit registry allows us to have a reference to the
        # dataset and thus we will always get the correct units after applying
        # the cross product.
        return ucross(r_vec, v_vec, registry=data.ds.unit_registry)

    registry.add_field(
        (ptype, "particle_specific_angular_momentum"),
        sampling_type="particle",
        function=_particle_specific_angular_momentum,
        units=unit_system["specific_angular_momentum"],
        validators=[ValidateParameter("center")],
    )

    def _get_spec_ang_mom_comp(axi, ax, _ptype):
        def _particle_specific_angular_momentum_component(field, data):
            return data[_ptype, "particle_specific_angular_momentum"][:, axi]

        def _particle_angular_momentum_component(field, data):
            return (
                data[_ptype, "particle_mass"]
                * data[ptype, f"particle_specific_angular_momentum_{ax}"]
            )

        return (
            _particle_specific_angular_momentum_component,
            _particle_angular_momentum_component,
        )

    for axi, ax in enumerate("xyz"):
        f, v = _get_spec_ang_mom_comp(axi, ax, ptype)
        registry.add_field(
            (ptype, f"particle_specific_angular_momentum_{ax}"),
            sampling_type="particle",
            function=f,
            units=unit_system["specific_angular_momentum"],
            validators=[ValidateParameter("center")],
        )
        registry.add_field(
            (ptype, f"particle_angular_momentum_{ax}"),
            sampling_type="particle",
            function=v,
            units=unit_system["angular_momentum"],
            validators=[ValidateParameter("center")],
        )

    def _particle_angular_momentum(field, data):
        am = (
            data[ptype, "particle_mass"]
            * data[ptype, "particle_specific_angular_momentum"].T
        )
        return am.T

    registry.add_field(
        (ptype, "particle_angular_momentum"),
        sampling_type="particle",
        function=_particle_angular_momentum,
        units=unit_system["angular_momentum"],
        validators=[ValidateParameter("center")],
    )

    create_magnitude_field(
        registry,
        "particle_angular_momentum",
        unit_system["angular_momentum"],
        sampling_type="particle",
        ftype=ptype,
    )

    def _particle_radius(field, data):
        """The spherical radius component of the particle positions

        Relative to the coordinate system defined by the *normal* vector,
        and *center* field parameters.
        """
        return get_radius(data, "particle_position_", field.name[0])

    registry.add_field(
        (ptype, "particle_radius"),
        sampling_type="particle",
        function=_particle_radius,
        units=unit_system["length"],
        validators=[ValidateParameter("center")],
    )

    def _relative_particle_position(field, data):
        """The cartesian particle positions in a rotated reference frame

        Relative to the coordinate system defined by *center* field parameter.

        Note that the orientation of the x and y axes are arbitrary.
        """
        field_names = [(ptype, f"particle_position_{ax}") for ax in "xyz"]
        return obtain_position_vector(data, field_names=field_names).T

    registry.add_field(
        (ptype, "relative_particle_position"),
        sampling_type="particle",
        function=_relative_particle_position,
        units=unit_system["length"],
        validators=[ValidateParameter("normal"), ValidateParameter("center")],
    )

    def _relative_particle_velocity(field, data):
        """The vector particle velocities in an arbitrary coordinate system

        Relative to the coordinate system defined by the *bulk_velocity*
        vector field parameter.

        Note that the orientation of the x and y axes are arbitrary.
        """
        field_names = [(ptype, f"particle_velocity_{ax}") for ax in "xyz"]
        return obtain_relative_velocity_vector(data, field_names=field_names).T

    registry.add_field(
        (ptype, "relative_particle_velocity"),
        sampling_type="particle",
        function=_relative_particle_velocity,
        units=unit_system["velocity"],
        validators=[ValidateParameter("normal"), ValidateParameter("center")],
    )

    def _get_coord_funcs_relative(axi, _ptype):
        def _particle_pos_rel(field, data):
            return data[_ptype, "relative_particle_position"][:, axi]

        def _particle_vel_rel(field, data):
            return data[_ptype, "relative_particle_velocity"][:, axi]

        return _particle_vel_rel, _particle_pos_rel

    for axi, ax in enumerate("xyz"):
        v, p = _get_coord_funcs_relative(axi, ptype)
        registry.add_field(
            (ptype, f"particle_velocity_relative_{ax}"),
            sampling_type="particle",
            function=v,
            units="code_velocity",
        )
        registry.add_field(
            (ptype, f"particle_position_relative_{ax}"),
            sampling_type="particle",
            function=p,
            units="code_length",
        )
        registry.add_field(
            (ptype, f"relative_particle_velocity_{ax}"),
            sampling_type="particle",
            function=v,
            units="code_velocity",
        )
        registry.add_field(
            (ptype, f"relative_particle_position_{ax}"),
            sampling_type="particle",
            function=p,
            units="code_length",
        )

    # this is just particle radius but we add it with an alias for the sake of
    # consistent naming
    registry.add_field(
        (ptype, "particle_position_spherical_radius"),
        sampling_type="particle",
        function=_particle_radius,
        units=unit_system["length"],
        validators=[ValidateParameter("normal"), ValidateParameter("center")],
    )

    def _particle_position_spherical_theta(field, data):
        """The spherical theta coordinate of the particle positions.

        Relative to the coordinate system defined by the *normal* vector
        and *center* field parameters.
        """
        normal = data.get_field_parameter("normal")
        pos = data[(ptype, "relative_particle_position")].T
        return data.ds.arr(get_sph_theta(pos, normal), "")

    registry.add_field(
        (ptype, "particle_position_spherical_theta"),
        sampling_type="particle",
        function=_particle_position_spherical_theta,
        units="",
        validators=[ValidateParameter("center"), ValidateParameter("normal")],
    )

    def _particle_position_spherical_phi(field, data):
        """The spherical phi component of the particle positions

        Relative to the coordinate system defined by the *normal* vector
        and *center* field parameters.
        """
        normal = data.get_field_parameter("normal")
        pos = data[(ptype, "relative_particle_position")].T
        return data.ds.arr(get_sph_phi(pos, normal), "")

    registry.add_field(
        (ptype, "particle_position_spherical_phi"),
        sampling_type="particle",
        function=_particle_position_spherical_phi,
        units="",
        validators=[ValidateParameter("normal"), ValidateParameter("center")],
    )

    def _particle_velocity_spherical_radius(field, data):
        """The spherical radius component of the particle velocities in an
         arbitrary coordinate system

        Relative to the coordinate system defined by the *normal* vector,
        *bulk_velocity* vector and *center* field parameters.
        """
        normal = data.get_field_parameter("normal")
        pos = data[(ptype, "relative_particle_position")].T
        vel = data[(ptype, "relative_particle_velocity")].T
        theta = get_sph_theta(pos, normal)
        phi = get_sph_phi(pos, normal)
        sphr = get_sph_r_component(vel, theta, phi, normal)
        return sphr

    registry.add_field(
        (ptype, "particle_velocity_spherical_radius"),
        sampling_type="particle",
        function=_particle_velocity_spherical_radius,
        units=unit_system["velocity"],
        validators=[ValidateParameter("normal"), ValidateParameter("center")],
    )

    registry.alias(
        (ptype, "particle_radial_velocity"),
        (ptype, "particle_velocity_spherical_radius"),
    )

    def _particle_velocity_spherical_theta(field, data):
        """The spherical theta component of the particle velocities in an
         arbitrary coordinate system

        Relative to the coordinate system defined by the *normal* vector,
        *bulk_velocity* vector and *center* field parameters.
        """
        normal = data.get_field_parameter("normal")
        pos = data[(ptype, "relative_particle_position")].T
        vel = data[(ptype, "relative_particle_velocity")].T
        theta = get_sph_theta(pos, normal)
        phi = get_sph_phi(pos, normal)
        spht = get_sph_theta_component(vel, theta, phi, normal)
        return spht

    registry.add_field(
        (ptype, "particle_velocity_spherical_theta"),
        sampling_type="particle",
        function=_particle_velocity_spherical_theta,
        units=unit_system["velocity"],
        validators=[ValidateParameter("normal"), ValidateParameter("center")],
    )

    def _particle_velocity_spherical_phi(field, data):
        """The spherical phi component of the particle velocities

        Relative to the coordinate system defined by the *normal* vector,
        *bulk_velocity* vector and *center* field parameters.
        """
        normal = data.get_field_parameter("normal")
        pos = data[(ptype, "relative_particle_position")].T
        vel = data[(ptype, "relative_particle_velocity")].T
        phi = get_sph_phi(pos, normal)
        sphp = get_sph_phi_component(vel, phi, normal)
        return sphp

    registry.add_field(
        (ptype, "particle_velocity_spherical_phi"),
        sampling_type="particle",
        function=_particle_velocity_spherical_phi,
        units=unit_system["velocity"],
        validators=[ValidateParameter("normal"), ValidateParameter("center")],
    )

    def _particle_position_cylindrical_radius(field, data):
        """The cylindrical radius component of the particle positions

        Relative to the coordinate system defined by the *normal* vector
        and *center* field parameters.
        """
        normal = data.get_field_parameter("normal")
        pos = data[(ptype, "relative_particle_position")].T
        pos.convert_to_units("code_length")
        return data.ds.arr(get_cyl_r(pos, normal), "code_length")

    registry.add_field(
        (ptype, "particle_position_cylindrical_radius"),
        sampling_type="particle",
        function=_particle_position_cylindrical_radius,
        units=unit_system["length"],
        validators=[ValidateParameter("normal"), ValidateParameter("center")],
    )

    def _particle_position_cylindrical_theta(field, data):
        """The cylindrical theta component of the particle positions

        Relative to the coordinate system defined by the *normal* vector
        and *center* field parameters.
        """
        normal = data.get_field_parameter("normal")
        pos = data[(ptype, "relative_particle_position")].T
        return data.ds.arr(get_cyl_theta(pos, normal), "")

    registry.add_field(
        (ptype, "particle_position_cylindrical_theta"),
        sampling_type="particle",
        function=_particle_position_cylindrical_theta,
        units="",
        validators=[ValidateParameter("center"), ValidateParameter("normal")],
    )

    def _particle_position_cylindrical_z(field, data):
        """The cylindrical z component of the particle positions

        Relative to the coordinate system defined by the *normal* vector
        and *center* field parameters.
        """
        normal = data.get_field_parameter("normal")
        pos = data[(ptype, "relative_particle_position")].T
        pos.convert_to_units("code_length")
        return data.ds.arr(get_cyl_z(pos, normal), "code_length")

    registry.add_field(
        (ptype, "particle_position_cylindrical_z"),
        sampling_type="particle",
        function=_particle_position_cylindrical_z,
        units=unit_system["length"],
        validators=[ValidateParameter("normal"), ValidateParameter("center")],
    )

    def _particle_velocity_cylindrical_radius(field, data):
        """The cylindrical radius component of the particle velocities

        Relative to the coordinate system defined by the *normal* vector,
        *bulk_velocity* vector and *center* field parameters.
        """
        normal = data.get_field_parameter("normal")
        pos = data[(ptype, "relative_particle_position")].T
        vel = data[(ptype, "relative_particle_velocity")].T
        theta = get_cyl_theta(pos, normal)
        cylr = get_cyl_r_component(vel, theta, normal)
        return cylr

    registry.add_field(
        (ptype, "particle_velocity_cylindrical_radius"),
        sampling_type="particle",
        function=_particle_velocity_cylindrical_radius,
        units=unit_system["velocity"],
        validators=[ValidateParameter("normal"), ValidateParameter("center")],
    )

    def _particle_velocity_cylindrical_theta(field, data):
        """The cylindrical theta component of the particle velocities

        Relative to the coordinate system defined by the *normal* vector,
        *bulk_velocity* vector and *center* field parameters.
        """
        normal = data.get_field_parameter("normal")
        pos = data[(ptype, "relative_particle_position")].T
        vel = data[(ptype, "relative_particle_velocity")].T
        theta = get_cyl_theta(pos, normal)
        cylt = get_cyl_theta_component(vel, theta, normal)
        return cylt

    registry.add_field(
        (ptype, "particle_velocity_cylindrical_theta"),
        sampling_type="particle",
        function=_particle_velocity_cylindrical_theta,
        units=unit_system["velocity"],
        validators=[ValidateParameter("normal"), ValidateParameter("center")],
    )

    def _particle_velocity_cylindrical_z(field, data):
        """The cylindrical z component of the particle velocities

        Relative to the coordinate system defined by the *normal* vector,
        *bulk_velocity* vector and *center* field parameters.
        """
        normal = data.get_field_parameter("normal")
        vel = data[(ptype, "relative_particle_velocity")].T
        cylz = get_cyl_z_component(vel, normal)
        return cylz

    registry.add_field(
        (ptype, "particle_velocity_cylindrical_z"),
        sampling_type="particle",
        function=_particle_velocity_cylindrical_z,
        units=unit_system["velocity"],
        validators=[ValidateParameter("normal"), ValidateParameter("center")],
    )


def add_particle_average(registry, ptype, field_name, weight=None, density=True):
    if weight is None:
        weight = (ptype, "particle_mass")
    field_units = registry[ptype, field_name].units

    def _pfunc_avg(field, data):
        pos = data[ptype, "particle_position"]
        f = data[ptype, field_name]
        wf = data[ptype, weight]
        f *= wf
        v = data.deposit(pos, [f], method="sum")
        w = data.deposit(pos, [wf], method="sum")
        v /= w
        if density:
            v /= data["index", "cell_volume"]
        v[np.isnan(v)] = 0.0
        return v

    fn = ("deposit", f"{ptype}_avg_{field_name}")
    registry.add_field(
        fn,
        sampling_type="cell",
        function=_pfunc_avg,
        validators=[ValidateSpatial(0)],
        units=field_units,
    )
    return fn


def add_nearest_neighbor_field(ptype, coord_name, registry, nneighbors=64):
    field_name = (ptype, f"nearest_neighbor_distance_{nneighbors}")

    def _nth_neighbor(field, data):
        pos = data[ptype, coord_name]
        pos.convert_to_units("code_length")
        distances = 0.0 * pos[:, 0]
        data.particle_operation(
            pos, [distances], method="nth_neighbor", nneighbors=nneighbors
        )
        # Now some quick unit conversions.
        return distances

    registry.add_field(
        field_name,
        sampling_type="particle",
        function=_nth_neighbor,
        validators=[ValidateSpatial(0)],
        units="code_length",
    )
    return [field_name]


def add_nearest_neighbor_value_field(ptype, coord_name, sampled_field, registry):
    """
    This adds a nearest-neighbor field, where values on the mesh are assigned
    based on the nearest particle value found.  This is useful, for instance,
    with voronoi-tesselations.
    """
    field_name = ("deposit", f"{ptype}_nearest_{sampled_field}")
    field_units = registry[ptype, sampled_field].units
    unit_system = registry.ds.unit_system

    def _nearest_value(field, data):
        pos = data[ptype, coord_name]
        pos = pos.convert_to_units("code_length")
        value = data[ptype, sampled_field].in_base(unit_system.name)
        rv = data.smooth(
            pos, [value], method="nearest", create_octree=True, nneighbors=1
        )
        rv = data.apply_units(rv, field_units)
        return rv

    registry.add_field(
        field_name,
        sampling_type="cell",
        function=_nearest_value,
        validators=[ValidateSpatial(0)],
        units=field_units,
    )
    return [field_name]


def add_union_field(registry, ptype, field_name, units):
    """
    Create a field that is the concatenation of multiple particle types.
    This allows us to create fields for particle unions using alias names.
    """

    def _cat_field(field, data):
        return uconcatenate(
            [data[dep_type, field_name] for dep_type in data.ds.particle_types_raw]
        )

    registry.add_field(
        (ptype, field_name), sampling_type="particle", function=_cat_field, units=units
    )
