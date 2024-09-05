import numpy as np

from .derived_field import ValidateParameter
from .field_plugin_registry import register_field_plugin
from .vector_operations import create_magnitude_field


@register_field_plugin
def setup_astro_fields(registry, ftype="gas", slice_info=None):
    unit_system = registry.ds.unit_system
    pc = registry.ds.units.physical_constants
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

    def _dynamical_time(field, data):
        """
        sqrt(3 pi / (16 G rho))
        """
        return np.sqrt(3.0 * np.pi / (16.0 * pc.G * data[ftype, "density"]))

    registry.add_field(
        (ftype, "dynamical_time"),
        sampling_type="local",
        function=_dynamical_time,
        units=unit_system["time"],
    )

    def _jeans_mass(field, data):
        MJ_constant = (((5.0 * pc.kboltz) / (pc.G * pc.mh)) ** 1.5) * (
            3.0 / (4.0 * np.pi)
        ) ** 0.5
        u = (
            MJ_constant
            * (
                (data[ftype, "temperature"] / data[ftype, "mean_molecular_weight"])
                ** 1.5
            )
            * (data[ftype, "density"] ** (-0.5))
        )
        return u

    registry.add_field(
        (ftype, "jeans_mass"),
        sampling_type="local",
        function=_jeans_mass,
        units=unit_system["mass"],
    )

    def _emission_measure(field, data):
        dV = data[ftype, "mass"] / data[ftype, "density"]
        nenhdV = data[ftype, "H_nuclei_density"] * dV
        nenhdV *= data[ftype, "El_number_density"]
        return nenhdV

    registry.add_field(
        (ftype, "emission_measure"),
        sampling_type="local",
        function=_emission_measure,
        units=unit_system["number_density"],
    )

    def _mazzotta_weighting(field, data):
        # Spectroscopic-like weighting field for galaxy clusters
        # Only useful as a weight_field for temperature, metallicity, velocity
        ret = data[ftype, "El_number_density"].d ** 2
        ret *= data[ftype, "kT"].d ** -0.75
        return ret

    registry.add_field(
        (ftype, "mazzotta_weighting"),
        sampling_type="local",
        function=_mazzotta_weighting,
        units="",
    )

    def _optical_depth(field, data):
        return data[ftype, "El_number_density"] * pc.sigma_thompson

    registry.add_field(
        (ftype, "optical_depth"),
        sampling_type="local",
        function=_optical_depth,
        units=unit_system["length"] ** -1,
    )

    def _sz_kinetic(field, data):
        # minus sign is because radial velocity is WRT viewer
        # See issue #1225
        return -data[ftype, "velocity_los"] * data[ftype, "optical_depth"] / pc.clight

    registry.add_field(
        (ftype, "sz_kinetic"),
        sampling_type="local",
        function=_sz_kinetic,
        units=unit_system["length"] ** -1,
        validators=[ValidateParameter("axis", {"axis": [0, 1, 2]})],
    )

    def _szy(field, data):
        kT = data[ftype, "kT"] / (pc.me * pc.clight * pc.clight)
        return data[ftype, "optical_depth"] * kT

    registry.add_field(
        (ftype, "szy"),
        sampling_type="local",
        function=_szy,
        units=unit_system["length"] ** -1,
    )

    def _entropy(field, data):
        mgammam1 = -2.0 / 3.0
        tr = data[ftype, "kT"] * data[ftype, "El_number_density"] ** mgammam1
        return data.apply_units(tr, field.units)

    registry.add_field(
        (ftype, "entropy"), sampling_type="local", units="keV*cm**2", function=_entropy
    )

    def _lorentz_factor(field, data):
        b2 = data[ftype, "velocity_magnitude"].to_value("c")
        b2 *= b2
        return 1.0 / np.sqrt(1.0 - b2)

    registry.add_field(
        (ftype, "lorentz_factor"),
        sampling_type="local",
        units="",
        function=_lorentz_factor,
    )

    # 4-velocity spatial components
    def four_velocity_xyz(u):
        def _four_velocity(field, data):
            return data["gas", f"velocity_{u}"] * data["gas", "lorentz_factor"]

        return _four_velocity

    for u in registry.ds.coordinates.axis_order:
        registry.add_field(
            ("gas", f"four_velocity_{u}"),
            sampling_type="local",
            function=four_velocity_xyz(u),
            units=unit_system["velocity"],
        )

    # 4-velocity t-component
    def _four_velocity_t(field, data):
        return data["gas", "lorentz_factor"] * pc.clight

    registry.add_field(
        ("gas", "four_velocity_t"),
        sampling_type="local",
        function=_four_velocity_t,
        units=unit_system["velocity"],
    )

    create_magnitude_field(
        registry,
        "four_velocity",
        unit_system["velocity"],
        ftype=ftype,
    )
