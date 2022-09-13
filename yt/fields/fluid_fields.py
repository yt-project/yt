import numpy as np

from yt.units.dimensions import current_mks  # type: ignore
from yt.units.unit_object import Unit  # type: ignore
from yt.utilities.chemical_formulas import compute_mu
from yt.utilities.lib.misc_utilities import obtain_relative_velocity_vector

from .derived_field import ValidateParameter, ValidateSpatial
from .field_plugin_registry import register_field_plugin
from .vector_operations import (
    create_averaged_field,
    create_magnitude_field,
    create_vector_fields,
)


@register_field_plugin
def setup_fluid_fields(registry, ftype="gas", slice_info=None):
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

    unit_system = registry.ds.unit_system

    if unit_system.base_units[current_mks] is None:
        mag_units = "magnetic_field_cgs"
    else:
        mag_units = "magnetic_field_mks"

    create_vector_fields(
        registry, "velocity", unit_system["velocity"], ftype, slice_info
    )

    create_vector_fields(
        registry, "magnetic_field", unit_system[mag_units], ftype, slice_info
    )

    def _cell_mass(field, data):
        return data[ftype, "density"] * data[ftype, "cell_volume"]

    registry.add_field(
        (ftype, "cell_mass"),
        sampling_type="cell",
        function=_cell_mass,
        units=unit_system["mass"],
    )
    registry.alias((ftype, "mass"), (ftype, "cell_mass"))

    # momentum
    def momentum_xyz(v):
        def _momm(field, data):
            return data["gas", "mass"] * data["gas", f"velocity_{v}"]

        def _momd(field, data):
            return data["gas", "density"] * data["gas", f"velocity_{v}"]

        return _momm, _momd

    for v in registry.ds.coordinates.axis_order:
        _momm, _momd = momentum_xyz(v)
        registry.add_field(
            ("gas", f"momentum_{v}"),
            sampling_type="local",
            function=_momm,
            units=unit_system["momentum"],
        )
        registry.add_field(
            ("gas", f"momentum_density_{v}"),
            sampling_type="local",
            function=_momd,
            units=unit_system["density"] * unit_system["velocity"],
        )

    def _sound_speed(field, data):
        tr = data.ds.gamma * data[ftype, "pressure"] / data[ftype, "density"]
        return np.sqrt(tr)

    registry.add_field(
        (ftype, "sound_speed"),
        sampling_type="local",
        function=_sound_speed,
        units=unit_system["velocity"],
    )

    def _radial_mach_number(field, data):
        """Radial component of M{|v|/c_sound}"""
        tr = data[ftype, "radial_velocity"] / data[ftype, "sound_speed"]
        return np.abs(tr)

    registry.add_field(
        (ftype, "radial_mach_number"),
        sampling_type="local",
        function=_radial_mach_number,
        units="",
    )

    def _kinetic_energy_density(field, data):
        v = obtain_relative_velocity_vector(data)
        return 0.5 * data[ftype, "density"] * (v**2).sum(axis=0)

    registry.add_field(
        (ftype, "kinetic_energy_density"),
        sampling_type="local",
        function=_kinetic_energy_density,
        units=unit_system["pressure"],
        validators=[ValidateParameter("bulk_velocity")],
    )

    def _mach_number(field, data):
        """M{|v|/c_sound}"""
        return data[ftype, "velocity_magnitude"] / data[ftype, "sound_speed"]

    registry.add_field(
        (ftype, "mach_number"), sampling_type="local", function=_mach_number, units=""
    )

    def _courant_time_step(field, data):
        t1 = data[ftype, "dx"] / (
            data[ftype, "sound_speed"] + np.abs(data[ftype, "velocity_x"])
        )
        t2 = data[ftype, "dy"] / (
            data[ftype, "sound_speed"] + np.abs(data[ftype, "velocity_y"])
        )
        t3 = data[ftype, "dz"] / (
            data[ftype, "sound_speed"] + np.abs(data[ftype, "velocity_z"])
        )
        tr = np.minimum(np.minimum(t1, t2), t3)
        return tr

    registry.add_field(
        (ftype, "courant_time_step"),
        sampling_type="cell",
        function=_courant_time_step,
        units=unit_system["time"],
    )

    def _pressure(field, data):
        """M{(Gamma-1.0)*rho*E}"""
        tr = (data.ds.gamma - 1.0) * (
            data[ftype, "density"] * data[ftype, "specific_thermal_energy"]
        )
        return tr

    registry.add_field(
        (ftype, "pressure"),
        sampling_type="local",
        function=_pressure,
        units=unit_system["pressure"],
    )

    def _kT(field, data):
        return (pc.kboltz * data[ftype, "temperature"]).in_units("keV")

    registry.add_field(
        (ftype, "kT"),
        sampling_type="local",
        function=_kT,
        units="keV",
        display_name="Temperature",
    )

    def _metallicity(field, data):
        return data[ftype, "metal_density"] / data[ftype, "density"]

    registry.add_field(
        (ftype, "metallicity"),
        sampling_type="local",
        function=_metallicity,
        units="Zsun",
    )

    def _metal_mass(field, data):
        Z = data[ftype, "metallicity"].to("dimensionless")
        return Z * data[ftype, "mass"]

    registry.add_field(
        (ftype, "metal_mass"),
        sampling_type="local",
        function=_metal_mass,
        units=unit_system["mass"],
    )

    if len(registry.ds.field_info.species_names) > 0:

        def _number_density(field, data):
            field_data = np.zeros_like(
                data["gas", f"{data.ds.field_info.species_names[0]}_number_density"]
            )
            for species in data.ds.field_info.species_names:
                field_data += data["gas", f"{species}_number_density"]
            return field_data

    else:

        def _number_density(field, data):
            mu = getattr(data.ds, "mu", compute_mu(data.ds.default_species_fields))
            return data[ftype, "density"] / (pc.mh * mu)

    registry.add_field(
        (ftype, "number_density"),
        sampling_type="local",
        function=_number_density,
        units=unit_system["number_density"],
    )

    def _mean_molecular_weight(field, data):
        return data[ftype, "density"] / (pc.mh * data[ftype, "number_density"])

    registry.add_field(
        (ftype, "mean_molecular_weight"),
        sampling_type="local",
        function=_mean_molecular_weight,
        units="",
    )

    setup_gradient_fields(
        registry, (ftype, "pressure"), unit_system["pressure"], slice_info
    )

    setup_gradient_fields(
        registry, (ftype, "density"), unit_system["density"], slice_info
    )

    create_averaged_field(
        registry,
        "density",
        unit_system["density"],
        ftype=ftype,
        slice_info=slice_info,
        weight="mass",
    )


def setup_gradient_fields(registry, grad_field, field_units, slice_info=None):

    assert isinstance(grad_field, tuple)
    ftype, fname = grad_field
    if slice_info is None:
        sl_left = slice(None, -2, None)
        sl_right = slice(2, None, None)
        div_fac = 2.0
    else:
        sl_left, sl_right, div_fac = slice_info
    slice_3d = (slice(1, -1), slice(1, -1), slice(1, -1))

    def grad_func(axi, ax):
        slice_3dl = slice_3d[:axi] + (sl_left,) + slice_3d[axi + 1 :]
        slice_3dr = slice_3d[:axi] + (sl_right,) + slice_3d[axi + 1 :]

        def func(field, data):
            block_order = getattr(data, "_block_order", "C")
            if block_order == "F":
                # Fortran-ordering: we need to swap axes here and
                # reswap below
                field_data = data[grad_field].swapaxes(0, 2)
            else:
                field_data = data[grad_field]
            dx = div_fac * data[ftype, f"d{ax}"]
            if ax == "theta":
                dx *= data[ftype, "r"]
            if ax == "phi":
                dx *= data[ftype, "r"] * np.sin(data[ftype, "theta"])
            f = field_data[slice_3dr] / dx[slice_3d]
            f -= field_data[slice_3dl] / dx[slice_3d]
            new_field = np.zeros_like(data[grad_field], dtype=np.float64)
            new_field = data.ds.arr(new_field, field_data.units / dx.units)
            new_field[slice_3d] = f

            if block_order == "F":
                new_field = new_field.swapaxes(0, 2)

            return new_field

        return func

    field_units = Unit(field_units, registry=registry.ds.unit_registry)
    grad_units = field_units / registry.ds.unit_system["length"]

    for axi, ax in enumerate(registry.ds.coordinates.axis_order):
        f = grad_func(axi, ax)
        registry.add_field(
            (ftype, f"{fname}_gradient_{ax}"),
            sampling_type="local",
            function=f,
            validators=[ValidateSpatial(1, [grad_field])],
            units=grad_units,
        )

    create_magnitude_field(
        registry,
        f"{fname}_gradient",
        grad_units,
        ftype=ftype,
        validators=[ValidateSpatial(1, [grad_field])],
    )
