from typing import Any, Callable, Optional, Tuple, Union

import numpy as np
from mypy_extensions import NoReturn
from unyt.array import unyt_array
from unyt.unit_object import Unit

from yt.fields.derived_field import DerivedField
from yt.fields.field_detector import FieldDetector
from yt.frontends.ramses.data_structures import RAMSESDomainSubset
from yt.frontends.ramses.fields import RAMSESFieldInfo
from yt.funcs import mylog
from yt.geometry.geometry_handler import is_curvilinear
from yt.utilities.chemical_formulas import default_mu
from yt.utilities.lib.misc_utilities import obtain_relative_velocity_vector

from .derived_field import ValidateParameter, ValidateSpatial
from .field_plugin_registry import register_field_plugin
from .vector_operations import (
    create_averaged_field,
    create_magnitude_field,
    create_vector_fields,
)


@register_field_plugin
def setup_fluid_fields(
    registry: RAMSESFieldInfo, ftype: str = "gas", slice_info: Optional[Any] = None
) -> None:
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

    if unit_system.name == "cgs":
        mag_units = "magnetic_field_cgs"
    else:
        mag_units = "magnetic_field_mks"

    create_vector_fields(
        registry, "velocity", unit_system["velocity"], ftype, slice_info
    )
    create_vector_fields(
        registry, "magnetic_field", unit_system[mag_units], ftype, slice_info
    )

    def _cell_mass(field: DerivedField, data: FieldDetector) -> unyt_array:
        return data[ftype, "density"] * data[ftype, "cell_volume"]

    registry.add_field(
        (ftype, "cell_mass"),
        sampling_type="cell",
        function=_cell_mass,
        units=unit_system["mass"],
    )
    registry.alias((ftype, "mass"), (ftype, "cell_mass"))

    def _sound_speed(field: DerivedField, data: FieldDetector) -> unyt_array:
        tr = data.ds.gamma * data[ftype, "pressure"] / data[ftype, "density"]
        return np.sqrt(tr)

    registry.add_field(
        (ftype, "sound_speed"),
        sampling_type="local",
        function=_sound_speed,
        units=unit_system["velocity"],
    )

    def _radial_mach_number(field: DerivedField, data: FieldDetector) -> unyt_array:
        """ Radial component of M{|v|/c_sound} """
        tr = data[ftype, "radial_velocity"] / data[ftype, "sound_speed"]
        return np.abs(tr)

    registry.add_field(
        (ftype, "radial_mach_number"),
        sampling_type="local",
        function=_radial_mach_number,
        units="",
    )

    def _kin_energy(field: DerivedField, data: FieldDetector) -> unyt_array:
        v = obtain_relative_velocity_vector(data)
        return 0.5 * data[ftype, "density"] * (v ** 2).sum(axis=0)

    registry.add_field(
        (ftype, "kinetic_energy"),
        sampling_type="local",
        function=_kin_energy,
        units=unit_system["pressure"],
        validators=[ValidateParameter("bulk_velocity")],
    )

    def _mach_number(field: DerivedField, data: FieldDetector) -> unyt_array:
        """ M{|v|/c_sound} """
        return data[ftype, "velocity_magnitude"] / data[ftype, "sound_speed"]

    registry.add_field(
        (ftype, "mach_number"), sampling_type="local", function=_mach_number, units=""
    )

    def _courant_time_step(field: DerivedField, data: FieldDetector) -> unyt_array:
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

    def _pressure(field: DerivedField, data: FieldDetector) -> NoReturn:
        """ M{(Gamma-1.0)*rho*E} """
        tr = (data.ds.gamma - 1.0) * (
            data[ftype, "density"] * data[ftype, "thermal_energy"]
        )
        return tr

    registry.add_field(
        (ftype, "pressure"),
        sampling_type="local",
        function=_pressure,
        units=unit_system["pressure"],
    )

    def _kT(field: DerivedField, data: FieldDetector) -> unyt_array:
        return (pc.kboltz * data[ftype, "temperature"]).in_units("keV")

    registry.add_field(
        (ftype, "kT"),
        sampling_type="local",
        function=_kT,
        units="keV",
        display_name="Temperature",
    )

    def _metallicity(field: DerivedField, data: FieldDetector) -> NoReturn:
        return data[ftype, "metal_density"] / data[ftype, "density"]

    registry.add_field(
        (ftype, "metallicity"),
        sampling_type="local",
        function=_metallicity,
        units="Zsun",
    )

    def _metal_mass(field: DerivedField, data: FieldDetector) -> unyt_array:
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
                data["gas", "%s_number_density" % data.ds.field_info.species_names[0]]
            )
            for species in data.ds.field_info.species_names:
                field_data += data["gas", "%s_number_density" % species]
            return field_data

    else:

        def _number_density(field: DerivedField, data: FieldDetector) -> unyt_array:
            mu = getattr(data.ds, "mu", default_mu)
            return data[ftype, "density"] / (pc.mh * mu)

    registry.add_field(
        (ftype, "number_density"),
        sampling_type="local",
        function=_number_density,
        units=unit_system["number_density"],
    )

    def _mean_molecular_weight(field: DerivedField, data: FieldDetector) -> unyt_array:
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


def setup_gradient_fields(
    registry: RAMSESFieldInfo,
    grad_field: Tuple[str, str],
    field_units: Unit,
    slice_info: Optional[Any] = None,
) -> None:

    geom = registry.ds.geometry
    if is_curvilinear(geom):
        mylog.warning(
            "In %s geometry, gradient fields may contain artifacts near cartesian axes.",
            geom,
        )

    assert isinstance(grad_field, tuple)
    ftype, fname = grad_field
    if slice_info is None:
        sl_left = slice(None, -2, None)
        sl_right = slice(2, None, None)
        div_fac = 2.0
    else:
        sl_left, sl_right, div_fac = slice_info
    slice_3d = (slice(1, -1), slice(1, -1), slice(1, -1))

    def grad_func(axi: int, ax: str) -> Callable:
        slice_3dl = slice_3d[:axi] + (sl_left,) + slice_3d[axi + 1 :]
        slice_3dr = slice_3d[:axi] + (sl_right,) + slice_3d[axi + 1 :]

        def func(
            field: DerivedField, data: Union[FieldDetector, RAMSESDomainSubset]
        ) -> unyt_array:
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
            (ftype, "%s_gradient_%s" % (fname, ax)),
            sampling_type="local",
            function=f,
            validators=[ValidateSpatial(1, [grad_field])],
            units=grad_units,
        )

    create_magnitude_field(
        registry,
        "%s_gradient" % fname,
        grad_units,
        ftype=ftype,
        validators=[ValidateSpatial(1, [grad_field])],
    )
