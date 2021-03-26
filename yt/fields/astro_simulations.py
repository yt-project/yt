from .domain_context import DomainContext

# Here's how this all works:
#
#   1. We have a mapping here that defines fields we might expect to find in an
#      astrophysical simulation to units (not necessarily definitive) that we
#      may want to use them in.
#   2. Simulations and frontends will register aliases from fields (which can
#      utilize units) to the fields enumerated here.
#   3. This plugin can call additional plugins on the registry.


class AstroSimulation(DomainContext):

    # This is an immutable of immutables.  Note that we do not specify the
    # fluid type here, although in most cases we expect it to be "gas".
    _known_fluid_fields = (
        ("density", "g/cm**3"),
        ("number_density", "1/cm**3"),
        ("pressure", "dyne / cm**2"),
        ("specific_thermal_energy", "erg / g"),
        ("temperature", "K"),
        ("velocity_x", "cm / s"),
        ("velocity_y", "cm / s"),
        ("velocity_z", "cm / s"),
        ("magnetic_field_x", "gauss"),
        ("magnetic_field_y", "gauss"),
        ("magnetic_field_z", "gauss"),
        ("radiation_acceleration_x", "cm / s**2"),
        ("radiation_acceleration_y", "cm / s**2"),
        ("radiation_acceleration_z", "cm / s**2"),
    )

    # This set of fields can be applied to any particle type.
    _known_particle_fields = (
        ("particle_position_x", "cm"),
        ("particle_position_y", "cm"),
        ("particle_position_z", "cm"),
        ("particle_velocity_x", "cm / s"),
        ("particle_velocity_y", "cm / s"),
        ("particle_velocity_z", "cm / s"),
        ("particle_mass", "g"),
        ("particle_index", ""),
    )
