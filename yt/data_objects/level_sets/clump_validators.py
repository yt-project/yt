import numpy as np

from yt.utilities.lib.misc_utilities import gravitational_binding_energy
from yt.utilities.operator_registry import OperatorRegistry
from yt.utilities.physical_constants import gravitational_constant_cgs as G

clump_validator_registry = OperatorRegistry()


def add_validator(name, function):
    clump_validator_registry[name] = ClumpValidator(function)


class ClumpValidator:
    r"""
    A ClumpValidator is a function that takes a clump and returns
    True or False as to whether the clump is valid and shall be kept.
    """

    def __init__(self, function, args=None, kwargs=None):
        self.function = function
        self.args = args
        if self.args is None:
            self.args = []
        self.kwargs = kwargs
        if self.kwargs is None:
            self.kwargs = {}

    def __call__(self, clump):
        return self.function(clump, *self.args, **self.kwargs)


def _gravitationally_bound(
    clump, use_thermal_energy=True, use_particles=True, truncate=True, num_threads=0
):
    "True if clump is gravitationally bound."

    use_particles &= ("all", "particle_mass") in clump.data.ds.field_info

    bulk_velocity = clump.quantities.bulk_velocity(use_particles=use_particles)

    kinetic = (
        0.5
        * (
            clump["gas", "mass"]
            * (
                (bulk_velocity[0] - clump["gas", "velocity_x"]) ** 2
                + (bulk_velocity[1] - clump["gas", "velocity_y"]) ** 2
                + (bulk_velocity[2] - clump["gas", "velocity_z"]) ** 2
            )
        ).sum()
    )

    if use_thermal_energy:
        kinetic += (
            clump["gas", "mass"] * clump["gas", "specific_thermal_energy"]
        ).sum()

    if use_particles:
        kinetic += (
            0.5
            * (
                clump["all", "particle_mass"]
                * (
                    (bulk_velocity[0] - clump["all", "particle_velocity_x"]) ** 2
                    + (bulk_velocity[1] - clump["all", "particle_velocity_y"]) ** 2
                    + (bulk_velocity[2] - clump["all", "particle_velocity_z"]) ** 2
                )
            ).sum()
        )

    if use_particles:
        m = np.concatenate(
            [clump["gas", "mass"].in_cgs(), clump["all", "particle_mass"].in_cgs()]
        )
        px = np.concatenate(
            [clump["index", "x"].in_cgs(), clump["all", "particle_position_x"].in_cgs()]
        )
        py = np.concatenate(
            [clump["index", "y"].in_cgs(), clump["all", "particle_position_y"].in_cgs()]
        )
        pz = np.concatenate(
            [clump["index", "z"].in_cgs(), clump["all", "particle_position_z"].in_cgs()]
        )
    else:
        m = clump["gas", "mass"].in_cgs()
        px = clump["index", "x"].in_cgs()
        py = clump["index", "y"].in_cgs()
        pz = clump["index", "z"].in_cgs()

    potential = clump.data.ds.quan(
        G
        * gravitational_binding_energy(
            m, px, py, pz, truncate, (kinetic / G).in_cgs(), num_threads=num_threads
        ),
        kinetic.in_cgs().units,
    )

    if truncate and potential >= kinetic:
        return True

    return potential >= kinetic


add_validator("gravitationally_bound", _gravitationally_bound)


def _min_cells(clump, n_cells):
    "True if clump has a minimum number of cells."
    return clump["index", "ones"].size >= n_cells


add_validator("min_cells", _min_cells)
