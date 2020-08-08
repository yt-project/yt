"""
This module gathers all user-facing functions with a `load_` prefix.

"""
# note: in the future, functions could be moved here instead
# in which case, this file should be removed from flake8 ignore list in setup.cfg

# note: simulation() should be renamed load_simulation()
from .convenience import load, simulation
from .frontends.stream.api import (
    load_amr_grids,
    load_hexahedral_mesh,
    load_octree,
    load_particles,
    load_uniform_grid,
    load_unstructured_mesh,
)
from .utilities import load_sample
