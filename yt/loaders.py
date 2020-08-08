"""
This module gathers all user-facing functions with a `load_` prefix.

"""
# note: in the future, functions could be moved here instead
# in which case, this file should be removed from flake8 ignore list in setup.cfg

from .convenience import load, load_simulation
from .frontends.stream.api import (
    load_amr_grids,
    load_hexahedral_mesh,
    load_octree,
    load_particles,
    load_uniform_grid,
    load_unstructured_mesh,
)
from .utilities import load_sample
