"""
Dyablo frontend module for yt.

Dyablo is a regular-sized-block AMR hydrodynamical simulation code
with outputs in HDF5 format, ordered along a Morton curve.
"""

from .api import DyabloDataset

__all__ = ["DyabloDataset"]
