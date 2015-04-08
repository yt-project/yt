"""
=============================================================
Spatial algorithms and data structures (:mod:`scipy.spatial`)
=============================================================

Nearest-neighbor queries:

.. autosummary::
   :toctree: generated/

   KDTree      -- class for efficient nearest-neighbor queries
   cKDTree     -- class for efficient nearest-neighbor queries (faster impl.)
   distance    -- module containing many different distance measures

Delaunay triangulation:

.. autosummary::
   :toctree: generated/

   Delaunay
   tsearch

"""
from __future__ import absolute_import

from .kdtree import *
from .ckdtree import *
#from qhull import *

__all__ = list(filter(lambda s: not s.startswith('_'), dir()))
__all__ += ['distance']

from . import distance
from numpy.testing import Tester
test = Tester().test
