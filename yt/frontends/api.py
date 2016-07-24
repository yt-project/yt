"""
API for yt.frontends



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import sys, types, os, glob, time, importlib
from yt.extern.six.moves import cPickle as pickle

_frontends = [
    'art',
    'artio',
    'athena',
    'boxlib',
    'chombo',
    'eagle',
    'enzo',
    'exodus_ii',
    'fits',
    'flash',
    'gadget',
    'gadget_fof',
    'gamer',
    'gdf',
    'gizmo',
    'halo_catalog',
    'http_stream',
    'moab',
    'owls',
    'owls_subfind',
    'ramses',
    'rockstar',
    'sdf',
    'stream',
    'tipsy',
    'ytdata',
]

class _frontend_container:
    def __init__(self):
        for frontend in _frontends:
            _mod = "yt.frontends.%s.api" % frontend
            setattr(self, frontend, importlib.import_module(_mod))
        setattr(self, 'api', importlib.import_module('yt.frontends.api'))
        setattr(self, '__name__', 'yt.frontends.api')
