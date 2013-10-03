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

import sys, types, os, glob, cPickle, time, importlib

_frontends = [
    'enzo',
    'castro',
    'nyx',
    'orion',
    'flash',
    'tiger',
    'artio',
    'ramses',
    'chombo',
    'gdf',
    'moab',
    'athena',
    'art',
    'pluto',
    'stream',
    'sph',
]

class _frontend_container:
    def __init__(self):
        for frontend in _frontends:
            _mod = "yt.frontends.%s.api" % frontend
            setattr(self, frontend, importlib.import_module(_mod))
