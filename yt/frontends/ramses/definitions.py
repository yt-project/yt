"""
Definitions for RAMSES files




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

# These functions are RAMSES-specific
from yt.config import ytcfg
from yt.funcs import mylog
import re

def ramses_header(hvals):
    header = ( ('ncpu', 1, 'i'),
               ('ndim', 1, 'i'),
               ('nx', 3, 'i'),
               ('nlevelmax', 1, 'i'),
               ('ngridmax', 1, 'i'),
               ('nboundary', 1, 'i'),
               ('ngrid_current', 1, 'i'),
               ('boxlen', 1, 'd'),
               ('nout', 3, 'i')
              )
    yield header
    # TODO: REMOVE
    noutput, iout, ifout = hvals['nout']
    next_set = ( ('tout', noutput, 'd'),
                 ('aout', noutput, 'd'),
                 ('t', 1, 'd'),
                 ('dtold', hvals['nlevelmax'], 'd'),
                 ('dtnew', hvals['nlevelmax'], 'd'),
                 ('nstep',  2, 'i'),
                 ('stat', 3, 'd'),
                 ('cosm', 7, 'd'),
                 ('timing', 5, 'd'),
                 ('mass_sph', 1, 'd', True)
                 )
    yield next_set

field_aliases = {
    'standard_five':     ('Density',
                          'x-velocity',
                          'y-velocity',
                          'z-velocity',
                          'Pressure'),
    'standard_six':      ('Density',
                          'x-velocity',
                          'y-velocity',
                          'z-velocity',
                          'Pressure',
                          'Metallicity'),

}

## Regular expressions used to parse file descriptors
VERSION_RE = re.compile(r'# version: *(\d+)')
# This will match comma-separated strings, discarding whitespaces
# on the left hand side
VAR_DESC_RE = re.compile(r'\s*([^\s]+),\s*([^\s]+),\s*([^\s]+)')


## Configure family mapping
particle_families = {
    'DM': 1,
    'star': 2,
    'cloud': 3,
    'dust': 4,
    'star_tracer': -2,
    'cloud_tracer': -3,
    'dust_tracer': -4,
    'gas_tracer': 0
}

if ytcfg.has_section('ramses-families'):
    for key in particle_families.keys():
        val = ytcfg.getint('ramses-families', key, fallback=None)
        if val is not None:
            mylog.info('Changing family %s from %s to %s' % (key, particle_families[key], val))
            particle_families[key] = val
