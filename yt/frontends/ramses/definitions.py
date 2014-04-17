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

def ramses_header(hvals):
    header = ( ('ncpu', 1, 'i'),
               ('ndim', 1, 'i'),
               ('nx', 3, 'i'),
               ('nlevelmax', 1, 'i'),
               ('ngridmax', 1, 'i'),
               ('nboundary', 1, 'i'),
               ('ngrid_current', 1, 'i'),
               ('boxlen', 1, 'd'),
               ('nout', 3, 'I')
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
                 ('mass_sph', 1, 'd') )
    yield next_set
    tree_header = ( ('headl', hvals['nlevelmax'] * hvals['ncpu'], 'i'),
                    ('taill', hvals['nlevelmax'] * hvals['ncpu'], 'i'),
                    ('numbl', hvals['nlevelmax'] * hvals['ncpu'], 'i'),
                  )
    yield tree_header

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
