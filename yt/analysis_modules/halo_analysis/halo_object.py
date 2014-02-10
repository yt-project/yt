"""
Halo object.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

class Halo(object):
    particles = None
    def __init__(self, halo_catalog):
        self.halo_catalog = halo_catalog
        self.quantities = {}
