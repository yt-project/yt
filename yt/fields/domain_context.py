"""
Domain context base class

Currently we largely apply this to the fields that get loaded.  Presumably
different analysis operations could be the subject of this type of examination
as well.
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

domain_context_registry = {}

class DomainContext(object):
    class __metaclass__(type):
        def __init__(cls, name, b, d):
            type.__init__(cls, name, b, d)
            domain_context_registry[name] = cls

    _known_fluid_fields     =  ()
    _known_particle_fields  =  ()

    def __init__(self, ds):
        self.ds = ds

