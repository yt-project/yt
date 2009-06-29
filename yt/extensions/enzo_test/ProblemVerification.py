"""
This is a definition of a set of classes for defining problem verification
methods, along with storage mechanisms.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2009 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from yt.mods import *
verification_registry = {}

class VerificationMechanism(object):
    class __metaclass__(type):
        def __init__(cls, name, b, d):
            type.__init__(cls, name, b, d)
            if not hasattr(cls,'_vtype_name'): return
            verification_registry[cls._vtype_name] = cls

    def __init__(self, name, pf):
        self.name = name
        self.pf = pf

class ProfileVerification(VerificationMechanism):
    _vtype_name = "profile"
    def __init__(self, name, pf, q1, q2,
                 q1_limits = None, q2_limits = None,
                 q1_nbins = 64, q2_nbins = 64):
        VerificationMechanism.__init__(self, name, pf)
        self.q1 = q1
        self.q2 = q2
        self.q1_limits = q1_limits
        self.q2_limits = q2_limits
        self.q1_nbins = q1_nbins
        self.q2_nbins = q2_nbins

    def _setup_profile(self):
        pass

    def run(self, pf):
        pass

class PhaseVerification(VerificationMechanism):
    _name = "phase"
