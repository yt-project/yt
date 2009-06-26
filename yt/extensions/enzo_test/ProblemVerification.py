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
            callback_registry[cls._vtype_name] = cls

    def verify_identical(self, identifier = 'reference'):
        other = self._load_reference_result(identifier)
        return na.all(self.output == other.output)

    def _load_result(self, identifier, repo=None):
        pass

    def _store_result(self, identifier, repo=None):
        pass

class ProfileVerificiation(VerificationMechanism):
    _vtype_name = "profile"
    def __init__(self, name, data, q1, q2,
                 q1_limits = None, q2_limits = None,
                 q1_nbins = None, q2_nbins = None):
        self.name = name
        self.data = data
        self.q1 = q1
        self.q2 = q2
        self.q1_limits = q1_limits
        self.q2_limits = q2_limits
        self.q1_nbins = q1_nbins
        self.q2_nbins = q2_nbins

    def _setup_profile(self):
        pass

    def run(self):
        pass

class PhaseVerification(VerificationMechanism):
    _name = "phase"
