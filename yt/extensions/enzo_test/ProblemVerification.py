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

class VerificationMechanism(object):
    def __init__(self):
        pass

    def __call__(self, pf):
        self.run(pf)

class ProfileVerification(VerificationMechanism):
    def __init__(self, q1, q2,
                 q1_limits = None, q1_nbins = 64):
        VerificationMechanism.__init__(self)
        self.q1 = q1
        self.q2 = q2
        self.q1_limits = q1_limits
        self.q1_nbins = q1_nbins

    def _setup_profile(self):
        pass

    def run(self, pf):
        data = pf.h.all_data()
        limits = self.q1_limits
        if limits is None:
            limits = data.quantities["Extrema"](
                self.q1, lazy_reader=True)[0]
        prof = BinnedProfile1D(
            data, self.q1_nbins, self.q1,
            limits[0], limits[1], lazy_reader=True)
        prof.add_fields(self.q2)
        return prof[self.q2].copy()

class RadialProfileVerification(VerificationMechanism):
    def __init__(self, radius, unit, q2,
                 r_limits = None, q2_limits = None,
                 r_nbins = 64, q2_nbins = 64):
        VerificationMechanism.__init__(self)
        self.radius = radius
        self.unit = unit
        self.q2 = q2
        self.q1_limits = q1_limits
        self.q2_limits = q2_limits
        self.q1_nbins = q1_nbins
        self.q2_nbins = q2_nbins

class TotalMassVerification(VerificationMechanism):
    def __init__(self):
        pass

    def run(self, pf):
        data = pf.h.all_data()
        return data.quantities["TotalQuantity"](
                "CellMassMsun", lazy_reader=True)
