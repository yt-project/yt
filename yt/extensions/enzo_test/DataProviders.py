"""
This is a definition of a class for providing data to a simulation
verification method

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

class TestDataSource(object):
    def provide(self, pf):
        pass

    def __call__(self, pf):
        return self.provide(pf)

class MaximumDensitySphere(TestDataSource):
    def __init__(self, radius, unit):
        self.radius = radius
        self.unit = unit

    def provide(self, pf):
        v, c = pf.h.find_max("Density")
        yield pf.h.sphere(c, self.radius/pf[self.unit])

class HaloSpheres(TestDataSource):
    def __init__(self, n_part):
        self.n_part = n_part

    def provide(self, pf):
       halo_list = HaloFinder(pf)
       for halo in halo_list:
            if halo.get_size() < self.n_part: continue
            yield halo.get_sphere()

class Halos(TestDataSource):
    def __init__(self):
        pass

    def provide(self, pf):
       halo_list = HaloFinder(pf)
       for halo in halo_list:
            yield halo

class EntireDataset(TestDataSource):
    def __init__(self):
        pass

    def provide(self, pf):
        yield pf.h.all_data()
