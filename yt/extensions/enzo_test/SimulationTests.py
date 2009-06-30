"""
This is a definition of a class for actually testing a result

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
from ProblemVerification import *
from DataProviders import *
test_registry = {}

class TestSimulation(object):
    class __metaclass__(type):
        def __init__(cls, name, b, d):
            type.__init__(cls, name, b, d)
            if not hasattr(cls,'_ttype_name'): return
            test_registry[cls._ttype_name] = cls

    def __init__(self, name, source, verifier):
        self.name = name
        self.source = source
        self.verifier = verifier

    def __call__(self, pf, plot=True):
        results = [self.verifier(d, plot)
            for d in self.source(pf)]
        return results

class HaloProfiles(TestSimulation):
    _ttype_name = "halo_profiles"
    def __init__(self, name, n_part=200):
        source = HaloSpheres(n_part=n_part)
        verifier = ProfileVerification(
            "Density", "Temperature", q1_nbins=32)
        TestSimulation.__init__(self, name, source, verifier)

class HaloBaryonMasses(TestSimulation):
    _ttype_name = "halo_baryon_masses"
    def __init__(self, name, n_part=200):
        source = HaloSpheres(n_part=n_part)
        verifier = TotalMassVerification()
        TestSimulation.__init__(self, name, source, verifier)

class RhoTempPDF(TestSimulation):
    _ttype_name = "rho_temp_pdf"
    def __init__(self, name):
        source = EntireDataset()
        verifier = ProfileVerification(
            "Density", "Temperature", q1_nbins=32)
        TestSimulation.__init__(self, name, source, verifier)

