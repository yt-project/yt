"""
Halo Finding methods

Author: Britton Smith <brittonsmith@gmail.com>
Affiliation: MSU
Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2013 Britton Smith, Matthew Turk.  All Rights Reserved.

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

from .operator_registry import \
    hf_methods
from .halo_object import \
    Halo

import numpy as np

class HaloFindingMethod(object):
    pass


class HOPFindingMethod(HaloFindingMethod):
    def __init__(self, threshold = 160):
        self.threshold = 160

    def __call__(self, halo_catalog):
        from yt.analysis_modules.halo_finding.hop.EnzoHop import \
            RunHOP
        ds = halo_catalog.data_source
        densities, tags = RunHOP(
            ds["particle_position_x"], ds["particle_position_y"],
            ds["particle_position_z"], self.threshold)
        tids = np.unique(tags)
        pi = np.zeros(tags.shape, dtype="bool")
        for tid in tids:
            # In-place equal op
            pi = np.equal(tags, tid, pi)
            halo = Halo(pi)
            yield halo

hf_methods["hop"] = HOPFindingMethod

