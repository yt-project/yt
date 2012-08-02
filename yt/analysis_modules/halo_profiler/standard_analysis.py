"""
This module contains a near-replacement for enzo_anyl

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2011 Matthew Turk.  All Rights Reserved.

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

import os

from yt.data_objects.profiles import BinnedProfile1D
from yt.funcs import *

class StandardRadialAnalysis(object):
    def __init__(self, pf, center, radius, n_bins = 128, inner_radius = None):
        self.pf = pf
        # We actually don't want to replicate the handling of setting the
        # center here, so we will pass it to the sphere creator.
        # Note also that the sphere can handle (val, unit) for radius, so we
        # will grab that from the sphere as well
        self.obj = pf.h.sphere(center, radius)
        if inner_radius is None:
            inner_radius = pf.h.get_smallest_dx() * pf['cm']
        self.inner_radius = inner_radius
        self.outer_radius = self.obj.radius * pf['cm']
        self.n_bins = n_bins

    def setup_field_parameters(self):
        # First the bulk velocity
        bv = self.obj.quantities["BulkVelocity"]()
        self.obj.set_field_parameter("bulk_velocity", bv)

    def construct_profile(self):
        # inner_bound in cm, outer_bound in same
        # Note that in some cases, we will need to massage this object.
        prof = BinnedProfile1D(self.obj, self.n_bins, "Radius",
                               self.inner_radius, self.outer_radius)
        by_weights = defaultdict(list)
        for fspec in analysis_field_list:
            if isinstance(fspec, types.TupleType) and len(fspec) == 2:
                field, weight = fspec
            else:
                field, weight = fspec, "CellMassMsun"
            by_weights[weight].append(field)
        known_fields = set(self.pf.h.field_list + self.pf.h.derived_field_list)
        for weight, fields in by_weights.items():
            fields = set(fields)
            fields.intersection_update(known_fields)
            prof.add_fields(list(fields), weight=weight)
        self.prof = prof

    def plot_everything(self, dirname = None):
        if dirname is None: dirname = "%s_profile_plots/" % (self.pf)
        if not os.path.isdir(dirname): os.makedirs(dirname)
        import matplotlib; matplotlib.use("Agg")
        import pylab
        for field in self.prof.keys():
            if field in ("UsedBins", "Radius"): continue
            pylab.clf()
            pylab.loglog(self.prof["Radius"], self.prof[field], '-x')
            pylab.xlabel("Radius [cm]")
            pylab.ylabel("%s" % field)
            pylab.savefig(os.path.join(
                dirname, "Radius_%s.png" % (field.replace(" ","_"))))
