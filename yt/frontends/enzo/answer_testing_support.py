"""
Answer Testing support for Enzo.

Author: Britton Smith <brittonsmith@gmail.com>
Affiliation: Michigan State University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2012 Britton Smith.  All Rights Reserved.

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

from yt.testing import *
from yt.config import ytcfg
from yt.mods import *

from yt.utilities.answer_testing.framework import \
     AnswerTestingTest, \
     can_run_pf, \
     FieldValuesTest, \
     GridHierarchyTest, \
     GridValuesTest, \
     ProjectionValuesTest, \
     ParentageRelationshipsTest, \
     temp_cwd

def requires_outputlog(path = ".", prefix = ""):
    def ffalse(func):
        return lambda: None
    def ftrue(func):
        @wraps(func)
        def fyielder(*args, **kwargs):
            with temp_cwd(path):
                for t in func(*args, **kwargs):
                    if isinstance(t, AnswerTestingTest):
                        t.prefix = prefix
                    yield t
        return fyielder
    if os.path.exists("OutputLog"):
        return ftrue
    with temp_cwd(path):
        if os.path.exists("OutputLog"):
            return ftrue
    return ffalse
     
def standard_small_simulation(pf_fn, fields):
    if not can_run_pf(pf_fn): return
    dso = [None]
    yield GridHierarchyTest(pf_fn)
    yield ParentageRelationshipsTest(pf_fn)
    for field in fields:
        yield GridValuesTest(pf_fn, field)
        if 'particle' in field: continue
        for axis in [0, 1, 2]:
            for ds in dso:
                for weight_field in [None, "Density"]:
                    yield ProjectionValuesTest(
                        pf_fn, axis, field, weight_field,
                        ds)
                yield FieldValuesTest(
                        pf_fn, field, ds)
                    
class ShockTubeTest(object):
    def __init__(self, data_file, solution_file, fields, 
                 left_edges, right_edges, rtol, atol):
        self.solution_file = solution_file
        self.data_file = data_file
        self.fields = fields
        self.left_edges = left_edges
        self.right_edges = right_edges
        self.rtol = rtol
        self.atol = atol

    def __call__(self):
        # Read in the pf
        pf = load(self.data_file)  
        exact = self.get_analytical_solution() 

        ad = pf.h.all_data()
        position = ad['x']
        for k in self.fields:
            field = ad[k]
            for xmin, xmax in zip(self.left_edges, self.right_edges):
                mask = (position >= xmin)*(position <= xmax)
                exact_field = np.interp(position[mask], exact['pos'], exact[k]) 
                # yield test vs analytical solution 
                yield assert_allclose, field[mask], exact_field, \
                    self.rtol, self.atol

    def get_analytical_solution(self):
        # Reads in from file 
        pos, dens, vel, pres, inte = \
                np.loadtxt(self.solution_file, unpack=True)
        exact = {}
        exact['pos'] = pos
        exact['Density'] = dens
        exact['x-velocity'] = vel
        exact['Pressure'] = pres
        exact['ThermalEnergy'] = inte
        return exact
