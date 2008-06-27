"""
Clump finding helper classes

@author: Britton Smith
@organization: University of Colorado at Boulder
@contact: U{Britton.Smith@colorado.edu<mailto:Britton.Smith@colorado.edu>
@license:
  Copyright (C) 2008 Britton Smith.  All Rights Reserved.

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


from yt.lagos import *
import numpy as na

class Clump(object):
    children = None
    def __init__(self, data, parent, field):
        self.parent = parent
        self.data = data
        self.field = field
        self.min = self.data[field].min()
        self.max = self.data[field].max()
        self.isBound = None

    def find_children(self, min, max = None):
        if self.children is not None:
            print "Wiping out existing children clumps."
        self.children = []
        if max is None: max = self.max
        contour_info = identify_contours(self.data, self.field, min, max)
        for cid in contour_info:
            new_clump = self.data.extract_region(contour_info[cid])
            self.children.append(Clump(new_clump, self, self.field))

    def get_IsBound(self):
        if self.isBound is None:
            self.isBound = self.data.quantities["IsBound"](truncate=True,include_thermal_energy=True)
        return self.isBound

def find_clumps(clump, min, max, d_clump):
    print "Finding clumps: min: %e, max: %e, step: %f" % (min, max, d_clump)
    if min >= max: return
    clump.find_children(min)

    if (len(clump.children) == 1):
        find_clumps(clump, min*d_clump, max, d_clump)

    elif (len(clump.children) > 0):
        these_children = []
        print "Investigating %d children." % len(clump.children)
        for child in clump.children:
            find_clumps(child, min*d_clump, max, d_clump)
            if ((child.children is not None) and (len(child.children) > 0)):
                these_children.append(child)
            elif (child.get_IsBound() > 1.0):
                these_children.append(child)
            else:
                print "Eliminating unbound, childless clump with %d cells." % len(child.data["CellMassMsun"])
        if (len(these_children) > 1):
            print "%d of %d children survived." % (len(these_children),len(clump.children))            
            clump.children = these_children
        elif (len(these_children) == 1):
            print "%d of %d children survived, linking its children to parent." % (len(these_children),len(clump.children))
            clump.children = these_children[0].children
        else:
            print "%d of %d children survived, erasing children." % (len(these_children),len(clump.children))
            clump.children = []


def write_clump_hierarchy(clump,level,f_ptr):
    for q in range(level):
        f_ptr.write("\t")
    f_ptr.write("Clump at level %d:\n" % level)
    write_clump_info(clump,level,f_ptr)
    f_ptr.write("\n")
    f_ptr.flush()
    if ((clump.children is not None) and (len(clump.children) > 0)):
        for child in clump.children:
            write_clump_hierarchy(child,(level+1),f_ptr)

def write_clumps(clump,level,f_ptr):
    if ((clump.children is None) or (len(clump.children) == 0)):
        f_ptr.write("%sClump:\n" % ("\t"*level))
        write_clump_info(clump,level,f_ptr)
        f_ptr.write("\n")
        f_ptr.flush()
    if ((clump.children is not None) and (len(clump.children) > 0)):
        for child in clump.children:
            write_clumps(child,0,f_ptr)

__clump_info_template = \
"""
%(tl)sCells: %(num_cells)s
%(tl)sMass: %(total_mass).6e Msolar
%(tl)sJeans Mass (vol-weighted): %(jeans_mass_vol).6e Msolar
%(tl)sJeans Mass (mass-weighted): %(jeans_mass_mass).6e Msolar
%(tl)sMax grid level: %(max_level)s
%(tl)sMin number density: %(min_density).6e cm^-3
%(tl)sMax number density: %(max_density).6e cm^-3

"""

def write_clump_info(clump,level,f_ptr):
    fmt_dict = {'tl':  "\t" * level}
    fmt_dict['num_cells'] = clump.data["CellMassMsun"].size,
    fmt_dict['total_mass'] = clump.data["CellMassMsun"].sum()
    fmt_dict['jeans_mass_vol'] = clump.data.quantities["WeightedAverageQuantity"]("JeansMassMsun","CellVolume")
    fmt_dict['jeans_mass_mass'] = clump.data.quantities["WeightedAverageQuantity"]("JeansMassMsun","CellMassMsun")
    fmt_dict['max_level'] =  clump.data["GridLevel"].max()
    fmt_dict['min_density'] =  clump.data["NumberDensity"].min()
    fmt_dict['max_density'] =  clump.data["NumberDensity"].max()
    f_ptr.write(__clump_info_template % fmt_dict)
