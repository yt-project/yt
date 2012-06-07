"""
Clump finding helper classes

Author: Britton Smith <Britton.Smith@colorado.edu>
Affiliation: University of Colorado at Boulder
License:
  Copyright (C) 2008-2011 Britton Smith.  All Rights Reserved.

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

import numpy as na
import copy

from yt.funcs import *

from .contour_finder import identify_contours

class Clump(object):
    children = None
    def __init__(self, data, parent, field, cached_fields = None, 
                 function=None, clump_info=None):
        self.parent = parent
        self.data = data
        self.quantities = data.quantities
        self.field = field
        self.min_val = self.data[field].min()
        self.max_val = self.data[field].max()
        self.cached_fields = cached_fields

        # List containing characteristics about clumps that are to be written 
        # out by the write routines.
        if clump_info is None:
            self.set_default_clump_info()
        else:
            # Clump info will act the same if add_info_item is called before or after clump finding.
            self.clump_info = copy.deepcopy(clump_info)

        # Function determining whether a clump is valid and should be kept.
        self.default_function = 'self.data.quantities["IsBound"](truncate=True,include_thermal_energy=True) > 1.0'
        if function is None:
            self.function = self.default_function
        else:
            self.function = function

        # Return value of validity function, saved so it does not have to be calculated again.
        self.function_value = None

    def add_info_item(self,quantity,format):
        "Adds an entry to clump_info list and tells children to do the same."

        self.clump_info.append({'quantity':quantity, 'format':format})
        if self.children is None: return
        for child in self.children:
            child.add_info_item(quantity,format)

    def set_default_clump_info(self):
        "Defines default entries in the clump_info array."

        # add_info_item is recursive so this function does not need to be.
        self.clump_info = []

        # Number of cells.
        self.add_info_item('self.data["CellMassMsun"].size','"Cells: %d" % value')
        # Gas mass in solar masses.
        self.add_info_item('self.data["CellMassMsun"].sum()','"Mass: %e Msolar" % value')
        # Volume-weighted Jeans mass.
        self.add_info_item('self.data.quantities["WeightedAverageQuantity"]("JeansMassMsun","CellVolume")',
                           '"Jeans Mass (vol-weighted): %.6e Msolar" % value')
        # Mass-weighted Jeans mass.
        self.add_info_item('self.data.quantities["WeightedAverageQuantity"]("JeansMassMsun","CellMassMsun")',
                           '"Jeans Mass (mass-weighted): %.6e Msolar" % value')
        # Max level.
        self.add_info_item('self.data["GridLevel"].max()','"Max grid level: %d" % value')
        # Minimum number density.
        self.add_info_item('self.data["NumberDensity"].min()','"Min number density: %.6e cm^-3" % value')
        # Maximum number density.
        self.add_info_item('self.data["NumberDensity"].max()','"Max number density: %.6e cm^-3" % value')

    def clear_clump_info(self):
        "Clears the clump_info array and passes the instruction to its children."

        self.clump_info = []
        if self.children is None: return
        for child in self.children:
            child.clear_clump_info()

    def write_info(self,level,f_ptr):
        "Writes information for clump using the list of items in clump_info."

        for item in self.clump_info:
            # Call if callable, otherwise do an eval.
            if callable(item['quantity']):
                value = item['quantity']()
            else:
                value = eval(item['quantity'])
            output = eval(item['format'])
            f_ptr.write("%s%s" % ('\t'*level,output))
            f_ptr.write("\n")

    def find_children(self, min_val, max_val = None):
        if self.children is not None:
            print "Wiping out existing children clumps."
        self.children = []
        if max_val is None: max_val = self.max_val
        contour_info = identify_contours(self.data, self.field, min_val, max_val,
                                         self.cached_fields)
        for cid in contour_info:
            new_clump = self.data.extract_region(contour_info[cid])
            self.children.append(Clump(new_clump, self, self.field,
                                       self.cached_fields,function=self.function,
                                       clump_info=self.clump_info))

    def pass_down(self,operation):
        "Performs an operation on a clump with an exec and passes the instruction down to clump children."

        # Call if callable, otherwise do an exec.
        if callable(operation):
            operation()
        else:
            exec(operation)

        for child in self.children:
            child.pass_down(operation)

    def _isValid(self):
        "Perform user specified function to determine if child clumps should be kept."

        # Only call function if it has not been already.
        if self.function_value is None:
            self.function_value = eval(self.function)

        return self.function_value

    def __reduce__(self):
        return (_reconstruct_clump, 
                (self.parent, self.field, self.min_val, self.max_val,
                 self.function_value, self.children, self.data, self.clump_info, self.function))

    def __getitem__(self,request):
        return self.data[request]

def _reconstruct_clump(parent, field, mi, ma, function_value, children, data, clump_info, 
        function=None):
    obj = object.__new__(Clump)
    if iterable(parent):
        try:
            parent = parent[1]
        except KeyError:
            parent = parent
    if children is None: children = []
    obj.parent, obj.field, obj.min_val, obj.max_val, obj.function_value, obj.children, obj.clump_info, obj.function = \
        parent, field, mi, ma, function_value, children, clump_info, function
    # Now we override, because the parent/child relationship seems a bit
    # unreliable in the unpickling
    for child in children: child.parent = obj
    obj.data = data[1] # Strip out the PF
    obj.quantities = obj.data.quantities
    if obj.parent is None: return (data[0], obj)
    return obj

def find_clumps(clump, min_val, max_val, d_clump):
    print "Finding clumps: min: %e, max: %e, step: %f" % (min_val, max_val, d_clump)
    if min_val >= max_val: return
    clump.find_children(min_val)

    if (len(clump.children) == 1):
        find_clumps(clump, min_val*d_clump, max_val, d_clump)

    elif (len(clump.children) > 0):
        these_children = []
        print "Investigating %d children." % len(clump.children)
        for child in clump.children:
            find_clumps(child, min_val*d_clump, max_val, d_clump)
            if ((child.children is not None) and (len(child.children) > 0)):
                these_children.append(child)
            elif (child._isValid()):
                these_children.append(child)
            else:
                print "Eliminating invalid, childless clump with %d cells." % len(child.data["Ones"])
        if (len(these_children) > 1):
            print "%d of %d children survived." % (len(these_children),len(clump.children))            
            clump.children = these_children
        elif (len(these_children) == 1):
            print "%d of %d children survived, linking its children to parent." % (len(these_children),len(clump.children))
            clump.children = these_children[0].children
        else:
            print "%d of %d children survived, erasing children." % (len(these_children),len(clump.children))
            clump.children = []

def get_lowest_clumps(clump, clump_list=None):
    "Return a list of all clumps at the bottom of the hierarchy."

    if clump_list is None: clump_list = []
    if clump.children is None or len(clump.children) == 0:
        clump_list.append(clump)
    if clump.children is not None and len(clump.children) > 0:
        for child in clump.children:
            get_lowest_clumps(child, clump_list=clump_list)

    return clump_list

def write_clump_hierarchy(clump,level,f_ptr):
    for q in range(level):
        f_ptr.write("\t")
    f_ptr.write("Clump at level %d:\n" % level)
    clump.write_info(level,f_ptr)
    f_ptr.write("\n")
    f_ptr.flush()
    if ((clump.children is not None) and (len(clump.children) > 0)):
        for child in clump.children:
            write_clump_hierarchy(child,(level+1),f_ptr)

def write_clumps(clump,level,f_ptr):
    if ((clump.children is None) or (len(clump.children) == 0)):
        f_ptr.write("%sClump:\n" % ("\t"*level))
        clump.write_info(level,f_ptr)
        f_ptr.write("\n")
        f_ptr.flush()
    if ((clump.children is not None) and (len(clump.children) > 0)):
        for child in clump.children:
            write_clumps(child,0,f_ptr)

# Old clump info writing routines.
def write_old_clump_hierarchy(clump,level,f_ptr):
    for q in range(level):
        f_ptr.write("\t")
    f_ptr.write("Clump at level %d:\n" % level)
    clump.write_info(level,f_ptr)
    write_old_clump_info(clump,level,f_ptr)
    f_ptr.write("\n")
    f_ptr.flush()
    if ((clump.children is not None) and (len(clump.children) > 0)):
        for child in clump.children:
            write_clump_hierarchy(child,(level+1),f_ptr)

def write_old_clumps(clump,level,f_ptr):
    if ((clump.children is None) or (len(clump.children) == 0)):
        f_ptr.write("%sClump:\n" % ("\t"*level))
        write_old_clump_info(clump,level,f_ptr)
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

def write_old_clump_info(clump,level,f_ptr):
    fmt_dict = {'tl':  "\t" * level}
    fmt_dict['num_cells'] = clump.data["CellMassMsun"].size,
    fmt_dict['total_mass'] = clump.data["CellMassMsun"].sum()
    fmt_dict['jeans_mass_vol'] = clump.data.quantities["WeightedAverageQuantity"]("JeansMassMsun","CellVolume")
    fmt_dict['jeans_mass_mass'] = clump.data.quantities["WeightedAverageQuantity"]("JeansMassMsun","CellMassMsun")
    fmt_dict['max_level'] =  clump.data["GridLevel"].max()
    fmt_dict['min_density'] =  clump.data["NumberDensity"].min()
    fmt_dict['max_density'] =  clump.data["NumberDensity"].max()
    f_ptr.write(__clump_info_template % fmt_dict)

# Recipes for various clump calculations.
recipes = {}

# Distance from clump center of mass to center of mass of top level object.
def _DistanceToMainClump(master,units='pc'):
    masterCOM = master.data.quantities['CenterOfMass']()
    pass_command = "self.masterCOM = [%.10f, %.10f, %.10f]" % (masterCOM[0],
                                                               masterCOM[1],
                                                               masterCOM[2])
    master.pass_down(pass_command)
    master.pass_down("self.com = self.data.quantities['CenterOfMass']()")

    quantity = "((self.com[0]-self.masterCOM[0])**2 + (self.com[1]-self.masterCOM[1])**2 + (self.com[2]-self.masterCOM[2])**2)**(0.5)*self.data.pf.units['%s']" % units
    format = "%s%s%s" % ("'Distance from center: %.6e ",units,"' % value")

    master.add_info_item(quantity,format)

recipes['DistanceToMainClump'] = _DistanceToMainClump
