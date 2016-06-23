"""
Clump finding helper classes



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import copy
import numpy as np
import uuid

from yt.fields.derived_field import \
    ValidateSpatial
from yt.funcs import mylog, iterable
from yt.extern.six import string_types

from .clump_info_items import \
    clump_info_registry
from .clump_validators import \
    clump_validator_registry
from .contour_finder import \
    identify_contours

def add_contour_field(ds, contour_key):
    def _contours(field, data):
        fd = data.get_field_parameter("contour_slices_%s" % contour_key)
        vals = data["index", "ones"] * -1
        if fd is None or fd == 0.0:
            return vals
        for sl, v in fd.get(data.id, []):
            vals[sl] = v
        return vals

    ds.add_field(("index", "contours_%s" % contour_key),
                 function=_contours,
                 validators=[ValidateSpatial(0)],
                 take_log=False,
                 display_field=False,
                 units='')

class Clump(object):
    children = None
    def __init__(self, data, field, parent=None,
                 clump_info=None, validators=None):
        self.data = data
        self.field = field
        self.parent = parent
        self.quantities = data.quantities
        self.min_val = self.data[field].min()
        self.max_val = self.data[field].max()

        if parent is not None:
            self.data.parent = self.parent.data

        # List containing characteristics about clumps that are to be written 
        # out by the write routines.
        if clump_info is None:
            self.set_default_clump_info()
        else:
            # Clump info will act the same if add_info_item is called 
            # before or after clump finding.
            self.clump_info = copy.deepcopy(clump_info)

        if validators is None:
            validators = []
        self.validators = validators
        # Return value of validity function.
        self.valid = None

    def add_validator(self, validator, *args, **kwargs):
        """
        Add a validating function to determine whether the clump should 
        be kept.
        """
        callback = clump_validator_registry.find(validator, *args, **kwargs)
        self.validators.append(callback)
        if self.children is None: return
        for child in self.children:
            child.add_validator(validator)
        
    def add_info_item(self, info_item, *args, **kwargs):
        "Adds an entry to clump_info list and tells children to do the same."

        callback = clump_info_registry.find(info_item, *args, **kwargs)
        self.clump_info.append(callback)
        if self.children is None: return
        for child in self.children:
            child.add_info_item(info_item)

    def set_default_clump_info(self):
        "Defines default entries in the clump_info array."

        # add_info_item is recursive so this function does not need to be.
        self.clump_info = []

        self.add_info_item("total_cells")
        self.add_info_item("cell_mass")

        if any("jeans" in f for f in self.data.pf.field_list):
            self.add_info_item("mass_weighted_jeans_mass")
            self.add_info_item("volume_weighted_jeans_mass")

        self.add_info_item("max_grid_level")

        if any("number_density" in f for f in self.data.pf.field_list):
            self.add_info_item("min_number_density")
            self.add_info_item("max_number_density")

    def clear_clump_info(self):
        """
        Clears the clump_info array and passes the instruction to its 
        children.
        """

        self.clump_info = []
        if self.children is None: return
        for child in self.children:
            child.clear_clump_info()

    def write_info(self, level, f_ptr):
        "Writes information for clump using the list of items in clump_info."

        for item in self.clump_info:
            value = item(self)
            f_ptr.write("%s%s\n" % ('\t'*level, value))

    def find_children(self, min_val, max_val = None):
        if self.children is not None:
            mylog.info("Wiping out existing children clumps: %d.",
                       len(self.children))
        self.children = []
        if max_val is None: max_val = self.max_val
        nj, cids = identify_contours(self.data, self.field, min_val, max_val)
        # Here, cids is the set of slices and values, keyed by the
        # parent_grid_id, that defines the contours.  So we can figure out all
        # the unique values of the contours by examining the list here.
        unique_contours = set([])
        for sl_list in cids.values():
            for sl, ff in sl_list:
                unique_contours.update(np.unique(ff))
        contour_key = uuid.uuid4().hex
        base_object = getattr(self.data, 'base_object', self.data)
        add_contour_field(base_object.ds, contour_key)
        for cid in sorted(unique_contours):
            if cid == -1: continue
            new_clump = base_object.cut_region(
                    ["obj['contours_%s'] == %s" % (contour_key, cid)],
                    {('contour_slices_%s' % contour_key): cids})
            if new_clump["ones"].size == 0:
                # This is to skip possibly duplicate clumps.
                # Using "ones" here will speed things up.
                continue
            self.children.append(Clump(new_clump, self.field, parent=self,
                                       clump_info=self.clump_info,
                                       validators=self.validators))

    def pass_down(self,operation):
        """
        Performs an operation on a clump with an exec and passes the 
        instruction down to clump children.
        """

        # Call if callable, otherwise do an exec.
        if callable(operation):
            operation()
        else:
            exec(operation)

        if self.children is None: return
        for child in self.children:
            child.pass_down(operation)

    def _validate(self):
        "Apply all user specified validator functions."

        # Only call functions if not done already.
        if self.valid is not None:
            return self.valid

        self.valid = True
        for validator in self.validators:
            self.valid &= validator(self)
            if not self.valid:
                break

        return self.valid

    def __reduce__(self):
        return (_reconstruct_clump, 
                (self.parent, self.field, self.min_val, self.max_val,
                 self.valid, self.children, self.data, self.clump_info, 
                 self.function))

    def __getitem__(self,request):
        return self.data[request]

def _reconstruct_clump(parent, field, mi, ma, valid, children, 
                       data, clump_info, 
        function=None):
    obj = object.__new__(Clump)
    if iterable(parent):
        try:
            parent = parent[1]
        except KeyError:
            parent = parent
    if children is None: children = []
    obj.parent, obj.field, obj.min_val, obj.max_val, \
      obj.valid, obj.children, obj.clump_info, obj.function = \
        parent, field, mi, ma, valid, children, clump_info, function
    # Now we override, because the parent/child relationship seems a bit
    # unreliable in the unpickling
    for child in children: child.parent = obj
    obj.data = data[1] # Strip out the PF
    obj.quantities = obj.data.quantities
    if obj.parent is None: return (data[0], obj)
    return obj

def find_clumps(clump, min_val, max_val, d_clump):
    mylog.info("Finding clumps: min: %e, max: %e, step: %f" % 
               (min_val, max_val, d_clump))
    if min_val >= max_val: return
    clump.find_children(min_val)

    if (len(clump.children) == 1):
        find_clumps(clump, min_val*d_clump, max_val, d_clump)

    elif (len(clump.children) > 0):
        these_children = []
        mylog.info("Investigating %d children." % len(clump.children))
        for child in clump.children:
            find_clumps(child, min_val*d_clump, max_val, d_clump)
            if ((child.children is not None) and (len(child.children) > 0)):
                these_children.append(child)
            elif (child._validate()):
                these_children.append(child)
            else:
                mylog.info(("Eliminating invalid, childless clump with " +
                            "%d cells.") % len(child.data["ones"]))
        if (len(these_children) > 1):
            mylog.info("%d of %d children survived." %
                       (len(these_children),len(clump.children)))
            clump.children = these_children
        elif (len(these_children) == 1):
            mylog.info(("%d of %d children survived, linking its " +
                        "children to parent.") % 
                        (len(these_children),len(clump.children)))
            clump.children = these_children[0].children
        else:
            mylog.info("%d of %d children survived, erasing children." %
                       (len(these_children),len(clump.children)))
            clump.children = []

def get_lowest_clumps(clump, clump_list=None):
    "Return a list of all clumps at the bottom of the index."

    if clump_list is None: clump_list = []
    if clump.children is None or len(clump.children) == 0:
        clump_list.append(clump)
    if clump.children is not None and len(clump.children) > 0:
        for child in clump.children:
            get_lowest_clumps(child, clump_list=clump_list)

    return clump_list

def write_clump_index(clump, level, fh):
    top = False
    if isinstance(fh, string_types):
        fh = open(fh, "w")
        top = True
    for q in range(level):
        fh.write("\t")
    fh.write("Clump at level %d:\n" % level)
    clump.write_info(level, fh)
    fh.write("\n")
    fh.flush()
    if ((clump.children is not None) and (len(clump.children) > 0)):
        for child in clump.children:
            write_clump_index(child, (level+1), fh)
    if top:
        fh.close()

def write_clumps(clump, level, fh):
    top = False
    if isinstance(fh, string_types):
        fh = open(fh, "w")
        top = True
    if ((clump.children is None) or (len(clump.children) == 0)):
        fh.write("%sClump:\n" % ("\t"*level))
        clump.write_info(level, fh)
        fh.write("\n")
        fh.flush()
    if ((clump.children is not None) and (len(clump.children) > 0)):
        for child in clump.children:
            write_clumps(child, 0, fh)
    if top:
        fh.close()
