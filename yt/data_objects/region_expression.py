"""
An object that can live on the dataset to facilitate data access.

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import weakref

from yt.extern.six import string_types
from yt.utilities.exceptions import YTDimensionalityError

class RegionExpression(object):
    _all_data = None
    def __init__(self, ds):
        self.ds = weakref.proxy(ds)

    @property
    def all_data(self):
        if self._all_data is None:
            self._all_data = self.ds.all_data()
        return self._all_data

    def __getitem__(self, item):
        # At first, we will only implement this as accepting a slice that is
        # (optionally) unitful corresponding to a specific set of coordinates
        # that result in a rectangular prism or a slice.
        if isinstance(item, string_types):
            # This is some field; we will instead pass this back to the
            # all_data object.
            return self.all_data[item]
        if isinstance(item, tuple) and isinstance(item[1], string_types):
            return self.all_data[item]
        if isinstance(item, slice):
            item = (item, item, item)
        if len(item) != self.ds.dimensionality:
            # Not the right specification, and we don't want to do anything
            # implicitly.
            raise YTDimensionalityError(len(item), self.ds.dimensionality)
        if self.ds.dimensionality != 3:
            # We'll pass on this for the time being.
            raise RuntimeError

        # OK, now we need to look at our slices.  How many are a specific
        # coordinate?
        
        if not all(isinstance(v, slice) for v in item):
            return self._create_slice(item)
        else:
            if all(s.start is s.stop is s.step is None for s in item):
                return self.all_data
            return self._create_region(item)
            
    def _spec_to_value(self, input_tuple):
        if not isinstance(input_tuple, tuple):
            # We now assume that it's in code_length
            return self.ds.quan(input_tuple, 'code_length')
        v, u = input_tuple
        value = self.ds.quan(v, u)
        return value

    def _create_slice(self, slice_tuple):
        axis = None
        new_slice = []
        for ax, v in enumerate(slice_tuple):
            if not isinstance(v, slice):
                if axis is not None: raise RuntimeError
                axis = ax
                coord = self._spec_to_value(v)
                new_slice.append(slice(None, None, None))
            else:
                new_slice.append(v)
        # This new slice doesn't need to be a tuple
        source = self._create_region(new_slice)
        sl = self.ds.slice(axis, coord, data_source = source)
        return sl

    def _slice_to_edges(self, ax, val):
        if val.start is None:
            l = self.ds.domain_left_edge[ax]
        else:
            l = self._spec_to_value(val.start)
        if val.stop is None:
            r = self.ds.domain_right_edge[ax]
        else:
            r = self._spec_to_value(val.stop)
        if r < l:
            raise RuntimeError
        return l, r

    def _create_region(self, bounds_tuple):
        left_edge = []
        right_edge = []
        dims = []
        for ax, b in enumerate(bounds_tuple):
            l, r = self._slice_to_edges(ax, b)
            left_edge.append(l)
            right_edge.append(r)
            dims.append(getattr(b.step, "imag", None))
        center = [ (cl + cr)/2.0 for cl, cr in zip(left_edge, right_edge)]
        if all(d is not None for d in dims):
            return self.ds.arbitrary_grid(left_edge, right_edge, dims)
        return self.ds.region(center, left_edge, right_edge)
