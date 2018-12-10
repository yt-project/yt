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
from yt.funcs import obj_length
from yt.units.yt_array import YTQuantity
from yt.utilities.exceptions import YTDimensionalityError
from yt.visualization.line_plot import LineBuffer

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
            if obj_length(item.start) == 3 and obj_length(item.stop) == 3:
                # This is for a ray that is not orthogonal to an axis.
                # it's straightforward to do this, so we create a ray
                # and drop out here.
                return self._create_ray(item)
            else:
                # This is for the case where we give a slice as an index; one
                # possible use case of this would be where we supply something
                # like ds.r[::256j] .  This would be expanded, implicitly into
                # ds.r[::256j, ::256j, ::256j].  Other cases would be if we do
                # ds.r[0.1:0.9] where it will be expanded along three dimensions.
                item = (item, item, item)
        if len(item) != self.ds.dimensionality:
            # Not the right specification, and we don't want to do anything
            # implicitly.  Note that this happens *after* the implicit expansion
            # of a single slice.
            raise YTDimensionalityError(len(item), self.ds.dimensionality)
        if self.ds.dimensionality != 3:
            # We'll pass on this for the time being.
            raise YTDimensionalityError(self.ds.dimensionality, '3')

        # OK, now we need to look at our slices.  How many are a specific
        # coordinate?

        nslices = sum(isinstance(v, slice) for v in item)
        if nslices == 0:
            return self._create_point(item)
        elif nslices == 1:
            return self._create_ortho_ray(item)
        elif nslices == 2:
            return self._create_slice(item)
        else:
            if all(s.start is s.stop is s.step is None for s in item):
                return self.all_data
            return self._create_region(item)

    def _spec_to_value(self, input):
        if isinstance(input, tuple):
            v = self.ds.quan(input[0], input[1]).to("code_length")
        elif isinstance(input, YTQuantity):
            v = self.ds.quan(input).to('code_length')
        else:
            v = self.ds.quan(input, "code_length")
        return v

    def _create_slice(self, slice_tuple):
        # This is somewhat more complex because we want to allow for slicing
        # in one dimension but also *not* using the entire domain; for instance
        # this means we allow something like ds.r[0.5, 0.1:0.4, 0.1:0.4].
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
        # Now, there's the possibility that what we're also seeing here
        # includes some steps, which would be for getting back a fixed
        # resolution buffer.  We check for that by checking if we have
        # exactly two imaginary steps.
        xax = self.ds.coordinates.x_axis[axis]
        yax = self.ds.coordinates.y_axis[axis]
        if getattr(new_slice[xax].step, "imag", 0.0) != 0.0 and \
           getattr(new_slice[yax].step, "imag", 0.0) != 0.0:
            # We now need to convert to a fixed res buffer.
            # We'll do this by getting the x/y axes, and then using that.
            width = source.right_edge[xax] - source.left_edge[xax]
            height = source.right_edge[yax] - source.left_edge[yax]
            # Make a resolution tuple with
            resolution = (int(new_slice[xax].step.imag),
                          int(new_slice[yax].step.imag))
            sl = sl.to_frb(width = width, resolution = resolution,
                           height = height)
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

    def _create_point(self, point_tuple):
        coord = [self._spec_to_value(p) for p in point_tuple]
        return self.ds.point(coord)

    def _create_ray(self, ray_slice):
        start_point = [self._spec_to_value(v) for v in ray_slice.start]
        end_point = [self._spec_to_value(v) for v in ray_slice.stop]
        if getattr(ray_slice.step, "imag", 0.0) != 0.0:
            return LineBuffer(self.ds, start_point, end_point, 
                              int(ray_slice.step.imag))
        else:
            return self.ds.ray(start_point, end_point)

    def _create_ortho_ray(self, ray_tuple):
        axis = None
        new_slice = []
        coord = []
        npoints = 0
        start_point = []
        end_point = []
        for ax, v in enumerate(ray_tuple):
            if not isinstance(v, slice):
                val = self._spec_to_value(v)
                coord.append(val)
                new_slice.append(slice(None, None, None))
                start_point.append(val)
                end_point.append(val)
            else:
                if axis is not None: raise RuntimeError
                if getattr(v.step, "imag", 0.0) != 0.0:
                    npoints = int(v.step.imag)
                    xi = self._spec_to_value(v.start)
                    xf = self._spec_to_value(v.stop)
                    dx = (xf-xi)/npoints
                    start_point.append(xi+0.5*dx)
                    end_point.append(xf-0.5*dx)
                else:
                    axis = ax
                    new_slice.append(v)
        if npoints > 0:
            ray = LineBuffer(self.ds, start_point, end_point, npoints)
        else:
            if axis == 1:
                coord = [coord[1], coord[0]]
            source = self._create_region(new_slice)
            ray = self.ds.ortho_ray(axis, coord, data_source=source)
        return ray
