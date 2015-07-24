"""
Hydro tests



"""
from __future__ import absolute_import

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import matplotlib
import pylab
from yt.mods import *
from .output_tests import SingleOutputTest, YTDatasetTest, create_test

class TestProjection(YTDatasetTest):

    field = None
    axis = None
    weight_field = None

    def run(self):
        # First we get our flattened projection -- this is the
        # Density, px, py, pdx, and pdy
        proj = self.ds.proj(self.field, self.axis, 
                              weight_field=self.weight_field)
        # Now let's stick it in a buffer
        pixelized_proj = self.pixelize(proj, self.field)
        # We just want the values, so this can be stored
        # independently of the dataset.
        # The .field_data attributes strip out everything other than the actual array
        # values.
        self.result = (proj.field_data, pixelized_proj.data)

    def compare(self, old_result):
        proj, pixelized_proj = self.result
        oproj, opixelized_proj = old_result

        self.compare_data_arrays(proj, oproj)
        self.compare_array_delta(
            pixelized_proj[self.field],
            opixelized_proj[self.field],
            1e-7)

    def plot(self):
        pylab.clf()
        pylab.imshow(self.result[1][self.field],
            interpolation='nearest', origin='lower')
        fn = "%s_%s_%s_projection.png" % (self.ds, self.field,
                                          self.weight_field)
        pylab.savefig(fn)
        return [fn]

class TestOffAxisProjection(YTDatasetTest):

    field = None
    weight_field = None

    def run(self):
        # Here proj will just be the data array.
        proj = off_axis_projection(self.ds, 
                                   (0.5 * (self.ds.domain_left_edge + 
                                           self.ds.domain_right_edge)),
                                   [1., 1., 1.], 1., 400,
                                   self.field, weight=self.weight_field)

        # values.
        self.result = proj

    def compare(self, old_result):
        proj  = self.result
        oproj = old_result

        self.compare_array_delta(proj, oproj, 1e-7)

    def plot(self):
        fn = "%s_%s_%s_off-axis_projection.png" % \
            (self.ds, self.field, self.weight_field)
        write_image(self.result, fn)
        return [fn]

class TestRay(YTDatasetTest):

    field = None

    def run(self):
        np.random.seed(4333)
        start_point = np.random.random(self.ds.dimensionality) * \
            (self.ds.domain_right_edge - self.ds.domain_left_edge) + \
            self.ds.domain_left_edge
        end_point   = np.random.random(self.ds.dimensionality) * \
            (self.ds.domain_right_edge - self.ds.domain_left_edge) + \
            self.ds.domain_left_edge

        # Here proj will just be the data array.
        ray = self.ds.ray(start_point, end_point, field=self.field)

        # values.
        self.result = ray[self.field]

    def compare(self, old_result):
        ray  = self.result
        oray = old_result

        self.compare_array_delta(ray, oray, 1e-7)

    def plot(self):
        return

class TestSlice(YTDatasetTest):

    field = None
    axis = None

    def run(self):
        # Here proj will just be the data array.
        slice = self.ds.slice(self.axis, 
                                (0.5 * (self.ds.domain_left_edge + 
                                        self.ds.domain_right_edge))[self.axis],
                                fields=self.field)
        # values.
        self.result = slice.field_data

    def compare(self, old_result):
        slice  = self.result
        oslice = old_result

        self.compare_data_arrays(slice, oslice)

    def plot(self):
        fn = "%s_%s_slice.png" % (self.ds, self.field)
        write_image(self.result[self.field], fn)
        return [fn]

# Now we create all our tests.  We are using the create_test
# function, which is a relatively simple function that takes the base class,
# a name, and any parameters that the test requires.
for axis in range(3):
    for field in ["density", "temperature"]:
        create_test(TestProjection, "projection_test_%s_%s" % (axis, field),
                    field = field, axis = axis)

class TestGasDistribution(YTDatasetTest):
    field_x = None
    field_y = None
    weight = "cell_mass"
    n_bins = 32

    def run(self):
        # We're NOT going to use the low-level profiling API here,
        # because we are avoiding the calculations of min/max,
        # as those should be tested in another test.
        pc = PlotCollection(self.ds, center=self.sim_center)
        p = pc.add_profile_object(self.entire_simulation,
            [self.field_x, self.field_y], x_bins = self.n_bins,
            weight=self.weight)
        # The arrays are all stored in a dictionary hanging off the profile
        # object
        self.result = p.data.field_data
                    
    def compare(self, old_result):
        self.compare_data_arrays(
            self.result, old_result)

    def plot(self):
        return []

# Now we create all our tests, but we're only going to check the binning
# against Density for now.
for field in ["temperature", "velocity_x"]:
    create_test(TestGasDistribution, "profile_density_test_%s" % field,
                field_x = "density", field_y = field)

class Test2DGasDistribution(TestGasDistribution):
    x_bins = 128
    y_bins = 128
    field_z = "cell_mass"
    weight = None
    def run(self):
        # We're NOT going to use the low-level profiling API here,
        # because we are avoiding the calculations of min/max,
        # as those should be tested in another test.
        pc = PlotCollection(self.ds, center=self.sim_center)
        p = pc.add_phase_object(self.entire_simulation,
            [self.field_x, self.field_y, self.field_z], x_bins = self.x_bins, y_bins = self.y_bins,
            weight=self.weight)
        # The arrays are all stored in a dictionary hanging off the profile
        # object
        self.result = p.data.field_data

