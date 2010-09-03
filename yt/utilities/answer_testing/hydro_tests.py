import matplotlib; matplotlib.use("Agg")
import pylab
from yt.mods import *
from output_tests import SingleOutputTest, YTStaticOutputTest, create_test

class TestProjection(YTStaticOutputTest):

    field = None
    axis = None

    def run(self):
        # First we get our flattened projection -- this is the
        # Density, px, py, pdx, and pdy
        proj = self.pf.h.proj(self.axis, self.field)
        # Now let's stick it in a buffer
        pixelized_proj = self.pixelize(proj, self.field)
        # We just want the values, so this can be stored
        # independently of the parameter file.
        # The .data attributes strip out everything other than the actual array
        # values.
        self.result = (proj.data, pixelized_proj.data)

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
        fn = "%s_%s_projection.png" % (self.pf, self.field)
        pylab.savefig(fn)
        return [fn]

# Now we create all our tests.  We are using the create_test
# function, which is a relatively simple function that takes the base class,
# a name, and any parameters that the test requires.
for axis in range(3):
    for field in ["Density", "Temperature"]:
        create_test(TestProjection, "projection_test_%s_%s" % (axis, field),
                    field = field, axis = axis)

class TestGasDistribution(YTStaticOutputTest):
    field_x = None
    field_y = None
    weight = "CellMassMsun"
    n_bins = 32

    def run(self):
        # We're NOT going to use the low-level profiling API here,
        # because we are avoiding the calculations of min/max,
        # as those should be tested in another test.
        pc = PlotCollection(self.pf, center=self.sim_center)
        p = pc.add_profile_object(self.entire_simulation,
            [self.field_x, self.field_y], x_bins = self.n_bins,
            weight=self.weight)
        # The arrays are all stored in a dictionary hanging off the profile
        # object
        self.result = p.data._data
                    
    def compare(self, old_result):
        self.compare_data_arrays(
            self.result, old_result)

    def plot(self):
        return []

# Now we create all our tests, but we're only going to check the binning
# against Density for now.
for field in ["Temperature", "x-velocity"]:
    create_test(TestGasDistribution, "profile_density_test_%s" % field,
                field_x = "Density", field_y = field)
