"""
Base classes for answer testing

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2010-2011 Matthew Turk.  All Rights Reserved.

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

import matplotlib
from yt.mods import *

# We first create our dictionary of tests to run.  This starts out empty, and
# as tests are imported it will be filled.
if "TestRegistry" not in locals():
    print "Initializing TestRegistry"
    class TestRegistry(dict):
        def __new__(cls, *p, **k):
            if not '_the_instance' in cls.__dict__:
                cls._the_instance = dict.__new__(cls)
                return cls._the_instance
if "test_registry" not in locals():
    print "Initializing test_registry"
    test_registry = TestRegistry()

# The exceptions we raise, related to the character of the failure.

class RegressionTestException(Exception):
    pass

class ValueDelta(RegressionTestException):
    def __init__(self, delta, acceptable):
        self.delta = delta
        self.acceptable = acceptable

    def __repr__(self):
        return "ValueDelta: Delta %s, max of %s" % (
            self.delta, self.acceptable)

class ArrayDelta(ValueDelta):
    def __repr__(self):
        nabove = len(na.where(self.delta > self.acceptable)[0])
        return "ArrayDelta: Delta max of %s, acceptable of %s.\n" \
               "%d of %d points above the acceptable limit" % \
               (na.nanmax(self.delta), self.acceptable, nabove,
                self.delta.size)

class ShapeMismatch(RegressionTestException):
    def __init__(self, old_shape, current_shape):
        self.old_shape = old_shape
        self.current_shape = current_shape

    def __repr__(self):
        return "Shape Mismatch: old_buffer %s, current_buffer %s" % (
            self.old_shape, self.current_shape)

class RegressionTest(object):
    name = None
    result = None
    output_type = None

    class __metaclass__(type):
        # This ensures that all the tests are auto-registered if they have a
        # name.  If they do not have a name, they are considered to be base
        # classes to be overridden and implemented by someone else.
        def __init__(cls, name, b, d):
            type.__init__(cls, name, b, d)
            if cls.name is not None:
                test_registry[cls.name] = cls

    def setup(self):
        """
        This function must be defined if the problem requires additional setup.
        Note that for the most part this will be defined in base classes where
        subclasses will only implement 'run'.
        """
        pass

    def run(self):
        """
        This function must generate a result value, of any type, and store it
        in self.result.
        """
        pass

    def compare(self, old_result):
        """
        This function must accept `old_result` and compare it somehow against
        the value stored in `self.result`.  If the result is a failure, it must
        raise an exception.  Otherwise it is considered to be a success.
        """
        pass

    def plot(self):
        """
        This function can optionally plot the contents of `self.result`.
        """
        pass

    def compare_array_delta(self, a1, a2, acceptable):
        """
        This is a helper function.  It accepts two numpy arrays and compares
        the maximum relative difference.  If the maximum relative difference is
        greater than `acceptable` it is considered a failure and an appropriate
        exception is raised.
        """
        if a1.shape != a2.shape:
            raise ShapeMismatch(a1, a2)
        delta = na.abs(a1 - a2).astype("float64")/(a1 + a2)
        if na.nanmax(delta) > acceptable:
            raise ArrayDelta(delta, acceptable)
        return True

    def compare_value_delta(self, v1, v2, acceptable):
        """
        This is a helper function.  It accepts two floating point values and
        calculates the maximum relative different.  If the maximum relative
        difference is greater than `acceptable` it is considered a failure and
        an appropriate exception is raised.
        """
        delta = na.abs(v1 - v2)/(v1 + v2)
        if delta > acceptable:
            raise ValueDelta(delta, acceptable)
        return True

class SingleOutputTest(RegressionTest):
    output_type = 'single'

    def __init__(self, filename):
        """
        This test mechanism is designed to accept a single filename and
        evaluate it, not necessarily utilizing yt's functionality to do so.
        """
        self.filename = filename

class MultipleOutputTest(RegressionTest):
    output_type = 'multiple'

    io_log_header = "DATASET WRITTEN"

    def __init__(self, io_log):
        """
        This test mechanism is designed to accept an OutputLog file and then
        iterate over it, evaluating every single dataset individually.
        """
        self.io_log = io_log

    def __iter__(self):
        if isinstance(self.io_log, types.StringTypes):
            for line in open(self.io_log):
                yield line[len(self.io_log_header):].split()[0].strip()
        elif isinstance(self.io_log, types.ListType):
            for line in self.io_log: yield line

def create_test(base, new_name, **attrs):
    """
    This function accepts a base class of a test, sets some attributes on it,
    and then registers a new test.  It's a fast way of registering multiple
    tests that share the same testing logic but that differ on a few parameters
    or combinations of parameters.
    """
    new_name = "%s_%s" % (base.__name__, new_name)
    attrs['name'] = new_name
    return type(new_name, (base,), attrs)

class YTStaticOutputTest(SingleOutputTest):

    def setup(self):
        self.pf = load(self.filename)

    def pixelize(self, data, field, edges = None, dims = (512, 512)):
        """
        This is a helper function that returns a 2D array of the specified
        source, in the specified field, at the specified spatial extent.
        """
        xax = x_dict[self.axis]
        yax = y_dict[self.axis]
        
        if edges is None:
            edges = (self.pf.domain_left_edge[xax],
                     self.pf.domain_right_edge[xax],
                     self.pf.domain_left_edge[yax],
                     self.pf.domain_right_edge[yax])
        frb = FixedResolutionBuffer( data, edges, dims)
        frb[field] # To make the pixelization
        return frb

    def compare_data_arrays(self, d1, d2, tol = 1e-7):
        """
        This is a helper function.  It accepts two dictionaries of numpy arrays
        and compares the maximum relative difference of every array.  If the
        maximum relative difference is greater than `acceptable` it is
        considered a failure and an appropriate exception is raised.
        """
        for field in d1.keys():
            self.compare_array_delta(d1[field], d2[field], tol)

    @property
    def sim_center(self):
        """
        This returns the center of the domain.
        """
        return 0.5*(self.pf.domain_right_edge + self.pf.domain_left_edge)

    @property
    def max_dens_location(self):
        """
        This is a helper function to return the location of the most dense
        point.
        """
        return self.pf.h.find_max("Density")[1]

    @property
    def entire_simulation(self):
        """
        Return an unsorted array of values that cover the entire domain.
        """
        return self.pf.h.all_data()
        

