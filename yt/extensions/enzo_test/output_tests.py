from yt.mods import *

test_registry = {}

class RegressionTestException(Exception):
    pass

class ValueDelta(RegressionTestException):
    def __init__(self, delta, acceptable):
        self.delta = delta
        self.acceptable = acceptable

    def __repr__(self):
        return "ValueDelta: Delta %0.5e, max of %0.5e" % (
            self.delta, self.acceptable)

class ArrayDelta(ValueDelta):
    def __repr__(self):
        return "ArrayDelta: Delta %0.5e, max of %0.5e" % (
            self.delta, self.acceptable)

class RegressionTest(object):
    name = None
    result = None
    output_type = None

    class __metaclass__(type):
        # This ensures that all the tests are auto-registered if they have a
        # name.

        def __init__(cls, name, b, d):
            type.__init__(cls, name, b, d)
            if cls.name is not None:
                test_registry[cls.name] = cls

    def setup(self):
        pass

    def run(self):
        pass

    def compare(self, old_result):
        pass

    def plot(self):
        pass

    def compare_array_delta(self, a1, a2, acceptable):
        delta = na.abs(a1 - a2)/(a1 + a2)
        if delta.max() > acceptable:
            raise ArrayDelta(delta, acceptable)
        return True

    def compare_value_delta(self, v1, v2, acceptable):
        delta = na.abs(v1 - v2)/(v1 + v2)
        if delta > acceptable:
            raise ValueDelta(delta, acceptable)
        return True

class SingleOutputTest(RegressionTest):
    output_type = 'single'

    def __init__(self, filename):
        self.filename = filename

class MultipleOutputTest(RegressionTest):
    output_type = 'multiple'

    io_log_header = "DATASET WRITTEN"

    def __init__(self, io_log):
        self.io_log = io_log

    def __iter__(self):
        for line in open(self.io_log):
            yield line[len(self.io_log_header):].strip()

def create_test(base, new_name, **attrs):
    new_name = "%s_%s" % (base.__name__, new_name)
    attrs['name'] = new_name
    return type(new_name, (base,), attrs)

class YTStaticOutputTest(SingleOutputTest):

    def setup(self):
        self.pf = load(self.filename)

    def pixelize(self, data, field, edges = None, dims = (512, 512)):
        xax = lagos.x_dict[self.axis]
        yax = lagos.y_dict[self.axis]
        
        if edges is None:
            edges = (self.pf["DomainLeftEdge"][xax],
                     self.pf["DomainRightEdge"][xax],
                     self.pf["DomainLeftEdge"][yax],
                     self.pf["DomainRightEdge"][yax])
        frb = raven.FixedResolutionBuffer( data, edges, dims)
        frb[field] # To make the pixelization
        return frb

    def compare_data_arrays(self, d1, d2, tol = 1e-7):
        for field in d1.keys():
            self.compare_array_delta(d1[field], d2[field], tol)

    @property
    def sim_center(self):
        return 0.5*(self.pf["DomainRightEdge"] + self.pf["DomainLeftEdge"])

    @property
    def max_dens_location(self):
        return self.pf.h.find_max("Density")[1]

    @property
    def entire_simulation(self):
        return self.pf.h.all_data()
        

