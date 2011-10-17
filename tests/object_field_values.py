import hashlib
import numpy as na

from yt.utilities.answer_testing.output_tests import \
    YTStaticOutputTest, RegressionTestException, create_test
from yt.funcs import ensure_list
from fields_to_test import field_list, particle_field_list

class FieldHashesDontMatch(RegressionTestException):
    pass

known_objects = {}

def register_object(func):
    known_objects[func.func_name] = func
    return func

@register_object
def centered_sphere(self):
    center = 0.5*(self.pf.domain_right_edge + self.pf.domain_left_edge)
    width = (self.pf.domain_right_edge - self.pf.domain_left_edge).max()
    self.data_object = self.pf.h.sphere(center, width/0.25)

@register_object
def off_centered_sphere(self):
    center = 0.5*(self.pf.domain_right_edge + self.pf.domain_left_edge)
    width = (self.pf.domain_right_edge - self.pf.domain_left_edge).max()
    self.data_object = self.pf.h.sphere(center - 0.25 * width, width/0.25)

@register_object
def corner_sphere(self):
    width = (self.pf.domain_right_edge - self.pf.domain_left_edge).max()
    self.data_object = self.pf.h.sphere(self.pf.domain_left_edge, width/0.25)

@register_object
def disk(self):
    center = (self.pf.domain_right_edge + self.pf.domain_left_edge)/2.
    radius = (self.pf.domain_right_edge - self.pf.domain_left_edge).max()/10.
    height = (self.pf.domain_right_edge - self.pf.domain_left_edge).max()/10.
    normal = na.array([1.]*3)
    self.data_object = self.pf.h.disk(center, normal, radius, height)
    
@register_object
def all_data(self):
    self.data_object = self.pf.h.all_data()

class YTFieldValuesTest(YTStaticOutputTest):
    def run(self):
        vals = self.data_object[self.field].copy()
        vals.sort()
        self.result = hashlib.sha256(vals.tostring()).hexdigest()

    def compare(self, old_result):
        if self.result != old_result: raise FieldHashesDontMatch

    def setup(self):
        YTStaticOutputTest.setup(self)
        known_objects[self.object_name](self)

for object_name in known_objects:
    for field in field_list + particle_field_list:
        create_test(YTFieldValuesTest, "%s_%s" % (object_name, field),
                    field = field, object_name = object_name)

class YTDerivedQuantityTest(YTStaticOutputTest):
    def setup(self):
        YTStaticOutputTest.setup(self)
        known_objects[self.object_name](self)

    def compare(self, old_results):
        if self.result != old_result: raise FieldHashesDontMatch

    def run(self):
        # This only works if it takes no arguments
        self.result = self.data_object.quantities[self.dq_name]()

dq_names = ["TotalMass", "AngularMomentumVector", "CenterOfMass",
            "BulkVelocity", "BaryonSpinParameter", "ParticleSpinParameter"]

# Extrema, WeightedAverageQuantity, TotalQuantity, MaxLocation,
# MinLocation

for object_name in known_objects:
    for dq in dq_names:
        create_test(YTDerivedQuantityTest, "%s_%s" % (object_name, dq),
                    dq_name = dq, object_name = object_name)

class YTDerivedQuantityTestField(YTDerivedQuantityTest):
    def run(self):
        self.result = self.data_object.quantities[self.dq_name](
            self.field_name)

for object_name in known_objects:
    for field in field_list:
        for dq in ["Extrema", "TotalQuantity", "MaxLocation", "MinLocation"]:
            create_test(YTDerivedQuantityTestField,
                        "%s_%s" % (object_name, field),
                        field_name = field, dq_name = dq,
                        object_name = object_name)

class YTDerivedQuantityTest_WeightedAverageQuantity(YTDerivedQuantityTest):
    def run(self):
        self.result = self.data_object.quantities["WeightedAverageQuantity"](
            self.field_name, weight="CellMassMsun")

for object_name in known_objects:
    for field in field_list:
        create_test(YTDerivedQuantityTest_WeightedAverageQuantity,
                    "%s_%s" % (object_name, field),
                    field_name = field, 
                    object_name = object_name)

