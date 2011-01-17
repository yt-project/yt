import hashlib
import numpy as na

from yt.utilities.answer_testing.output_tests import \
    YTStaticOutputTest, RegressionTestException, create_test
from yt.funcs import ensure_list
from fields_to_test import field_list, particle_field_list

class FieldHashesDontMatch(RegressionTestException):
    pass

class YTFieldValuesTest(YTStaticOutputTest):
    def run(self):
        vals = self.data_object[self.field].copy()
        vals.sort()
        self.result = hashlib.sha256(vals.tostring()).hexdigest()

    def compare(self, old_result):
        if self.result != old_result: raise FieldHashesDontMatch

class CenteredSphere(YTFieldValuesTest):

    def setup(self):
        YTFieldValuesTest.setup(self)
        center = 0.5*(self.pf.domain_right_edge + self.pf.domain_left_edge)
        width = (self.pf.domain_right_edge - self.pf.domain_left_edge).max()
        self.data_object = self.pf.h.sphere(center, width/0.25)

for field in field_list + particle_field_list:
    create_test(CenteredSphere, "centered_sphere_%s" % (field), field = field)

class OffCenteredSphere(YTFieldValuesTest):

    def setup(self):
        YTFieldValuesTest.setup(self)
        center = 0.5*(self.pf.domain_right_edge + self.pf.domain_left_edge)
        width = (self.pf.domain_right_edge - self.pf.domain_left_edge).max()
        self.data_object = self.pf.h.sphere(center - 0.25 * width, width/0.25)

for field in field_list + particle_field_list:
    create_test(OffCenteredSphere, "off_centered_sphere_%s" % (field), field = field)

class CornerSphere(YTFieldValuesTest):

    def setup(self):
        YTFieldValuesTest.setup(self)
        width = (self.pf.domain_right_edge - self.pf.domain_left_edge).max()
        self.data_object = self.pf.h.sphere(self.pf.domain_left_edge, width/0.25)

for field in field_list + particle_field_list:
    create_test(CornerSphere, "corner_sphere_%s" % (field), field = field)

class AllData(YTFieldValuesTest):
    def setup(self):
        YTFieldValuesTest.setup(self)
        self.data_object = self.pf.h.all_data()

for field in field_list + particle_field_list:
    create_test(AllData, "all_data_%s" % (field), field = field)
