import hashlib
import numpy as na

from yt.utilities.answer_testing.output_tests import \
    YTStaticOutputTest, RegressionTestException, create_test
from yt.funcs import ensure_list, iterable
from fields_to_test import field_list, particle_field_list


class FieldHashesDontMatch(RegressionTestException):
    pass

known_objects = {}


def register_object(func):
    known_objects[func.func_name] = func
    return func


@register_object
def centered_sphere(tobj):
    center = 0.5 * (tobj.pf.domain_right_edge + tobj.pf.domain_left_edge)
    width = (tobj.pf.domain_right_edge - tobj.pf.domain_left_edge).max()
    tobj.data_object = tobj.pf.h.sphere(center, width / 0.25)


@register_object
def off_centered_sphere(tobj):
    center = 0.5 * (tobj.pf.domain_right_edge + tobj.pf.domain_left_edge)
    width = (tobj.pf.domain_right_edge - tobj.pf.domain_left_edge).max()
    tobj.data_object = tobj.pf.h.sphere(center - 0.25 * width, width / 0.25)


@register_object
def corner_sphere(tobj):
    width = (tobj.pf.domain_right_edge - tobj.pf.domain_left_edge).max()
    tobj.data_object = tobj.pf.h.sphere(tobj.pf.domain_left_edge, width / 0.25)


@register_object
def disk(self):
    center = (self.pf.domain_right_edge + self.pf.domain_left_edge) / 2.
    radius = (self.pf.domain_right_edge - self.pf.domain_left_edge).max() / 10.
    height = (self.pf.domain_right_edge - self.pf.domain_left_edge).max() / 10.
    normal = na.array([1.] * 3)
    self.data_object = self.pf.h.disk(center, normal, radius, height)


@register_object
def all_data(self):
    self.data_object = self.pf.h.all_data()

_new_known_objects = {}
for field in ["Density"]:  # field_list:
    for object_name in known_objects:

        def _rfunc(oname, fname):

            def func(tobj):
                known_objects[oname](tobj)
                tobj.orig_data_object = tobj.data_object
                avg_value = tobj.orig_data_object.quantities[
                        "WeightedAverageQuantity"](fname, "Density")
                tobj.data_object = tobj.orig_data_object.cut_region(
                        ["grid['%s'] > %s" % (fname, avg_value)])
            return func
        _new_known_objects["%s_cut_region_%s" % (object_name, field)] = \
                _rfunc(object_name, field)
known_objects.update(_new_known_objects)


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


class YTExtractIsocontoursTest(YTFieldValuesTest):

    def run(self):
        val = self.data_object.quantities["WeightedAverageQuantity"](
            "Density", "Density")
        rset = self.data_object.extract_isocontours("Density",
            val, rescale=False, sample_values="Temperature")
        self.result = rset

    def compare(self, old_result):
        if self.result[0].size == 0 and old_result[0].size == 0:
            return True
        self.compare_array_delta(self.result[0].ravel(),
                                 old_result[0].ravel(), 1e-7)
        self.compare_array_delta(self.result[1], old_result[1], 1e-7)


class YTIsocontourFluxTest(YTFieldValuesTest):

    def run(self):
        val = self.data_object.quantities["WeightedAverageQuantity"](
            "Density", "Density")
        flux = self.data_object.calculate_isocontour_flux(
           "Density", val, "x-velocity", "y-velocity", "z-velocity")
        self.result = flux

    def compare(self, old_result):
        self.compare_value_delta(self.result, old_result, 1e-7)

for object_name in known_objects:
    for field in field_list + particle_field_list:
        if "cut_region" in object_name and field in particle_field_list:
            continue
        create_test(YTFieldValuesTest, "%s_%s" % (object_name, field),
                    field=field, object_name=object_name)
    create_test(YTExtractIsocontoursTest, "%s" % (object_name),
                object_name=object_name)
    create_test(YTIsocontourFluxTest, "%s" % (object_name),
                object_name=object_name)


class YTDerivedQuantityTest(YTStaticOutputTest):

    def setup(self):
        YTStaticOutputTest.setup(self)
        known_objects[self.object_name](self)

    def compare(self, old_result):
        if hasattr(self.result, 'tostring'):
            self.compare_array_delta(self.result, old_result, 1e-7)
            return
        elif iterable(self.result):
            a1 = na.array(self.result)
            a2 = na.array(old_result)
            self.compare_array_delta(a1, a2, 1e-7)
        else:
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
        # Some special exceptions
        if "cut_region" in object_name and (
            "SpinParameter" in dq or
            "TotalMass" in dq):
            continue
        create_test(YTDerivedQuantityTest, "%s_%s" % (object_name, dq),
                    dq_name=dq, object_name=object_name)


class YTDerivedQuantityTestField(YTDerivedQuantityTest):

    def run(self):
        self.result = self.data_object.quantities[self.dq_name](
            self.field_name)

for object_name in known_objects:
    for field in field_list:
        for dq in ["Extrema", "TotalQuantity", "MaxLocation", "MinLocation"]:
            create_test(YTDerivedQuantityTestField,
                        "%s_%s" % (object_name, field),
                        field_name=field, dq_name=dq,
                        object_name=object_name)


class YTDerivedQuantityTest_WeightedAverageQuantity(YTDerivedQuantityTest):

    def run(self):
        self.result = self.data_object.quantities["WeightedAverageQuantity"](
            self.field_name, weight="CellMassMsun")

for object_name in known_objects:
    for field in field_list:
        create_test(YTDerivedQuantityTest_WeightedAverageQuantity,
                    "%s_%s" % (object_name, field),
                    field_name=field,
                    object_name=object_name)
