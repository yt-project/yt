from yt.utilities.answer_testing.output_tests import \
    SingleOutputTest, create_test
from yt.utilities.answer_testing.boolean_region_tests import \
    TestBooleanANDGridQuantity, TestBooleanORGridQuantity, \
    TestBooleanNOTGridQuantity, TestBooleanANDParticleQuantity, \
    TestBooleanORParticleQuantity, TestBooleanNOTParticleQuantity

create_test(TestBooleanANDGridQuantity, "BooleanANDGrid")

create_test(TestBooleanORGridQuantity, "BooleanORGrid")

create_test(TestBooleanNOTGridQuantity, "BooleanNOTGrid")

create_test(TestBooleanANDParticleQuantity, "BooleanANDParticle")

create_test(TestBooleanORParticleQuantity, "BooleanORParticle")

create_test(TestBooleanNOTParticleQuantity, "BooleanNOTParticle")
