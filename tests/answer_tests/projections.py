from yt.utilities.answer_testing.output_tests import \
    SingleOutputTest, create_test
from yt.utilities.answer_testing.hydro_tests import \
    TestProjection, TestGasDistribution
from fields_to_test import field_list

for axis in range(3):
    for field in field_list:
        create_test(TestProjection, "projection_test_%s_%s" % (axis, field),
                    field = field, axis = axis)

for field in field_list:
    create_test(TestGasDistribution, "profile_density_test_%s" % field,
                field_x = "Density", field_y = field)


