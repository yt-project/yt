from yt.utilities.answer_testing.output_tests import \
    SingleOutputTest, create_test
from yt.utilities.answer_testing.hydro_tests import \
    TestProjection, TestOffAxisProjection, TestSlice, \
    TestRay, TestGasDistribution, Test2DGasDistribution

from fields_to_test import field_list

for field in field_list:
    create_test(TestRay, "%s" % field, field=field)

for axis in range(3):
    for field in field_list:
        create_test(TestSlice, "%s_%s" % (axis, field),
                    field=field, axis=axis)

for axis in range(3):
    for field in field_list:
        create_test(TestProjection, "%s_%s" % (axis, field),
                    field=field, axis=axis)
        create_test(TestProjection, "%s_%s_Density" % (axis, field),
                    field=field, axis=axis, weight_field="Density")

for field in field_list:
    create_test(TestOffAxisProjection, "%s_%s" % (axis, field),
                field=field, axis=axis)
    create_test(TestOffAxisProjection, "%s_%s_Density" % (axis, field),
                field=field, axis=axis, weight_field="Density")

for field in field_list:
    if field != "Density":
        create_test(TestGasDistribution, "density_%s" % field,
                    field_x="Density", field_y=field)
    if field not in ("x-velocity", "Density"):
        create_test(Test2DGasDistribution, "density_x-vel_%s" % field,
                    field_x="Density", field_y="x-velocity", field_z=field,
                    weight="CellMassMsun")
