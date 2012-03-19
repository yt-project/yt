from yt.utilities.answer_testing.output_tests import \
    SingleOutputTest, create_test
from yt.utilities.answer_testing.halo_tests import \
    TestHaloCountHOP, TestHaloCountFOF, TestHaloCountPHOP

create_test(TestHaloCountHOP, "halo_count_HOP", threshold=80.0)

create_test(TestHaloCountFOF, "halo_count_FOF", link=0.2, padding=0.02)

create_test(TestHaloCountPHOP, "halo_count_PHOP", threshold=80.0)
