from yt.utilities.answer_testing.output_tests import \
    SingleOutputTest, create_test
from yt.utilities.answer_testing.halo_tests import \
    TestHaloCompositionHashHOP, TestHaloCompositionHashFOF, \
    TestHaloCompositionHashPHOP 

create_test(TestHaloCompositionHashHOP, "halo_composition_test_hash_HOP", threshold=80.0)

create_test(TestHaloCompositionHashFOF, "halo_composition_test_hash_FOF", threshold=80.0)

create_test(TestHaloCompositionHashPHOP, "halo_composition_test_hash_PHOP", threshold=80.0)
