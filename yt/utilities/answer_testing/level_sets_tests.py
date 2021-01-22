import numpy as np

from yt.testing import assert_equal
from yt.utilities.answer_testing.framework import AnswerTestingTest


class ExtractConnectedSetsTest(AnswerTestingTest):
    _type_name = "ExtractConnectedSets"
    _attrs = ()

    def __init__(self, ds_fn, data_source, field, num_levels, min_val, max_val):
        super().__init__(ds_fn)
        self.data_source = data_source
        self.field = field
        self.num_levels = num_levels
        self.min_val = min_val
        self.max_val = max_val

    def run(self):
        n, all_sets = self.data_source.extract_connected_sets(
            self.field, self.num_levels, self.min_val, self.max_val
        )
        result = []
        for level in all_sets:
            for set_id in all_sets[level]:
                result.append(
                    [
                        all_sets[level][set_id]["cell_mass"].size,
                        all_sets[level][set_id]["cell_mass"].sum(),
                    ]
                )
        result = np.array(result)
        return result

    def compare(self, new_result, old_result):
        err_msg = f"Size and/or mass of connected sets do not agree for {self.ds_fn}."
        assert_equal(new_result, old_result, err_msg=err_msg, verbose=True)
