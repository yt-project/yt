import numpy as np

import yt
from yt.frontends.ramses.hilbert import get_cpu_list, hilbert3d
from yt.testing import assert_equal, requires_file


def test_hilbert3d():
    # 8 different cases, checked against RAMSES' own implementation
    inputs = [
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
        [1, 1, 0],
        [0, 0, 1],
        [1, 0, 1],
        [0, 1, 1],
        [1, 1, 1],
    ]
    outputs = [0, 1, 7, 6, 3, 2, 4, 5]

    for i, o in zip(inputs, outputs):
        assert_equal(int(hilbert3d(i, 3)), o)


output_00080 = "output_00080/info_00080.txt"


@requires_file(output_00080)
def test_get_cpu_list():
    ds = yt.load(output_00080)

    np.random.seed(16091992)
    # These are randomly generated outputs, checked against RAMSES' own implementation
    inputs = (
        [[0.27747276, 0.30018937, 0.17916189], [0.42656026, 0.40509483, 0.29927838]],
        [[0.90660856, 0.44201328, 0.22770587], [1.09175462, 0.58017918, 0.2836648]],
        [[0.98542323, 0.58543376, 0.45858327], [1.04441105, 0.62079207, 0.58919283]],
        [[0.42274841, 0.44887745, 0.87793679], [0.52066634, 0.58936331, 1.00666222]],
        [[0.69964803, 0.65893669, 0.03660775], [0.80565696, 0.67409752, 0.11434604]],
    )
    outputs = ([0, 15], [0, 15], [0, 1, 15], [0, 13, 14, 15], [0])

    for i, o in zip(inputs, outputs):
        bbox = i
        ls = get_cpu_list(ds, bbox)
        assert len(ls) > 0
        assert all(np.array(o) == np.array(ls))
