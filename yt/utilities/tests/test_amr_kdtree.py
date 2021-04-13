import itertools

import numpy as np

from yt.testing import assert_almost_equal, fake_amr_ds


def test_amr_kdtree_set_fields():
    ds = fake_amr_ds(fields=["density", "pressure"], units=["g/cm**3", "dyn/cm**2"])
    dd = ds.all_data()

    fields = ds.field_list
    dd.tiles.set_fields(fields, [True, True], False)
    gold = {}
    for i, block in enumerate(dd.tiles.traverse()):
        gold[i] = [data.copy() for data in block.my_data]

    for log_fields in itertools.product([True, False], [True, False]):
        dd.tiles.set_fields(fields, log_fields, False)
        for iblock, block in enumerate(dd.tiles.traverse()):
            for i in range(len(fields)):
                if log_fields[i]:
                    data = block.my_data[i]
                else:
                    data = np.log10(block.my_data[i])
                assert_almost_equal(gold[iblock][i], data)
