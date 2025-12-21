import numpy as np
from numpy.testing import assert_equal

from yt.testing import fake_random_ds

def test_chained_cut_region_with_locals():
    ds = fake_random_ds(16, fields=("density",), units=("g/cm**3",))
    ad = ds.all_data()
    cut_1 = ad.cut_region("obj['density'] > density_min", locals={"density_min": 0.2})

    # This chained cut region should inherit 'density_min' from cut_1's locals
    cut_2 = cut_1.cut_region("obj['density'] < density_max", locals={"density_max": 0.8})

    # Verify locals are merged
    assert "density_min" in cut_2.locals
    assert "density_max" in cut_2.locals
    assert cut_2.locals["density_min"] == 0.2
    assert cut_2.locals["density_max"] == 0.8

    mask = (ad["density"] > 0.2) & (ad["density"] < 0.8)
    expected = np.sort(ad["density"][mask])

    # This access triggers evaluation involving locals
    result = np.sort(cut_2["density"])

    assert_equal(expected, result)
