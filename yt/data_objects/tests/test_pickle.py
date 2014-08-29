"""
Testsuite for pickling yt objects.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
from yt.extern.six.moves import cPickle
import os
import tempfile
from yt.testing \
    import fake_random_ds, assert_equal


def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"


def test_save_load_pickle():
    """Main test for loading pickled objects"""
    return # Until boolean regions are implemented we can't test this
    test_ds = fake_random_ds(64)

    # create extracted region from boolean (fairly complex object)
    center = (test_ds.domain_left_edge + test_ds.domain_right_edge) / 2
    sp_outer = test_ds.sphere(center, test_ds.domain_width[0])
    sp_inner = test_ds.sphere(center, test_ds.domain_width[0] / 10.0)
    sp_boolean = test_ds.boolean([sp_outer, "NOT", sp_inner])

    minv, maxv = sp_boolean.quantities["Extrema"]("density")[0]
    contour_threshold = min(minv * 10.0, 0.9 * maxv)

    contours = sp_boolean.extract_connected_sets(
        "density", 1, contour_threshold, maxv + 1, log_space=True, cache=True)

    # save object
    cpklfile = tempfile.NamedTemporaryFile(delete=False)
    cPickle.dump(contours[1][0], cpklfile)
    cpklfile.close()

    # load object
    test_load = cPickle.load(open(cpklfile.name, "rb"))

    assert_equal.description = \
        "%s: File was pickle-loaded succesfully" % __name__
    yield assert_equal, test_load is not None, True
    assert_equal.description = \
        "%s: Length of pickle-loaded connected set object" % __name__
    yield assert_equal, len(contours[1][0]), len(test_load)

    os.remove(cpklfile.name)
