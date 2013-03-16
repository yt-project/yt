from yt.testing import fake_random_pf, assert_equal
from yt.analysis_modules.level_sets.api import identify_contours
import yt.data_objects.api 
import cPickle
import os

def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"


def test_save_load_pickle():

    pf = fake_random_pf(64)

    # create extracted region from boolean (fairly complex object)
    center = (pf.domain_left_edge + pf.domain_right_edge) / 2
    sp_outer = pf.h.sphere(center, pf.domain_width[0])
    sp_inner = pf.h.sphere(center, pf.domain_width[0]/10.0)
    sp_boolean = pf.h.boolean([sp_outer, "NOT", sp_inner])

    minv, maxv = sp_boolean.quantities["Extrema"]("Density")[0]
    contour_threshold = min(minv*10, 0.9*maxv)
    
    contours = sp_boolean.extract_connected_sets("Density", 1, contour_threshold, maxv+1, log_space=True, cache=True)

    # save object
    cPickle.dump(contours[1][0], open("myobject.cpkl", "wb"))
    
    # load object
    test_load = cPickle.load(open("myobject.cpkl", "rb"))

    yield assert_equal, test_load != None, True
    yield assert_equal, len(contours[1][0]), len(test_load)

    os.remove("myobject.cpkl")
