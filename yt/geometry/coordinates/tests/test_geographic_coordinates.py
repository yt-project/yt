# Some tests for the geographic coordinates handler

import numpy as np

from yt.testing import \
    fake_amr_ds, \
    assert_equal, \
    assert_rel_equal

# Our canonical tests are that we can access all of our fields and we can
# compute our volume correctly.

def test_geographic_coordinates():
    # We're going to load up a simple AMR grid and check its volume
    # calculations and path length calculations.

    # Note that we are setting it up to have an altitude of 1000 maximum, which
    # means our volume will be that of a shell 1000 wide, starting at r of
    # whatever our surface_height is set to.
    ds = fake_amr_ds(geometry="geographic")
    ds.surface_height = ds.quan(5000, "code_length")
    axes = ["latitude", "longitude", "altitude"]
    for i, axis in enumerate(axes):
        dd = ds.all_data()
        fi = ("index", axis)
        fd = ("index", "d%s" % axis)
        ma = np.argmax(dd[fi])
        yield assert_equal, dd[fi][ma] + dd[fd][ma] / 2.0, ds.domain_right_edge[i].d
        mi = np.argmin(dd[fi])
        yield assert_equal, dd[fi][mi] - dd[fd][mi] / 2.0, ds.domain_left_edge[i].d
        yield assert_equal, dd[fd].max(), (ds.domain_width/ds.domain_dimensions)[i].d
    inner_r = ds.surface_height
    outer_r = ds.surface_height + ds.domain_width[2]
    yield assert_equal, dd["index","dtheta"], \
                        dd["index","dlatitude"]*np.pi/180.0
    yield assert_equal, dd["index","dphi"], \
                        dd["index","dlongitude"]*np.pi/180.0
    # Note our terrible agreement here.
    yield assert_rel_equal, dd["cell_volume"].sum(dtype="float64"), \
                        (4.0/3.0) * np.pi * (outer_r**3 - inner_r**3), \
                        3
    yield assert_equal, dd["index", "path_element_altitude"], \
                        dd["index", "daltitude"]
    yield assert_equal, dd["index", "path_element_altitude"], \
                        dd["index", "dr"]
    # Note that latitude corresponds to theta, longitude to phi
    yield assert_equal, dd["index", "path_element_latitude"], \
                        dd["index", "r"] * \
                        dd["index", "dlatitude"] * np.pi/180.0
    yield assert_equal, dd["index", "path_element_longitude"], \
                        dd["index", "r"] * \
                        dd["index", "dlongitude"] * np.pi/180.0 * \
                        np.sin((dd["index", "latitude"] + 90.0) * np.pi/180.0)
    # We also want to check that our radius is correct
    yield assert_equal, dd["index","r"], \
                        dd["index","altitude"] + ds.surface_height

def test_internal_geographic_coordinates():
    # We're going to load up a simple AMR grid and check its volume
    # calculations and path length calculations.

    # Note that we are setting it up to have depth of 1000 maximum, which
    # means our volume will be that of a shell 1000 wide, starting at r of
    # outer_radius - 1000.
    ds = fake_amr_ds(geometry="internal_geographic")
    ds.outer_radius = ds.quan(5000, "code_length")
    axes = ["latitude", "longitude", "depth"]
    for i, axis in enumerate(axes):
        dd = ds.all_data()
        fi = ("index", axis)
        fd = ("index", "d%s" % axis)
        ma = np.argmax(dd[fi])
        yield assert_equal, dd[fi][ma] + dd[fd][ma] / 2.0, ds.domain_right_edge[i].d
        mi = np.argmin(dd[fi])
        yield assert_equal, dd[fi][mi] - dd[fd][mi] / 2.0, ds.domain_left_edge[i].d
        yield assert_equal, dd[fd].max(), (ds.domain_width/ds.domain_dimensions)[i].d
    inner_r = ds.outer_radius - ds.domain_right_edge[2]
    outer_r = ds.outer_radius
    yield assert_equal, dd["index","dtheta"], \
                        dd["index","dlatitude"]*np.pi/180.0
    yield assert_equal, dd["index","dphi"], \
                        dd["index","dlongitude"]*np.pi/180.0
    # Note our terrible agreement here.
    yield assert_rel_equal, dd["cell_volume"].sum(dtype="float64"), \
                        (4.0/3.0) * np.pi * (outer_r**3 - inner_r**3), \
                        3
    yield assert_equal, dd["index", "path_element_depth"], \
                        dd["index", "ddepth"]
    yield assert_equal, dd["index", "path_element_depth"], \
                        dd["index", "dr"]
    # Note that latitude corresponds to theta, longitude to phi
    yield assert_equal, dd["index", "path_element_latitude"], \
                        dd["index", "r"] * \
                        dd["index", "dlatitude"] * np.pi/180.0
    yield assert_equal, dd["index", "path_element_longitude"], \
                        dd["index", "r"] * \
                        dd["index", "dlongitude"] * np.pi/180.0 * \
                        np.sin((dd["index", "latitude"] + 90.0) * np.pi/180.0)
    # We also want to check that our radius is correct
    yield assert_equal, dd["index","r"], \
                        -1.0*dd["index","depth"] + ds.outer_radius
