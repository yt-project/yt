import itertools as it
import os
import shutil
import tempfile
import unittest

import numpy as np
import pytest
import unyt
from numpy.testing import assert_equal

from yt.testing import (
    assert_fname,
    assert_rel_equal,
    cubicspline_python,
    fake_octree_ds,
    fake_random_ds,
    fake_sph_flexible_grid_ds,
    integrate_kernel,
)
from yt.visualization.api import (
    OffAxisProjectionPlot, 
    OffAxisSlicePlot,
    ProjectionPlot)
from yt.visualization.image_writer import write_projection
from yt.visualization.volume_rendering.api import off_axis_projection


# TODO: replace this with pytest.mark.parametrize
def expand_keywords(keywords, full=False):
    """
    expand_keywords is a means for testing all possible keyword
    arguments in the nosetests.  Simply pass it a dictionary of all the
    keyword arguments and all of the values for these arguments that you
    want to test.

    It will return a list of kwargs dicts containing combinations of
    the various kwarg values you passed it.  These can then be passed
    to the appropriate function in nosetests.

    If full=True, then every possible combination of keywords is produced,
    otherwise, every keyword option is included at least once in the output
    list.  Be careful, by using full=True, you may be in for an exponentially
    larger number of tests!

    Parameters
    ----------

    keywords : dict
        a dictionary where the keys are the keywords for the function,
        and the values of each key are the possible values that this key
        can take in the function

    full : bool
        if set to True, every possible combination of given keywords is
        returned

    Returns
    -------

    array of dicts
        An array of dictionaries to be individually passed to the appropriate
        function matching these kwargs.

    Examples
    --------

    >>> keywords = {}
    >>> keywords["dpi"] = (50, 100, 200)
    >>> keywords["cmap"] = ("cmyt.arbre", "cmyt.kelp")
    >>> list_of_kwargs = expand_keywords(keywords)
    >>> print(list_of_kwargs)

    array([{'cmap': 'cmyt.arbre', 'dpi': 50},
           {'cmap': 'cmyt.kelp', 'dpi': 100},
           {'cmap': 'cmyt.arbre', 'dpi': 200}], dtype=object)

    >>> list_of_kwargs = expand_keywords(keywords, full=True)
    >>> print(list_of_kwargs)

    array([{'cmap': 'cmyt.arbre', 'dpi': 50},
           {'cmap': 'cmyt.arbre', 'dpi': 100},
           {'cmap': 'cmyt.arbre', 'dpi': 200},
           {'cmap': 'cmyt.kelp', 'dpi': 50},
           {'cmap': 'cmyt.kelp', 'dpi': 100},
           {'cmap': 'cmyt.kelp', 'dpi': 200}], dtype=object)

    >>> for kwargs in list_of_kwargs:
    ...     write_projection(*args, **kwargs)
    """
    # if we want every possible combination of keywords, use iter magic
    if full:
        keys = sorted(keywords)
        list_of_kwarg_dicts = np.array(
            [
                dict(zip(keys, prod))
                for prod in it.product(*(keywords[key] for key in keys))
            ]
        )

    # if we just want to probe each keyword, but not necessarily every
    # combination
    else:
        # Determine the maximum number of values any of the keywords has
        num_lists = 0
        for val in keywords.values():
            if isinstance(val, str):
                num_lists = max(1.0, num_lists)
            else:
                num_lists = max(len(val), num_lists)

        # Construct array of kwargs dicts, each element of the list is a different
        # **kwargs dict.  each kwargs dict gives a different combination of
        # the possible values of the kwargs

        # initialize array
        list_of_kwarg_dicts = np.array([{} for x in range(num_lists)])

        # fill in array
        for i in np.arange(num_lists):
            list_of_kwarg_dicts[i] = {}
            for key in keywords.keys():
                # if it's a string, use it (there's only one)
                if isinstance(keywords[key], str):
                    list_of_kwarg_dicts[i][key] = keywords[key]
                # if there are more options, use the i'th val
                elif i < len(keywords[key]):
                    list_of_kwarg_dicts[i][key] = keywords[key][i]
                # if there are not more options, use the 0'th val
                else:
                    list_of_kwarg_dicts[i][key] = keywords[key][0]

    return list_of_kwarg_dicts


class TestOffAxisProjection(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()
        cls.curdir = os.getcwd()
        os.chdir(cls.tmpdir)

    @classmethod
    def tearDownClass(cls):
        os.chdir(cls.curdir)
        shutil.rmtree(cls.tmpdir)

    def test_oap(self):
        """Tests functionality of off_axis_projection and write_projection."""

        # args for off_axis_projection
        test_ds = fake_random_ds(64)
        c = test_ds.domain_center
        norm = [0.5, 0.5, 0.5]
        W = test_ds.arr([0.5, 0.5, 1.0], "unitary")
        N = 256
        field = ("gas", "density")
        oap_args = [test_ds, c, norm, W, N, field]

        # kwargs for off_axis_projection
        oap_kwargs = {}
        oap_kwargs["weight"] = (None, "cell_mass")
        oap_kwargs["no_ghost"] = (True, False)
        oap_kwargs["interpolated"] = (False,)
        oap_kwargs["north_vector"] = ((1, 0, 0), (0, 0.5, 1.0))
        oap_kwargs_list = expand_keywords(oap_kwargs)

        # args or write_projection
        fn = "test_%d.png"

        # kwargs for write_projection
        wp_kwargs = {}
        wp_kwargs["colorbar"] = (True, False)
        wp_kwargs["colorbar_label"] = "test"
        wp_kwargs["title"] = "test"
        wp_kwargs["vmin"] = (None,)
        wp_kwargs["vmax"] = (1e3, 1e5)
        wp_kwargs["take_log"] = (True, False)
        wp_kwargs["figsize"] = ((8, 6), [1, 1])
        wp_kwargs["dpi"] = (100, 50)
        wp_kwargs["cmap_name"] = ("cmyt.arbre", "cmyt.kelp")
        wp_kwargs_list = expand_keywords(wp_kwargs)

        # test all off_axis_projection kwargs and write_projection kwargs
        # make sure they are able to be projected, then remove and try next
        # iteration
        for i, oap_kwargs in enumerate(oap_kwargs_list):
            image = off_axis_projection(*oap_args, **oap_kwargs)
            for wp_kwargs in wp_kwargs_list:
                write_projection(image, fn % i, **wp_kwargs)
                assert_fname(fn % i)

        # Test remaining parameters of write_projection
        write_projection(image, "test_2", xlabel="x-axis", ylabel="y-axis")
        assert_fname("test_2.png")

        write_projection(image, "test_3.pdf", xlabel="x-axis", ylabel="y-axis")
        assert_fname("test_3.pdf")

        write_projection(image, "test_4.eps", xlabel="x-axis", ylabel="y-axis")
        assert_fname("test_4.eps")


def test_field_cut_off_axis_octree():
    ds = fake_octree_ds()
    cut = ds.all_data().cut_region('obj["gas", "density"]>0.5')
    p1 = OffAxisProjectionPlot(ds, [1, 0, 0], ("gas", "density"))
    p2 = OffAxisProjectionPlot(ds, [1, 0, 0], ("gas", "density"), data_source=cut)
    assert_equal(p2.frb["gas", "density"].min() == 0.0, True)  # Lots of zeros
    assert_equal((p1.frb["gas", "density"] == p2.frb["gas", "density"]).all(), False)
    p3 = OffAxisSlicePlot(ds, [1, 0, 0], ("gas", "density"))
    p4 = OffAxisSlicePlot(ds, [1, 0, 0], ("gas", "density"), data_source=cut)
    assert_equal((p3.frb["gas", "density"] == p4.frb["gas", "density"]).all(), False)
    p4rho = p4.frb["gas", "density"]
    assert_equal(np.nanmin(p4rho[p4rho > 0.0]) >= 0.5, True)


def test_offaxis_moment():
    ds = fake_random_ds(64)

    def _vlos_sq(field, data):
        return data["gas", "velocity_los"] ** 2

    ds.add_field(
        ("gas", "velocity_los_squared"),
        _vlos_sq,
        sampling_type="local",
        units="cm**2/s**2",
    )
    p1 = OffAxisProjectionPlot(
        ds,
        [1, 1, 1],
        [("gas", "velocity_los"), ("gas", "velocity_los_squared")],
        weight_field=("gas", "density"),
        moment=1,
        buff_size=(400, 400),
    )
    p2 = OffAxisProjectionPlot(
        ds,
        [1, 1, 1],
        ("gas", "velocity_los"),
        weight_field=("gas", "density"),
        moment=2,
        buff_size=(400, 400),
    )
    ## this failed because some <v**2> - <v>**2 values come out
    ## marginally < 0, resulting in unmatched NaN values in the
    ## first assert_rel_equal argument. The compute_stddev_image
    ## function used in OffAxisProjectionPlot checks for and deals
    ## with these cases. 
    #assert_rel_equal(
    #    np.sqrt(
    #        p1.frb["gas", "velocity_los_squared"] - p1.frb["gas", "velocity_los"] ** 2
    #    ),
    #    p2.frb["gas", "velocity_los"],
    #    10,
    #)
    p1_expsq = p1.frb["gas", "velocity_los_squared"]
    p1_sqexp = p1.frb["gas", "velocity_los"] ** 2
    p1res = np.sqrt(p1_expsq - p1_sqexp)
    # set values to zero that have <v>**2 - <v>**2 < 0, but
    # the absolute values are much smaller than the smallest 
    # postive values of <v>**2 and <v>**2 
    # (i.e., the difference is pretty much zero)
    mindiff = 1e-10 * min(np.min(p1_expsq[p1_expsq > 0]), 
                          np.min(p1_sqexp[p1_sqexp > 0]))
    print(mindiff)
    setzero = np.logical_and(p1_expsq - p1_sqexp < 0,
                             p1_expsq - p1_sqexp > -1. * mindiff)
    p1res[setzero] = 0.
    p2res = p2.frb["gas", "velocity_los"]
    assert_rel_equal(p1res, p2res, 10)


@pytest.mark.parametrize("weighted", [True, False])
@pytest.mark.parametrize("periodic", [True, False])
@pytest.mark.parametrize("depth", [None, (1., "cm"), (0.5, "cm")])
@pytest.mark.parametrize("shiftcenter", [False, True])
@pytest.mark.parametrize("northvector", [None, (1.e-4, 1., 0.)])
def test_sph_proj_general_offaxis(
        northvector: tuple[float, float, float] | None,
        shiftcenter: bool,
        depth: tuple[float, str] | None,
        periodic: bool, 
        weighted: bool) -> None:
    '''
    Same as the on-axis projections, but we rotate the basis vectors
    to test whether roations are handled ok. the rotation is chosen to
    be small so that in/exclusion of particles within bboxes, etc. 
    works out the same way. 
    We just send lines of sight through pixel centers for convenience.
    Particles at [0.5, 1.5, 2.5] (in each coordinate)
    smoothing lengths 0.25
    all particles have mass 1., density 1.5, 
    except the single center particle, with mass 2., density 3.
    
    Parameters:
    -----------
    northvector: tuple
        y-axis direction in the final plot (direction vector)
    shiftcenter: bool
        shift the coordinates to center the projection on. 
        (The grid is offset to this same center)
    depth: float or None
        depth of the projection slice
    periodic: bool
        assume periodic boundary conditions, or not
    weighted: bool
        make a weighted projection (density-weighted density), or not
    
    Returns:
    --------
    None
    '''
    if shiftcenter:
        center = np.array((0.625, 0.625, 0.625)) # cm
    else:
        center = np.array((1.5, 1.5, 1.5)) # cm
    bbox = unyt.unyt_array(np.array([[0., 3.], [0., 3.], [0., 3.]]), 'cm')
    hsml_factor = 0.5
    unitrho = 1.5
    # test correct centering, particle selection
    def makemasses(i, j, k):
        if i == j == k == 1:
            return 2. 
        else:
            return 1.
    #dl_shift = _integrate_kernel(_cubicspline_python, 
    #                             np.sqrt(2) * 0.125, 0.25)

    # result shouldn't depend explicitly on the center if we re-center
    # the data, unless we get cut-offs in the non-periodic case 
    # *almost* the z-axis
    # try to make sure dl differences from periodic wrapping are small
    epsilon = 1e-4 
    projaxis = np.array([epsilon, 0.00, np.sqrt(1. - epsilon**2)])
    e1dir = projaxis / np.sqrt(np.sum(projaxis**2))
    # TODO: figure out other (default) axes for basis vectors here
    if northvector is None:
        e2dir = np.array([0., 1., 0.])
    else:
        e2dir = np.asarray(northvector)
    e2dir = e2dir - np.sum(e1dir * e2dir) * e2dir # orthonormalize
    e2dir /= np.sqrt(np.sum(e2dir**2))
    e3dir = np.cross(e1dir, e2dir)

    ds = fake_sph_flexible_grid_ds(hsml_factor=hsml_factor, nperside=3,
                                   periodic=periodic, 
                                   offsets=np.full(3, 0.5),
                                   massgenerator=makemasses,
                                   unitrho=unitrho,
                                   bbox=bbox.v, 
                                   recenter=center,
                                   e1hat=e1dir, e2hat=e2dir, e3hat=e3dir)

    source = ds.all_data()
    # couple to dataset -> right unit registry
    center = ds.arr(center, 'cm')
    #print('position:\n', source['gas','position'])

    # m / rho, factor 1. / hsml**2 is included in the kernel integral
    # (density is adjusted, so same for center particle)
    prefactor = 1. / unitrho #/ (0.5 * 0.5)**2
    dl_cen = integrate_kernel(cubicspline_python, 0., 0.25)

    if weighted: 
        toweight_field = ("gas", "density")
    else:
        toweight_field = None
    # we don't actually want a plot, it's just a straightforward, 
    # common way to get an frb / image array
    prj = ProjectionPlot(ds, projaxis, ("gas", "density"),
                         width=(2.5, "cm"),
                         weight_field=toweight_field,
                         buff_size=(5, 5),
                         center=center,
                         data_source=source,
                         north_vector=northvector,
                         depth=depth)
    img = prj.frb.data[('gas', 'density')]
    if weighted:
        # periodic shifts will modify the (relative) dl values a bit
        expected_out = np.zeros((5, 5,), dtype=img.v.dtype)
        expected_out[::2, ::2] = unitrho
        if depth is None: 
            expected_out[2, 2] *= 1.5
        else:
            # only 2 * unitrho element included
            expected_out[2, 2] *= 2.
    else:
        expected_out = np.zeros((5, 5,), dtype=img.v.dtype)
        expected_out[::2, ::2] = dl_cen * prefactor * unitrho
        if depth is None:
            # 3 particles per l.o.s., including the denser one
            expected_out *= 3.
            expected_out[2, 2] *= 4. / 3.
        else:
            # 1 particle per l.o.s., including the denser one
            expected_out[2, 2] *= 2.
    # grid is shifted to the left -> 'missing' stuff at the left
    if (not periodic) and shiftcenter: 
        expected_out[:1, :] = 0.
        expected_out[:, :1] = 0.
    #print(axis, shiftcenter, depth, periodic, weighted)
    print('expected:\n', expected_out)
    print('recovered:\n', img.v)
    assert_rel_equal(expected_out, img.v, 4)
