import numpy as np
import pytest
import unyt
from numpy.testing import assert_allclose

from yt.testing import (
    assert_rel_equal,
    cubicspline_python,
    fake_amr_ds,
    fake_sph_flexible_grid_ds,
    fake_sph_grid_ds,
    integrate_kernel,
    requires_module_pytest,
)
from yt.visualization.api import OffAxisProjectionPlot, ProjectionPlot


@pytest.mark.parametrize("weighted", [True, False])
@pytest.mark.parametrize("periodic", [True, False])
@pytest.mark.parametrize("depth", [None, (1.0, "cm"), (0.5, "cm")])
@pytest.mark.parametrize("shiftcenter", [False, True])
@pytest.mark.parametrize("northvector", [None, (1.0e-4, 1.0, 0.0)])
def test_sph_proj_general_offaxis(
    northvector: tuple[float, float, float] | None,
    shiftcenter: bool,
    depth: tuple[float, str] | None,
    periodic: bool,
    weighted: bool,
) -> None:
    """
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
    """
    if shiftcenter:
        center = np.array((0.625, 0.625, 0.625))  # cm
    else:
        center = np.array((1.5, 1.5, 1.5))  # cm
    bbox = unyt.unyt_array(np.array([[0.0, 3.0], [0.0, 3.0], [0.0, 3.0]]), "cm")
    hsml_factor = 0.5
    unitrho = 1.5

    # test correct centering, particle selection
    def makemasses(i, j, k):
        if i == j == k == 1:
            return 2.0
        else:
            return 1.0

    # result shouldn't depend explicitly on the center if we re-center
    # the data, unless we get cut-offs in the non-periodic case
    # *almost* the z-axis
    # try to make sure dl differences from periodic wrapping are small
    epsilon = 1e-4
    projaxis = np.array([epsilon, 0.00, np.sqrt(1.0 - epsilon**2)])
    e1dir = projaxis / np.sqrt(np.sum(projaxis**2))
    # TODO: figure out other (default) axes for basis vectors here
    if northvector is None:
        e2dir = np.array([0.0, 1.0, 0.0])
    else:
        e2dir = np.asarray(northvector)
    e2dir = e2dir - np.sum(e1dir * e2dir) * e2dir  # orthonormalize
    e2dir /= np.sqrt(np.sum(e2dir**2))
    e3dir = np.cross(e1dir, e2dir)

    ds = fake_sph_flexible_grid_ds(
        hsml_factor=hsml_factor,
        nperside=3,
        periodic=periodic,
        offsets=np.full(3, 0.5),
        massgenerator=makemasses,
        unitrho=unitrho,
        bbox=bbox.v,
        recenter=center,
        e1hat=e1dir,
        e2hat=e2dir,
        e3hat=e3dir,
    )

    source = ds.all_data()
    # couple to dataset -> right unit registry
    center = ds.arr(center, "cm")
    # print('position:\n', source['gas','position'])

    # m / rho, factor 1. / hsml**2 is included in the kernel integral
    # (density is adjusted, so same for center particle)
    prefactor = 1.0 / unitrho  # / (0.5 * 0.5)**2
    dl_cen = integrate_kernel(cubicspline_python, 0.0, 0.25)

    if weighted:
        toweight_field = ("gas", "density")
    else:
        toweight_field = None
    # we don't actually want a plot, it's just a straightforward,
    # common way to get an frb / image array
    prj = ProjectionPlot(
        ds,
        projaxis,
        ("gas", "density"),
        width=(2.5, "cm"),
        weight_field=toweight_field,
        buff_size=(5, 5),
        center=center,
        data_source=source,
        north_vector=northvector,
        depth=depth,
    )
    img = prj.frb.data[("gas", "density")]
    if weighted:
        # periodic shifts will modify the (relative) dl values a bit
        expected_out = np.zeros(
            (
                5,
                5,
            ),
            dtype=img.v.dtype,
        )
        expected_out[::2, ::2] = unitrho
        if depth is None:
            expected_out[2, 2] *= 1.5
        else:
            # only 2 * unitrho element included
            expected_out[2, 2] *= 2.0
    else:
        expected_out = np.zeros(
            (
                5,
                5,
            ),
            dtype=img.v.dtype,
        )
        expected_out[::2, ::2] = dl_cen * prefactor * unitrho
        if depth is None:
            # 3 particles per l.o.s., including the denser one
            expected_out *= 3.0
            expected_out[2, 2] *= 4.0 / 3.0
        else:
            # 1 particle per l.o.s., including the denser one
            expected_out[2, 2] *= 2.0
    # grid is shifted to the left -> 'missing' stuff at the left
    if (not periodic) and shiftcenter:
        expected_out[:1, :] = 0.0
        expected_out[:, :1] = 0.0
    # print(axis, shiftcenter, depth, periodic, weighted)
    # print("expected:\n", expected_out)
    # print("recovered:\n", img.v)
    assert_rel_equal(expected_out, img.v, 4)


_diag_dist = np.sqrt(3.0)  # diagonal distance of a cube with length 1.
# each case is depth, center, expected integrated distance
_cases_to_test = [
    (_diag_dist / 3.0, "domain_left_edge", _diag_dist / 3.0 / 2.0),
    (_diag_dist * 2.0, "domain_left_edge", _diag_dist),
    (_diag_dist * 4.0, "domain_left_edge", _diag_dist),
    (None, "domain_center", _diag_dist),
]


@pytest.mark.parametrize("depth,proj_center,expected", _cases_to_test)
def test_offaxisprojection_depth(depth, proj_center, expected):
    # this checks that the depth keyword argument works as expected.
    # in all cases, it integrates the (index, ones) field for a normal
    # pointing to the right edge corner of the domain.
    #
    # For the tests where the projection is centered on the left edge,
    # the integrate distance will scale as depth / 2.0. When centered
    # on the origin, it will scale with depth. The integrated distance
    # should max out at the diagonal distance of the cube (when the depth
    # exceeds the cube diagonal distance).
    #
    # Also note that the accuracy will depend on the buffer dimensions:
    # using the default (800,800) results in accuracy of about 1 percent

    ds = fake_amr_ds()

    n = [1.0, 1.0, 1.0]
    c = getattr(ds, proj_center)
    field = ("index", "ones")

    p = ProjectionPlot(ds, n, field, depth=depth, weight_field=None, center=c)

    maxval = p.frb[field].max().d
    assert_allclose(expected, maxval, atol=1e-2)


_sph_test_cases = [
    ([1.0, 1.0, 0.7], 27),
    ([1.0, 1.0, 1], 19),
    ([0.0, 0.0, 1], 9),
]


@requires_module_pytest("contourpy")
@pytest.mark.parametrize("normal_vec, n_particles", _sph_test_cases)
def test_offaxisprojection_sph_defaultdepth(normal_vec, n_particles):
    # checks that particles are picked up as expected for a range of
    # depths and normal vectors. Certain viewing angles will result
    # in overlapping particles (since fake_sph_grid_ds aligns particles
    # on a grid): n_particles is the expected number of circles in the
    # resulting image for the given normal vector. Circle counts are
    # calculated here using a contour generator.
    from contourpy import contour_generator

    ds = fake_sph_grid_ds()
    c = ds.domain_center
    diag_dist = np.linalg.norm(ds.domain_width)
    field = ("gas", "mass")
    p = OffAxisProjectionPlot(
        ds, normal_vec, field, weight_field=None, center=c, width=diag_dist
    )
    p.render()

    # get the number of circles in the plot
    cg = contour_generator(z=p.frb[("gas", "mass")].d)
    assert n_particles == len(cg.lines(1.0))
