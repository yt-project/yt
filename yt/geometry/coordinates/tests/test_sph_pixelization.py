import numpy as np

import yt
from yt.testing import (
    assert_rel_equal,
    cubicspline_python,
    fake_sph_flexible_grid_ds,
    integrate_kernel,
    requires_file,
)
from yt.utilities.math_utils import compute_stddev_image

## off-axis projection tests for SPH data are in
## yt/visualization/tests/test_offaxisprojection.py


magneticum = "MagneticumCluster/snap_132"

mag_kwargs = {
    "long_ids": True,
    "field_spec": "magneticum_box2_hr",
}


@requires_file(magneticum)
def test_sph_moment():
    ds = yt.load(magneticum, **mag_kwargs)

    def _vysq(field, data):
        return data["gas", "velocity_y"] ** 2

    ds.add_field(("gas", "vysq"), _vysq, sampling_type="local", units="cm**2/s**2")
    prj1 = yt.ProjectionPlot(
        ds,
        "y",
        [("gas", "velocity_y"), ("gas", "vysq")],
        weight_field=("gas", "density"),
        moment=1,
        buff_size=(400, 400),
    )
    prj2 = yt.ProjectionPlot(
        ds,
        "y",
        ("gas", "velocity_y"),
        moment=2,
        weight_field=("gas", "density"),
        buff_size=(400, 400),
    )
    sigy = compute_stddev_image(prj1.frb["gas", "vysq"], prj1.frb["gas", "velocity_y"])
    assert_rel_equal(sigy, prj2.frb["gas", "velocity_y"].d, 10)


def test_sph_projection_basic1():
    """
    small, uniform grid: expected values for given dl?
    pixel centers at 0.5, 1., 1.5, 2., 2.5
    particles at 0.5, 1.5, 2.5
    """
    bbox = np.array([[0.0, 3.0]] * 3)
    ds = fake_sph_flexible_grid_ds(hsml_factor=1.0, nperside=3, bbox=bbox)
    # works, but no depth control (at least without specific filters)
    proj = ds.proj(("gas", "density"), 2)
    frb = proj.to_frb(
        width=(2.5, "cm"),
        resolution=(5, 5),
        height=(2.5, "cm"),
        center=np.array([1.5, 1.5, 1.5]),
        periodic=False,
    )
    out = frb.get_image(("gas", "density"))

    expected_out = np.zeros((5, 5), dtype=np.float64)
    dl_1part = integrate_kernel(cubicspline_python, 0.0, 0.5)
    linedens_1part = dl_1part * 1.0  # unit mass, density
    linedens = 3.0 * linedens_1part
    expected_out[::2, ::2] = linedens

    assert_rel_equal(expected_out, out.v, 5)
    # return out


def test_sph_projection_basic2():
    """
    small, uniform grid: expected values for given dl?
    pixel centers at 0.5, 1., 1.5, 2., 2.5
    particles at 0.5, 1.5, 2.5
    but hsml radii are 0.25 -> try non-zero impact parameters,
    other pixels are still zero.
    """
    bbox = np.array([[0.0, 3.0]] * 3)
    ds = fake_sph_flexible_grid_ds(hsml_factor=0.5, nperside=3, bbox=bbox)
    proj = ds.proj(("gas", "density"), 2)
    frb = proj.to_frb(
        width=(2.5, "cm"),
        resolution=(5, 5),
        height=(2.5, "cm"),
        center=np.array([1.375, 1.375, 1.5]),
        periodic=False,
    )
    out = frb.get_image(("gas", "density"))

    expected_out = np.zeros((5, 5), dtype=np.float64)
    dl_1part = integrate_kernel(cubicspline_python, np.sqrt(2) * 0.125, 0.25)
    linedens_1part = dl_1part * 1.0  # unit mass, density
    linedens = 3.0 * linedens_1part
    expected_out[::2, ::2] = linedens

    # print(expected_out)
    # print(out.v)
    assert_rel_equal(expected_out, out.v, 4)
    # return out


def get_dataset_sphrefine(reflevel: int = 1):
    """
    constant density particle grid,
    with increasing particle sampling
    """
    lenfact = (1.0 / 3.0) ** (reflevel - 1)
    massfact = lenfact**3
    nperside = 3**reflevel

    e1hat = np.array([lenfact, 0, 0])
    e2hat = np.array([0, lenfact, 0])
    e3hat = np.array([0, 0, lenfact])
    hsml_factor = lenfact
    bbox = np.array([[0.0, 3.0]] * 3)
    offsets = np.ones(3, dtype=np.float64) * 0.5  # in units of ehat

    def refmass(i: int, j: int, k: int) -> float:
        return massfact

    unitrho = 1.0 / massfact  # want density 1 for decreasing mass

    ds = fake_sph_flexible_grid_ds(
        hsml_factor=hsml_factor,
        nperside=nperside,
        periodic=True,
        e1hat=e1hat,
        e2hat=e2hat,
        e3hat=e3hat,
        offsets=offsets,
        massgenerator=refmass,
        unitrho=unitrho,
        bbox=bbox,
    )
    return ds


def getdata_test_gridproj2():
    # initial pixel centers at 0.5, 1., 1.5, 2., 2.5
    # particles at 0.5, 1.5, 2.5
    # refine particle grid, check if pixel values remain the
    # same in the pixels passing through initial particle centers
    outlist = []
    dss = []
    for rl in range(1, 4):
        ds = get_dataset_sphrefine(reflevel=rl)
        proj = ds.proj(("gas", "density"), 2)
        frb = proj.to_frb(
            width=(2.5, "cm"),
            resolution=(5, 5),
            height=(2.5, "cm"),
            center=np.array([1.5, 1.5, 1.5]),
            periodic=False,
        )
        out = frb.get_image(("gas", "density"))
        outlist.append(out)
        dss.append(ds)
    return outlist, dss


def test_sph_gridproj_reseffect1():
    """
    Comparing same pixel centers with higher particle resolution.
    The pixel centers are at x/y coordinates [0.5, 1., 1.5, 2., 2.5]
    at the first level, the spacing halves at each level.
    Checking the pixels at [0.5, 1.5, 2.5],
    which should have the same values at each resolution.
    """
    imgs, _ = getdata_test_gridproj2()
    ref = imgs[-1]
    for img in imgs:
        assert_rel_equal(
            img[:: img.shape[0] // 2, :: img.shape[1] // 2],
            ref[:: ref.shape[0] // 2, :: ref.shape[1] // 2],
            4,
        )


def test_sph_gridproj_reseffect2():
    """
    refine the pixel grid instead of the particle grid
    """
    ds = get_dataset_sphrefine(reflevel=2)
    proj = ds.proj(("gas", "density"), 2)
    imgs = {}
    maxrl = 5
    for rl in range(1, maxrl + 1):
        npix = 1 + 2 ** (rl + 1)
        margin = 0.5 - 0.5 ** (rl + 1)
        frb = proj.to_frb(
            width=(3.0 - 2.0 * margin, "cm"),
            resolution=(npix, npix),
            height=(3.0 - 2.0 * margin, "cm"),
            center=np.array([1.5, 1.5, 1.5]),
            periodic=False,
        )
        out = frb.get_image(("gas", "density"))
        imgs[rl] = out
    ref = imgs[maxrl]
    pixspace_ref = 2 ** (maxrl)
    for rl in imgs:
        img = imgs[rl]
        pixspace = 2 ** (rl)
        # print(f'Grid refinement level {rl}:')
        assert_rel_equal(
            img[::pixspace, ::pixspace], ref[::pixspace_ref, ::pixspace_ref], 4
        )
