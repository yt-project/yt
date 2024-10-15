import numpy as np
from numpy.testing import assert_equal

from yt import load
from yt.testing import (
    assert_rel_equal,
    cubicspline_python,
    fake_random_ds,
    fake_sph_grid_ds,
    integrate_kernel,
    requires_file,
    requires_module,
)
from yt.units._numpy_wrapper_functions import uconcatenate


def test_ray():
    for nproc in [1, 2, 4, 8]:
        ds = fake_random_ds(64, nprocs=nproc)
        dx = (ds.domain_right_edge - ds.domain_left_edge) / ds.domain_dimensions
        # Three we choose, to get varying vectors, and ten random
        pp1 = np.random.random((3, 13))
        pp2 = np.random.random((3, 13))
        pp1[:, 0] = [0.1, 0.2, 0.3]
        pp2[:, 0] = [0.8, 0.1, 0.4]
        pp1[:, 1] = [0.9, 0.2, 0.3]
        pp2[:, 1] = [0.8, 0.1, 0.4]
        pp1[:, 2] = [0.9, 0.2, 0.9]
        pp2[:, 2] = [0.8, 0.1, 0.4]
        unitary = ds.arr(1.0, "")
        for i in range(pp1.shape[1]):
            p1 = ds.arr(pp1[:, i] + 1e-8 * np.random.random(3), "code_length")
            p2 = ds.arr(pp2[:, i] + 1e-8 * np.random.random(3), "code_length")

            my_ray = ds.ray(p1, p2)
            assert_rel_equal(my_ray["dts"].sum(), unitary, 14)
            ray_cells = my_ray["dts"] > 0

            # find cells intersected by the ray
            my_all = ds.all_data()

            dt = np.abs(dx / (p2 - p1))
            tin = uconcatenate(
                [
                    [(my_all["index", "x"] - p1[0]) / (p2 - p1)[0] - 0.5 * dt[0]],
                    [(my_all["index", "y"] - p1[1]) / (p2 - p1)[1] - 0.5 * dt[1]],
                    [(my_all["index", "z"] - p1[2]) / (p2 - p1)[2] - 0.5 * dt[2]],
                ]
            )
            tout = uconcatenate(
                [
                    [(my_all["index", "x"] - p1[0]) / (p2 - p1)[0] + 0.5 * dt[0]],
                    [(my_all["index", "y"] - p1[1]) / (p2 - p1)[1] + 0.5 * dt[1]],
                    [(my_all["index", "z"] - p1[2]) / (p2 - p1)[2] + 0.5 * dt[2]],
                ]
            )
            tin = tin.max(axis=0)
            tout = tout.min(axis=0)
            my_cells = (tin < tout) & (tin < 1) & (tout > 0)

            assert_equal(ray_cells.sum(), my_cells.sum())
            assert_rel_equal(
                my_ray["gas", "density"][ray_cells].sum(),
                my_all["gas", "density"][my_cells].sum(),
                14,
            )
            assert_rel_equal(my_ray["dts"].sum(), unitary, 14)


@requires_module("h5py")
@requires_file("GadgetDiskGalaxy/snapshot_200.hdf5")
def test_ray_particle():
    ds = load("GadgetDiskGalaxy/snapshot_200.hdf5")
    ray = ds.ray(ds.domain_left_edge, ds.domain_right_edge)
    assert_equal(ray["t"].shape, (1451,))
    assert ray["dts"].sum(dtype="f8") > 0


## test that kernels integrate correctly
# (1) including the right particles
# (2) correct t and dts values for those particles
## fake_sph_grid_ds:
# This dataset should have 27 particles with the particles arranged
# uniformly on a 3D grid. The bottom left corner is (0.5,0.5,0.5) and
# the top right corner is (2.5,2.5,2.5). All particles will have
# non-overlapping smoothing regions with a radius of 0.05, masses of
# 1, and densities of 1, and zero velocity.


def test_ray_particle2():
    kernelfunc = cubicspline_python
    ds = fake_sph_grid_ds(hsml_factor=1.0)
    ds.kernel_name = "cubic"

    ## Ray through the one particle at (0.5, 0.5, 0.5):
    ## test basic kernel integration
    eps = 0.0  # 1e-7
    start0 = np.array((1.0 + eps, 0.0, 0.5))
    end0 = np.array((0.0, 1.0 + eps, 0.5))
    ray0 = ds.ray(start0, end0)
    b0 = np.array([np.sqrt(2.0) * eps])
    hsml0 = np.array([0.05])
    len0 = np.sqrt(np.sum((end0 - start0) ** 2))
    # for a ParticleDataset like this one, the Ray object attempts
    # to generate the 't' and 'dts' fields using the grid method
    ray0.field_data["t"] = ray0.ds.arr(ray0._generate_container_field_sph("t"))
    ray0.field_data["dts"] = ray0.ds.arr(ray0._generate_container_field_sph("dts"))
    # not demanding too much precision;
    # from kernel volume integrals, the linear interpolation
    # restricts you to 4 -- 5 digits precision
    assert_equal(ray0["t"].shape, (1,))
    assert_rel_equal(ray0["t"], np.array([0.5]), 5)
    assert_rel_equal(ray0[("gas", "position")].v, np.array([[0.5, 0.5, 0.5]]), 5)
    dl0 = integrate_kernel(kernelfunc, b0, hsml0)
    dl0 *= ray0[("gas", "mass")].v / ray0[("gas", "density")].v
    assert_rel_equal(ray0[("dts")].v, dl0 / len0, 4)

    ## Ray in the middle of the box:
    ## test end points, >1 particle
    start1 = np.array((1.53, 0.53, 1.0))
    end1 = np.array((1.53, 0.53, 3.0))
    ray1 = ds.ray(start1, end1)
    b1 = np.array([np.sqrt(2.0) * 0.03] * 2)
    hsml1 = np.array([0.05] * 2)
    len1 = np.sqrt(np.sum((end1 - start1) ** 2))
    # for a ParticleDataset like this one, the Ray object attempts
    # to generate the 't' and 'dts' fields using the grid method
    ray1.field_data["t"] = ray1.ds.arr(ray1._generate_container_field_sph("t"))
    ray1.field_data["dts"] = ray1.ds.arr(ray1._generate_container_field_sph("dts"))
    # not demanding too much precision;
    # from kernel volume integrals, the linear interpolation
    # restricts you to 4 -- 5 digits precision
    assert_equal(ray1["t"].shape, (2,))
    assert_rel_equal(ray1["t"], np.array([0.25, 0.75]), 5)
    assert_rel_equal(
        ray1[("gas", "position")].v, np.array([[1.5, 0.5, 1.5], [1.5, 0.5, 2.5]]), 5
    )
    dl1 = integrate_kernel(kernelfunc, b1, hsml1)
    dl1 *= ray1[("gas", "mass")].v / ray1[("gas", "density")].v
    assert_rel_equal(ray1[("dts")].v, dl1 / len1, 4)

    ## Ray missing all particles:
    ## test handling of size-0 selections
    start2 = np.array((1.0, 2.0, 0.0))
    end2 = np.array((1.0, 2.0, 3.0))
    ray2 = ds.ray(start2, end2)
    # for a ParticleDataset like this one, the Ray object attempts
    # to generate the 't' and 'dts' fields using the grid method
    ray2.field_data["t"] = ray2.ds.arr(ray2._generate_container_field_sph("t"))
    ray2.field_data["dts"] = ray2.ds.arr(ray2._generate_container_field_sph("dts"))
    assert_equal(ray2["t"].shape, (0,))
    assert_equal(ray2["dts"].shape, (0,))
    assert_equal(ray2[("gas", "position")].v.shape, (0, 3))
