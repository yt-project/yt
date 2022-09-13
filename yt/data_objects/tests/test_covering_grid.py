import numpy as np

from yt.fields.derived_field import ValidateParameter
from yt.loaders import load, load_particles
from yt.testing import (
    assert_almost_equal,
    assert_array_equal,
    assert_equal,
    fake_octree_ds,
    fake_random_ds,
    requires_file,
    requires_module,
)
from yt.units import kpc

# cylindrical data for covering_grid test
cyl_2d = "WDMerger_hdf5_chk_1000/WDMerger_hdf5_chk_1000.hdf5"
cyl_3d = "MHD_Cyl3d_hdf5_plt_cnt_0100/MHD_Cyl3d_hdf5_plt_cnt_0100.hdf5"


def setup():
    from yt.config import ytcfg

    ytcfg["yt", "internals", "within_testing"] = True


@requires_file(cyl_2d)
@requires_file(cyl_3d)
def test_covering_grid():
    # We decompose in different ways
    for level in [0, 1, 2]:
        for nprocs in [1, 2, 4, 8]:
            ds = fake_random_ds(16, nprocs=nprocs)
            axis_name = ds.coordinates.axis_name
            dn = ds.refine_by**level
            cg = ds.covering_grid(level, [0.0, 0.0, 0.0], dn * ds.domain_dimensions)
            # Test coordinate generation
            assert_equal(np.unique(cg[("index", f"d{axis_name[0]}")]).size, 1)
            xmi = cg[("index", axis_name[0])].min()
            xma = cg[("index", axis_name[0])].max()
            dx = cg[("index", f"d{axis_name[0]}")].flat[0:1]
            edges = ds.arr([[0, 1], [0, 1], [0, 1]], "code_length")
            assert_equal(xmi, edges[0, 0] + dx / 2.0)
            assert_equal(xmi, cg[("index", axis_name[0])][0, 0, 0])
            assert_equal(xmi, cg[("index", axis_name[0])][0, 1, 1])
            assert_equal(xma, edges[0, 1] - dx / 2.0)
            assert_equal(xma, cg[("index", axis_name[0])][-1, 0, 0])
            assert_equal(xma, cg[("index", axis_name[0])][-1, 1, 1])
            assert_equal(np.unique(cg[("index", f"d{axis_name[1]}")]).size, 1)
            ymi = cg[("index", axis_name[1])].min()
            yma = cg[("index", axis_name[1])].max()
            dy = cg[("index", f"d{axis_name[1]}")][0]
            assert_equal(ymi, edges[1, 0] + dy / 2.0)
            assert_equal(ymi, cg[("index", axis_name[1])][0, 0, 0])
            assert_equal(ymi, cg[("index", axis_name[1])][1, 0, 1])
            assert_equal(yma, edges[1, 1] - dy / 2.0)
            assert_equal(yma, cg[("index", axis_name[1])][0, -1, 0])
            assert_equal(yma, cg[("index", axis_name[1])][1, -1, 1])
            assert_equal(np.unique(cg[("index", f"d{axis_name[2]}")]).size, 1)
            zmi = cg[("index", axis_name[2])].min()
            zma = cg[("index", axis_name[2])].max()
            dz = cg[("index", f"d{axis_name[2]}")][0]
            assert_equal(zmi, edges[2, 0] + dz / 2.0)
            assert_equal(zmi, cg[("index", axis_name[2])][0, 0, 0])
            assert_equal(zmi, cg[("index", axis_name[2])][1, 1, 0])
            assert_equal(zma, edges[2, 1] - dz / 2.0)
            assert_equal(zma, cg[("index", axis_name[2])][0, 0, -1])
            assert_equal(zma, cg[("index", axis_name[2])][1, 1, -1])
            # Now we test other attributes
            assert_equal(cg[("index", "ones")].max(), 1.0)
            assert_equal(cg[("index", "ones")].min(), 1.0)
            assert_equal(cg[("index", "grid_level")], level)
            assert_equal(cg[("index", "cell_volume")].sum(), ds.domain_width.prod())
            for g in ds.index.grids:
                di = g.get_global_startindex()
                dd = g.ActiveDimensions
                for i in range(dn):
                    f = cg[("gas", "density")][
                        dn * di[0] + i : dn * (di[0] + dd[0]) + i : dn,
                        dn * di[1] + i : dn * (di[1] + dd[1]) + i : dn,
                        dn * di[2] + i : dn * (di[2] + dd[2]) + i : dn,
                    ]
                    assert_equal(f, g[("gas", "density")])

    # More tests for cylindrical geometry
    for fn in [cyl_2d, cyl_3d]:
        ds = load(fn)
        ad = ds.all_data()
        upper_ad = ad.cut_region(["obj[('index', 'z')] > 0"])
        sp = ds.sphere((0, 0, 0), 0.5 * ds.domain_width[0], data_source=upper_ad)
        sp.quantities.total_mass()


@requires_file(cyl_2d)
@requires_file(cyl_3d)
def test_covering_grid_data_source():
    # test the data_source kwarg (new with PR 4063)

    for level in [0, 1, 2]:
        for nprocs in [1, 2, 4, 8]:
            ds = fake_random_ds(16, nprocs=nprocs)
            dn = ds.refine_by**level
            cg = ds.covering_grid(level, [0.0, 0.0, 0.0], dn * ds.domain_dimensions)

            # check that the regularized sphere is centered where it should be,
            # to within the tolerance of the grid. Use a center offsect from
            # the domain center by a bit
            center = ds.domain_center + np.min(ds.domain_width) * 0.15
            sp = ds.sphere(center, (0.1, "code_length"))
            dims = dn * ds.domain_dimensions
            cg_sp = ds.covering_grid(level, [0.0, 0.0, 0.0], dims, data_source=sp)

            # find the discrete center of the sphere by averaging the position
            # along each dimension where values are non-zero
            cg_mask = cg_sp["gas", "density"] != 0
            discrete_c = ds.arr(
                [
                    np.mean(cg_sp[("index", dim)][cg_mask])
                    for dim in ds.coordinates.axis_order
                ]
            )

            # should be centered to within the tolerance of the grid at least
            grid_tol = ds.arr(
                [cg_sp[("index", dim)][0, 0, 0] for dim in ds.coordinates.axis_order]
            )
            assert np.all(np.abs(discrete_c - center) <= grid_tol)

            # check that a region covering the whole domain matches no data_source
            reg = ds.region(ds.domain_center, ds.domain_left_edge, ds.domain_right_edge)
            cg_reg = ds.covering_grid(level, [0.0, 0.0, 0.0], dims, data_source=reg)
            assert np.all(cg[("gas", "density")] == cg_reg[("gas", "density")])

            # check that a box covering a subset of the domain is the right volume
            right_edge = ds.domain_left_edge + ds.domain_width * 0.5
            c = (right_edge + ds.domain_left_edge) / 2
            reg = ds.region(c, ds.domain_left_edge, right_edge)
            cg_reg = ds.covering_grid(level, [0.0, 0.0, 0.0], dims, data_source=reg)
            box_vol = cg_reg[("index", "cell_volume")][
                cg_reg["gas", "density"] != 0
            ].sum()
            actual_vol = np.prod(right_edge - ds.domain_left_edge)
            assert box_vol == actual_vol


@requires_module("xarray")
def test_xarray_export():
    def _run_tests(cg):
        xarr = cg.to_xarray(fields=[("gas", "density"), ("gas", "temperature")])
        assert ("gas", "density") in xarr.variables
        assert ("gas", "temperature") in xarr.variables
        assert ("gas", "specific_thermal_energy") not in xarr.variables
        assert "x" in xarr.coords
        assert "y" in xarr.coords
        assert "z" in xarr.coords
        assert xarr.dims["x"] == dn * ds.domain_dimensions[0]
        assert xarr.dims["y"] == dn * ds.domain_dimensions[1]
        assert xarr.dims["z"] == dn * ds.domain_dimensions[2]
        assert_equal(xarr.x, cg[("index", "x")][:, 0, 0])
        assert_equal(xarr.y, cg[("index", "y")][0, :, 0])
        assert_equal(xarr.z, cg[("index", "z")][0, 0, :])

    fields = ("density", "temperature", "specific_thermal_energy")
    units = ("g/cm**3", "K", "erg/g")
    for level in [0, 1, 2]:
        ds = fake_random_ds(16, fields=fields, units=units)
        dn = ds.refine_by**level
        rcg = ds.covering_grid(level, [0.0, 0.0, 0.0], dn * ds.domain_dimensions)
        _run_tests(rcg)
        scg = ds.smoothed_covering_grid(
            level, [0.0, 0.0, 0.0], dn * ds.domain_dimensions
        )
        _run_tests(scg)
        ag1 = ds.arbitrary_grid(
            [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], dn * ds.domain_dimensions
        )
        _run_tests(ag1)
        ag2 = ds.arbitrary_grid(
            [0.1, 0.3, 0.2], [0.4, 1.0, 0.9], dn * ds.domain_dimensions
        )
        _run_tests(ag2)


def test_smoothed_covering_grid():
    # We decompose in different ways
    for level in [0, 1, 2]:
        for nprocs in [1, 2, 4, 8]:
            ds = fake_random_ds(16, nprocs=nprocs)
            dn = ds.refine_by**level
            cg = ds.smoothed_covering_grid(
                level, [0.0, 0.0, 0.0], dn * ds.domain_dimensions
            )
            assert_equal(cg[("index", "ones")].max(), 1.0)
            assert_equal(cg[("index", "ones")].min(), 1.0)
            assert_equal(cg[("index", "cell_volume")].sum(), ds.domain_width.prod())
            for g in ds.index.grids:
                if level != g.Level:
                    continue
                di = g.get_global_startindex()
                dd = g.ActiveDimensions
                for i in range(dn):
                    f = cg[("gas", "density")][
                        dn * di[0] + i : dn * (di[0] + dd[0]) + i : dn,
                        dn * di[1] + i : dn * (di[1] + dd[1]) + i : dn,
                        dn * di[2] + i : dn * (di[2] + dd[2]) + i : dn,
                    ]
                    assert_equal(f, g[("gas", "density")])


def test_arbitrary_grid():
    for ncells in [32, 64]:
        for px in [0.125, 0.25, 0.55519]:

            particle_data = {
                "particle_position_x": np.array([px]),
                "particle_position_y": np.array([0.5]),
                "particle_position_z": np.array([0.5]),
                "particle_mass": np.array([1.0]),
            }

            ds = load_particles(particle_data)

            for dims in ([ncells] * 3, [ncells, ncells / 2, ncells / 4]):
                LE = np.array([0.05, 0.05, 0.05])
                RE = np.array([0.95, 0.95, 0.95])
                dims = np.array(dims)

                dds = (RE - LE) / dims
                volume = ds.quan(np.product(dds), "cm**3")

                obj = ds.arbitrary_grid(LE, RE, dims)
                deposited_mass = obj[("deposit", "all_density")].sum() * volume

                assert_equal(deposited_mass, ds.quan(1.0, "g"))

                LE = np.array([0.00, 0.00, 0.00])
                RE = np.array([0.05, 0.05, 0.05])

                obj = ds.arbitrary_grid(LE, RE, dims)

                deposited_mass = obj[("deposit", "all_density")].sum()

                assert_equal(deposited_mass, 0)

    # Test that we get identical results to the covering grid for unigrid data.
    # Testing AMR data is much harder.
    for nprocs in [1, 2, 4, 8]:
        ds = fake_random_ds(32, nprocs=nprocs)
        for ref_level in [0, 1, 2]:
            cg = ds.covering_grid(
                ref_level, [0.0, 0.0, 0.0], 2**ref_level * ds.domain_dimensions
            )
            ag = ds.arbitrary_grid(
                [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 2**ref_level * ds.domain_dimensions
            )
            assert_almost_equal(cg[("gas", "density")], ag[("gas", "density")])


def test_octree_cg():
    ds = fake_octree_ds(num_zones=1, partial_coverage=0)
    cgrid = ds.covering_grid(
        0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions
    )
    density_field = cgrid[("gas", "density")]
    assert_equal((density_field == 0.0).sum(), 0)


def test_smoothed_covering_grid_2d_dataset():
    ds = fake_random_ds([32, 32, 1], nprocs=4)
    ds.force_periodicity()
    scg = ds.smoothed_covering_grid(1, [0.0, 0.0, 0.0], [32, 32, 1])
    assert_equal(scg[("gas", "density")].shape, [32, 32, 1])


def test_arbitrary_grid_derived_field():
    def custom_metal_density(field, data):
        # Calculating some random value
        return data["gas", "density"] * np.random.random_sample()

    ds = fake_random_ds(64, nprocs=8, particles=16**2)
    ds.add_field(
        ("gas", "Metal_Density"),
        units="g/cm**3",
        function=custom_metal_density,
        sampling_type="cell",
    )

    def _tracerf(field, data):
        return data[("gas", "Metal_Density")] / data["gas", "density"]

    ds.add_field(
        ("gas", "tracerf"),
        function=_tracerf,
        units="dimensionless",
        sampling_type="cell",
        take_log=False,
    )

    galgas = ds.arbitrary_grid([0.4, 0.4, 0.4], [0.99, 0.99, 0.99], dims=[32, 32, 32])
    galgas[("gas", "tracerf")]


def test_arbitrary_field_parameters():
    def _test_field(field, data):
        par = data.get_field_parameter("test_parameter")
        return par * data["all", "particle_mass"]

    ds = fake_random_ds(64, nprocs=8, particles=16**2)
    ds.add_field(
        ("all", "test_field"),
        units="g",
        function=_test_field,
        sampling_type="particle",
        validators=[ValidateParameter("test_parameter")],
    )

    agrid = ds.arbitrary_grid([0.4, 0.4, 0.4], [0.99, 0.99, 0.99], dims=[32, 32, 32])
    agrid.set_field_parameter("test_parameter", 2)
    assert_array_equal(2 * agrid["all", "particle_mass"], agrid["all", "test_field"])


def test_arbitrary_grid_edge():
    # Tests bug fix for issue #2087
    # Regardless of how left_edge and right_edge are passed, the result should be
    # a YTArray with a unit registry that matches that of the dataset.
    dims = [32, 32, 32]
    ds = fake_random_ds(dims)
    # Test when edge is a list, numpy array, YTArray with dataset units, and
    # YTArray with non-dataset units
    ledge = [
        [0.0, 0.0, 0.0],
        np.array([0.0, 0.0, 0.0]),
        [0.0, 0.0, 0.0] * ds.length_unit,
        [0.0, 0.0, 0.0] * kpc,
    ]

    redge = [
        [1.0, 1.0, 1.0],
        np.array([1.0, 1.0, 1.0]),
        [1.0, 1.0, 1.0] * ds.length_unit,
        [1.0, 1.0, 1.0] * kpc,
    ]

    ledge_ans = [
        [0.0, 0.0, 0.0] * ds.length_unit.to("code_length"),
        np.array([0.0, 0.0, 0.0]) * ds.length_unit.to("code_length"),
        [0.0, 0.0, 0.0] * ds.length_unit,
        [0.0, 0.0, 0.0] * kpc,
    ]

    redge_ans = [
        [1.0, 1.0, 1.0] * ds.length_unit.to("code_length"),
        np.array([1.0, 1.0, 1.0]) * ds.length_unit.to("code_length"),
        [1.0, 1.0, 1.0] * ds.length_unit,
        [1.0, 1.0, 1.0] * kpc,
    ]

    for le, re, le_ans, re_ans in zip(ledge, redge, ledge_ans, redge_ans):
        ag = ds.arbitrary_grid(left_edge=le, right_edge=re, dims=dims)
        assert np.array_equal(ag.left_edge, le_ans)
        assert np.array_equal(ag.right_edge, re_ans)
        assert ag.left_edge.units.registry == ds.unit_registry
        assert ag.right_edge.units.registry == ds.unit_registry
        ag[("gas", "density")]
