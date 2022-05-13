from numpy.testing import assert_allclose, assert_array_less, assert_raises

import yt
from yt.loaders import load
from yt.testing import fake_random_ds, requires_file
from yt.utilities.exceptions import YTBoundsDefinitionError


def test_cic_deposit():
    ds = fake_random_ds(64, nprocs=8, particles=64**3)
    my_reg = ds.arbitrary_grid(
        ds.domain_left_edge, ds.domain_right_edge, dims=[1, 800, 800]
    )
    f = ("deposit", "all_cic")
    assert_raises(YTBoundsDefinitionError, my_reg.__getitem__, f)


RAMSES = "output_00080/info_00080.txt"
RAMSES_small = "ramses_new_format/output_00002/info_00002.txt"
ISOGAL = "IsolatedGalaxy/galaxy0030/galaxy0030"


@requires_file(RAMSES)
def test_one_zone_octree_deposit():
    ds = load(RAMSES)

    # Get a sphere centred on the main halo
    hpos = ds.arr(
        [0.5215110772898429, 0.5215110772898429, 0.5215110772898429], "code_length"
    )
    hrvir = ds.quan(0.042307235300540924, "Mpc")

    sp = ds.sphere(hpos, hrvir * 10)
    assert sp["deposit", "io_cic"].shape == (1,)


@requires_file(RAMSES)
@requires_file(ISOGAL)
def test_mesh_sampling():
    for fn in (RAMSES, ISOGAL):
        ds = yt.load(fn)
        ds.add_mesh_sampling_particle_field(("index", "x"), ptype="all")
        ds.add_mesh_sampling_particle_field(("index", "dx"), ptype="all")

        dx = ds.r["all", "cell_index_dx"]
        xc = ds.r["all", "cell_index_x"]
        xp = ds.r["all", "particle_position_x"]

        dist = xp - xc

        assert_array_less(dist, dx)
        assert_array_less(-dist, dx)


@requires_file(RAMSES)
@requires_file(ISOGAL)
def test_mesh_sampling_for_filtered_particles():
    for fn in (RAMSES, ISOGAL):
        ds = yt.load(fn)

        @yt.particle_filter(requires=["particle_position_x"], filtered_type="io")
        def left(pfilter, data):
            return (
                data[(pfilter.filtered_type, "particle_position_x")].to("code_length")
                < 0.5
            )

        ds.add_particle_filter("left")

        for f in (("index", "x"), ("index", "dx"), ("gas", "density")):
            ds.add_mesh_sampling_particle_field(f, ptype="io")
            ds.add_mesh_sampling_particle_field(f, ptype="left")

        data_sources = (ds.all_data(), ds.box([0] * 3, [0.1] * 3))

        def test_source(ptype, src):
            # Test accessing
            src[ptype, "cell_index_x"]
            src[ptype, "cell_index_dx"]
            src[ptype, "cell_gas_density"]

        for ptype in ("io", "left"):
            for src in data_sources:
                test_source(ptype, src)


@requires_file(RAMSES)
def test_mesh_sampling_with_indexing():
    # Access with index caching
    ds = yt.load(RAMSES)
    ds.add_mesh_sampling_particle_field(("gas", "density"), ptype="all")

    ad = ds.all_data()
    ad["all", "cell_index"]
    v1 = ad["all", "cell_gas_density"]

    # Access with no index caching
    ds = yt.load(RAMSES)
    ds.add_mesh_sampling_particle_field(("gas", "density"), ptype="all")

    ad = ds.all_data()
    v2 = ad["all", "cell_gas_density"]

    # Check same answer is returned
    assert_allclose(v1, v2)


@requires_file(RAMSES_small)
def test_mesh_sampling_vs_field_value_at_point():
    all_ds = (fake_random_ds(ndims=3, particles=500), yt.load(RAMSES_small))

    for ds in all_ds:
        ds.add_mesh_sampling_particle_field(("gas", "density"), ptype="all")

        val = ds.r["all", "cell_gas_density"]
        ref = ds.find_field_values_at_points(
            ("gas", "density"), ds.r["all", "particle_position"]
        )

        assert_allclose(val, ref)
