from collections import defaultdict
from itertools import product

from yt.data_objects.profiles import create_profile
from yt.testing import assert_equal, assert_raises, fake_random_ds
from yt.utilities.exceptions import YTFieldNotFound
from yt.visualization.plot_window import (
    OffAxisProjectionPlot,
    ProjectionPlot,
    SlicePlot,
)
from yt.visualization.profile_plotter import PhasePlot, ProfilePlot


def test_field_access():
    ds = fake_random_ds(16)

    ad = ds.all_data()
    sp = ds.sphere(ds.domain_center, 0.25)
    cg = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    scg = ds.smoothed_covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    sl = ds.slice(0, ds.domain_center[0])
    proj = ds.proj("density", 0)
    prof = create_profile(ad, "radius", "density")

    for data_object in [ad, sp, cg, scg, sl, proj, prof]:
        assert_equal(data_object["gas", "density"], data_object[ds.fields.gas.density])

    for field in [("gas", "density"), ds.fields.gas.density]:
        ad = ds.all_data()
        prof = ProfilePlot(ad, "radius", field)
        phase = PhasePlot(ad, "radius", field, "cell_mass")
        s = SlicePlot(ds, 2, field)
        oas = SlicePlot(ds, [1, 1, 1], field)
        p = ProjectionPlot(ds, 2, field)
        oap = OffAxisProjectionPlot(ds, [1, 1, 1], field)

        for plot_object in [s, oas, p, oap, prof, phase]:
            plot_object._setup_plots()
            if hasattr(plot_object, "_frb"):
                plot_object._frb[field]


def test_unexisting_field_access():
    ds = fake_random_ds(16, particles=10)

    fname2ftype = defaultdict(list)

    for ft, fn in ds.derived_field_list:
        fname2ftype[fn].append(ft)

    ad = ds.all_data()

    ftypes = ("gas", "io")
    fnames = (
        "density",
        "particle_position_x",
        "particle_position_y",
        "particle_position_z",
    )

    # Try invalid ftypes, fnames combinations
    for ft, fn in product(ftypes, fnames):
        if (ft, fn) in ds.derived_field_list:
            continue
        with assert_raises(YTFieldNotFound) as ex:
            ad[(ft, fn)]

        # Make sure the existing field has been suggested
        for possible_ft in fname2ftype[fn]:
            assert (possible_ft, fn) in ex.exception.suggestions

    # Try typos
    for bad_field, good_field in (
        (("gas", "densi_y"), ("gas", "density")),
        (("oi", "particle_mass"), ("io", "particle_mass")),
        (("gas", "DENSITY"), ("gas", "density")),
    ):
        with assert_raises(YTFieldNotFound) as ex:
            ad[bad_field]

        assert good_field in ex.exception.suggestions
