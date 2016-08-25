from yt.testing import fake_random_ds, assert_equal
from yt.data_objects.profiles import create_profile
from yt.visualization.plot_window import \
    SlicePlot, \
    ProjectionPlot, \
    OffAxisProjectionPlot
from yt.visualization.profile_plotter import \
    ProfilePlot, \
    PhasePlot

def test_field_access():
    ds = fake_random_ds(16)

    ad = ds.all_data()
    sp = ds.sphere(ds.domain_center, 0.25)
    cg = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    scg = ds.smoothed_covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    sl = ds.slice(0, ds.domain_center[0])
    proj = ds.proj('density', 0)
    prof = create_profile(ad, 'radius', 'density')

    for data_object in [ad, sp, cg, scg, sl, proj, prof]:
        assert_equal(
            data_object['gas', 'density'],
            data_object[ds.fields.gas.density]
        )

    for field in [('gas', 'density'), ds.fields.gas.density]:
        ad = ds.all_data()
        prof = ProfilePlot(ad, 'radius', field)
        phase = PhasePlot(ad, 'radius', field, 'cell_mass')
        s = SlicePlot(ds, 2, field)
        oas = SlicePlot(ds, [1, 1, 1], field)
        p = ProjectionPlot(ds, 2, field)
        oap = OffAxisProjectionPlot(ds, [1, 1, 1], field)

        for plot_object in [s, oas, p, oap, prof, phase]:
            plot_object._setup_plots()
            if hasattr(plot_object, '_frb'):
                plot_object._frb[field]
