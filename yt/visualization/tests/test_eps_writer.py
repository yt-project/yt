import yt
from yt.testing import fake_amr_ds, requires_external_executable, requires_module


@requires_external_executable("tex")
@requires_module("pyx")
def test_eps_writer(tmp_path):
    import yt.visualization.eps_writer as eps

    fields = [
        ("gas", "density"),
        ("gas", "temperature"),
    ]
    units = [
        "g/cm**3",
        "K",
    ]
    ds = fake_amr_ds(fields=fields, units=units)
    slc = yt.SlicePlot(
        ds,
        "z",
        fields=fields,
    )
    eps_fig = eps.multiplot_yt(2, 1, slc, bare_axes=True)
    eps_fig.scale_line(0.2, "5 cm")
    savefile = tmp_path / "multi"
    eps_fig.save_fig(savefile, format="eps")
    assert savefile.with_suffix(".eps").exists()
