import os
import shutil
import tempfile
from unittest import TestCase

import numpy as np

from yt.testing import assert_almost_equal, assert_equal, fake_random_ds


def setup():
    from yt.config import ytcfg

    ytcfg["yt", "internals", "within_testing"] = True


def test_flux_calculation():
    ds = fake_random_ds(64, nprocs=4)
    dd = ds.all_data()
    surf = ds.surface(dd, ("index", "x"), 0.51)
    assert_equal(surf[("index", "x")], 0.51)
    flux = surf.calculate_flux(
        ("index", "ones"), ("index", "zeros"), ("index", "zeros"), ("index", "ones")
    )
    assert_almost_equal(flux.value, 1.0, 12)
    assert_equal(str(flux.units), "cm**2")
    flux2 = surf.calculate_flux(
        ("index", "ones"), ("index", "zeros"), ("index", "zeros")
    )
    assert_almost_equal(flux2.value, 1.0, 12)
    assert_equal(str(flux2.units), "cm**2")


def test_sampling():
    ds = fake_random_ds(64, nprocs=4)
    dd = ds.all_data()
    for i, ax in enumerate([("index", "x"), ("index", "y"), ("index", "z")]):
        surf = ds.surface(dd, ax, 0.51)
        surf.get_data(ax, sample_type="vertex")
        assert_equal(surf.vertex_samples[ax], surf.vertices[i, :])
        assert_equal(str(surf.vertices.units), "code_length")
        dens = surf[("gas", "density")]
        vert_shape = surf.vertices.shape
        assert_equal(dens.shape[0], vert_shape[1] // vert_shape[0])
        assert_equal(str(dens.units), "g/cm**3")


class ExporterTests(TestCase):
    def setUp(self):
        self.curdir = os.getcwd()
        self.tmpdir = tempfile.mkdtemp()
        os.chdir(self.tmpdir)

    def tearDown(self):
        os.chdir(self.curdir)
        shutil.rmtree(self.tmpdir)

    def test_export_ply(self):
        ds = fake_random_ds(64, nprocs=4)
        dd = ds.all_data()
        surf = ds.surface(dd, ("index", "x"), 0.51)
        surf.export_ply("my_ply.ply", bounds=[(0, 1), (0, 1), (0, 1)])
        assert os.path.exists("my_ply.ply")
        surf.export_ply(
            "my_ply2.ply",
            bounds=[(0, 1), (0, 1), (0, 1)],
            sample_type="vertex",
            color_field=("gas", "density"),
        )
        assert os.path.exists("my_ply2.ply")

    def test_export_obj(self):
        ds = fake_random_ds(
            16,
            nprocs=4,
            particles=16**3,
            fields=("density", "temperature"),
            units=("g/cm**3", "K"),
        )
        sp = ds.sphere("max", (1.0, "cm"))
        surf = ds.surface(sp, ("gas", "density"), 0.5)
        surf.export_obj("my_galaxy", transparency=1.0, dist_fac=1.0)
        assert os.path.exists("my_galaxy.obj")
        assert os.path.exists("my_galaxy.mtl")

        mi, ma = sp.quantities.extrema(("gas", "temperature"))
        rhos = [0.5, 0.25]
        trans = [0.5, 1.0]
        for i, r in enumerate(rhos):
            basename = "my_galaxy_color"
            surf = ds.surface(sp, ("gas", "density"), r)
            surf.export_obj(
                basename,
                transparency=trans[i],
                color_field=("gas", "temperature"),
                dist_fac=1.0,
                plot_index=i,
                color_field_max=ma,
                color_field_min=mi,
            )

            assert os.path.exists(f"{basename}.obj")
            assert os.path.exists(f"{basename}.mtl")

        def _Emissivity(field, data):
            return (
                data["gas", "density"]
                * data["gas", "density"]
                * np.sqrt(data["gas", "temperature"])
            )

        ds.add_field(
            ("gas", "emissivity"),
            sampling_type="cell",
            function=_Emissivity,
            units=r"g**2*sqrt(K)/cm**6",
        )
        for i, r in enumerate(rhos):
            basename = "my_galaxy_emis"
            surf = ds.surface(sp, ("gas", "density"), r)
            surf.export_obj(
                basename,
                transparency=trans[i],
                color_field=("gas", "temperature"),
                emit_field=("gas", "emissivity"),
                dist_fac=1.0,
                plot_index=i,
            )

            basename = "my_galaxy_emis"
            assert os.path.exists(f"{basename}.obj")
            assert os.path.exists(f"{basename}.mtl")


def test_correct_output_unit_fake_ds():
    # see issue #1368
    ds = fake_random_ds(64, nprocs=4, particles=16**3)
    x = y = z = 0.5
    sp1 = ds.sphere((x, y, z), (300, "kpc"))
    Nmax = sp1.max(("gas", "density"))
    sur = ds.surface(sp1, ("gas", "density"), 0.5 * Nmax)
    sur[("index", "x")][0]


def test_radius_surface():
    # see #1407
    ds = fake_random_ds(64, nprocs=4, particles=16**3, length_unit=10.0)
    reg = ds.all_data()
    sp = ds.sphere(ds.domain_center, (0.5, "code_length"))
    for obj in [reg, sp]:
        for rad in [0.05, 0.1, 0.4]:
            surface = ds.surface(obj, ("index", "radius"), (rad, "code_length"))
            assert_almost_equal(surface.surface_area.v, 4 * np.pi * rad**2, decimal=2)
            verts = surface.vertices
            for i in range(3):
                assert_almost_equal(verts[i, :].min().v, 0.5 - rad, decimal=2)
                assert_almost_equal(verts[i, :].max().v, 0.5 + rad, decimal=2)
