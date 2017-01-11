import os
import shutil
import tempfile

import numpy as np

from unittest import TestCase

from yt.testing import \
    fake_random_ds, \
    assert_almost_equal, \
    assert_equal, \
    requires_file

from yt.convenience import \
    load

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

def test_flux_calculation():
    ds = fake_random_ds(64, nprocs = 4)
    dd = ds.all_data()
    surf = ds.surface(dd, "x", 0.51)
    assert_equal(surf["x"], 0.51)
    flux = surf.calculate_flux("ones", "zeros", "zeros", "ones")
    assert_almost_equal(flux.value, 1.0, 12)
    assert_equal(str(flux.units), 'cm**2')

def test_sampling():
    ds = fake_random_ds(64, nprocs = 4)
    dd = ds.all_data()
    for i, ax in enumerate('xyz'):
        surf = ds.surface(dd, ax, 0.51)
        surf.get_data(ax, "vertex")
        assert_equal(surf.vertex_samples[ax], surf.vertices[i,:])
        assert_equal(str(surf.vertices.units), 'code_length')
        dens = surf['density']
        vert_shape = surf.vertices.shape
        assert_equal(dens.shape[0], vert_shape[1]//vert_shape[0])
        assert_equal(str(dens.units), 'g/cm**3')

ISOGAL = 'IsolatedGalaxy/galaxy0030/galaxy0030'

class ExporterTests(TestCase):

    def setUp(self):
        self.curdir = os.getcwd()
        self.tmpdir = tempfile.mkdtemp()
        os.chdir(self.tmpdir)

    def tearDown(self):
        os.chdir(self.curdir)
        shutil.rmtree(self.tmpdir)

    def test_export_ply(self):
        ds = fake_random_ds(64, nprocs = 4)
        dd = ds.all_data()
        surf = ds.surface(dd, 'x', 0.51)
        surf.export_ply('my_ply.ply', bounds=[(0, 1), (0, 1), (0, 1)])
        assert os.path.exists('my_ply.ply')

    @requires_file(ISOGAL)
    def test_export_obj(self):
        ds = load(ISOGAL)
        sp = ds.sphere("max", (10, "kpc"))
        trans = 1.0
        distf = 3.1e18*1e3 # distances into kpc
        surf = ds.surface(sp, "density", 5e-27)
        surf.export_obj("my_galaxy", transparency=trans, dist_fac = distf)
        assert os.path.exists('my_galaxy.obj')
        assert os.path.exists('my_galaxy.mtl')

        mi, ma = sp.quantities.extrema('temperature')
        rhos = [1e-24, 1e-25]
        trans = [0.5, 1.0]
        for i, r in enumerate(rhos):
            surf = ds.surface(sp,'density',r)
            surf.export_obj("my_galaxy_color".format(i),
                            transparency=trans[i],
                            color_field='temperature', dist_fac = distf,
                            plot_index = i, color_field_max = ma,
                            color_field_min = mi)

        assert os.path.exists('my_galaxy_color.obj')
        assert os.path.exists('my_galaxy_color.mtl')

        def _Emissivity(field, data):
            return (data['density']*data['density'] *
                    np.sqrt(data['temperature']))
        ds.add_field("emissivity", sampling_type='cell', function=_Emissivity,
                     units=r"g**2*sqrt(K)/cm**6")
        for i, r in enumerate(rhos):
            surf = ds.surface(sp,'density',r)
            surf.export_obj("my_galaxy_emis".format(i),
                            transparency=trans[i],
                            color_field='temperature',
                            emit_field='emissivity',
                            dist_fac = distf, plot_index = i)

        assert os.path.exists('my_galaxy_emis.obj')
        assert os.path.exists('my_galaxy_emis.mtl')
