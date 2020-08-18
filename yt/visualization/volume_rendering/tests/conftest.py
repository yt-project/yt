import pytest

from yt.testing import \
    fake_hexahedral_ds, \
    fake_vr_orientation_test_ds
from yt.visualization.volume_rendering.api import \
    Scene, \
    VolumeSource, \
    ColorTransferFunction


hex8_fields = [('connect1', 'diffused'), ('connect2', 'convected')]
tet4_fields = [("connect1", "u")]
hex20_fields = [('connect2', 'temp')]
wedge6_fields = [('connect1', 'diffused')]
tet10_fields = [('connect1', 'uz')]


def pytest_generate_tests(metafunc):
    if metafunc.function.__name__ == 'test_fake_hexahedral_ds_render':
        field_list = [('connect1', 'elem'), ('connect1', 'test')]
        metafunc.parametrize('field', field_list, ids=['elem', 'test'])
    if metafunc.function.__name__ in ['test_hex20_render_pyembree', 'test_hex20_render']:
        metafunc.parametrize('field', hex20_fields, ids=['temp'])
    if metafunc.function.__name__ in ['test_hex8_render_pyembree', 'test_hex8_render']:
        metafunc.parametrize('field', hex8_fields, ids=['diffused', 'convected'])
    if metafunc.function.__name__ in ['test_tet10_render_pyembree', 'test_tet10_render']:
        metafunc.parametrize('field', tet10_fields, ids=['uz'])
    if metafunc.function.__name__ in ['test_tet4_render_pyembree', 'test_tet4_render']:
        metafunc.parametrize('field', tet4_fields, ids=['u'])
    if metafunc.function.__name__ in ['test_wedge6_render_pyembree', 'test_wedge6_render']:
        metafunc.parametrize('field', wedge6_fields, ids=['diffused'])
    if metafunc.function.__name__ == 'test_vr_images':
        lens_types = ['plane-parallel', 'perspective']
        metafunc.parametrize('lens_type', lens_types, ids=['plane-parallel', 'perspective'])
    if metafunc.function.__name__ == 'test_orientations':
        orientations = [[-0.3, -0.1, 0.8]]
        metafunc.parametrize('orientation', orientations, ids=['(-0.3, -0.1, 0.8)'])


@pytest.fixture(scope='class')
def ds_hex():
    ds = fake_hexahedral_ds()
    return ds


@pytest.fixture(scope='class')
def ds_vr():
    ds = fake_vr_orientation_test_ds()
    return ds


@pytest.fixture(scope='class')
def sc(ds_vr):
    sc = Scene()
    vol = VolumeSource(ds_vr, field=('gas', 'density'))
    sc.add_source(vol)
    tf = vol.transfer_function
    tf = ColorTransferFunction((0.1, 1.0))
    tf.sample_colormap(1.0, 0.01, colormap="coolwarm")
    tf.sample_colormap(0.8, 0.01, colormap="coolwarm")
    tf.sample_colormap(0.6, 0.01, colormap="coolwarm")
    tf.sample_colormap(0.3, 0.01, colormap="coolwarm")
    return sc
