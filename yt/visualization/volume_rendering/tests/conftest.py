import pytest

from yt.testing import \
    fake_hexahedral_ds


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


@pytest.fixture(scope='class')
def ds_hex():
    ds = fake_hexahedral_ds()
    return ds
