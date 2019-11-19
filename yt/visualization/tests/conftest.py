from collections import OrderedDict

import pytest

from yt.testing import fake_amr_ds, fake_random_ds
from yt.visualization.geo_plot_utils import transform_list


def simple_contour(plot_field, plot):
    plot.annotate_contour(plot_field)

def simple_velocity(plot_field, plot):
    plot.annotate_velocity()


FPROPS = {'family': 'sans-serif', 'style': 'italic',
          'weight': 'bold', 'size': 24}

ATTR_ARGS = {"pan": [(((0.1, 0.1), ), {})],
             "pan_rel": [(((0.1, 0.1), ), {})],
             "set_axes_unit": [(("kpc", ), {}),
                               (("Mpc", ), {}),
                               ((("kpc", "kpc"),), {}),
                               ((("kpc", "Mpc"),), {})],
             "set_buff_size": [((1600, ), {}),
                               (((600, 800), ), {})],
             "set_center": [(((0.4, 0.3), ), {})],
             "set_cmap": [(('density', 'RdBu'), {}),
                          (('density', 'kamae'), {})],
             "set_font": [((OrderedDict(sorted(FPROPS.items(), key=lambda t: t[0])), ),
                           {})],
             "set_log": [(('density', False), {})],
             "set_window_size": [((7.0, ), {})],
             "set_zlim": [(('density', 1e-25, 1e-23), {}),
                          (('density', 1e-25, None), {'dynamic_range': 4})],
             "zoom": [((10, ), {})],
             "toggle_right_handed": [((),{})]}

# override some of the plotwindow ATTR_ARGS
PROJ_ATTR_ARGS = ATTR_ARGS.copy() 
PROJ_ATTR_ARGS["set_cmap"] = [(('particle_mass', 'RdBu'), {}), 
                                  (('particle_mass', 'kamae'), {})]
PROJ_ATTR_ARGS["set_log"] = [(('particle_mass', False), {})]
PROJ_ATTR_ARGS["set_zlim"] = [(('particle_mass', 1e-25, 1e-23), {}),
                                  (('particle_mass', 1e-25, None), 
                                   {'dynamic_range': 4})]

PHASE_ATTR_ARGS = {"annotate_text": [(((5e-29, 5e7), "Hello YT"), {}), 
                               (((5e-29, 5e7), "Hello YT"), {'color':'b'})],
                   "set_title": [(('particle_mass', 'A phase plot.'), {})],
                   "set_log": [(('particle_mass', False), {})],
                   "set_unit": [(('particle_mass', 'Msun'), {})],
                   "set_xlim": [((-4e7, 4e7), {})],
                   "set_ylim": [((-4e7, 4e7), {})]}

CALLBACK_TESTS = (
    ("simple_contour", (simple_contour,)),
    ("simple_velocity", (simple_velocity,)),
)

PROFILE_ATTR_ARGS = {"annotate_text": [(((5e-29, 5e7), "Hello YT"), {}),
                               (((5e-29, 5e7), "Hello YT"), {'color': 'b'})],

             "set_title": [(('cell_mass', 'A phase plot.'), {})],
             "set_log": [(('cell_mass', False), {})],
             "set_unit": [(('cell_mass', 'Msun'), {})],
             "set_xlim": [((1e-27, 1e-24), {})],
             "set_ylim": [((1e2, 1e6), {})]}


def get_attr_pairs(attr_args):
    pairs = []
    for k in attr_args:
        for v in attr_args[k]:
            pairs.append((k, v))
    return pairs

attr_pairs = get_attr_pairs(ATTR_ARGS)
part_attr_pairs = get_attr_pairs(PROJ_ATTR_ARGS)
phase_attr_pairs = get_attr_pairs(PHASE_ATTR_ARGS)
profile_pairs = get_attr_pairs(PROFILE_ATTR_ARGS)

attr_ids = ['pan', 'pan_rel', 'axis_unit_kpc', 'axis_unit_mpc', 'axis_unit_kpc_kpc',
    'axis_unit_kpc_mpc', 'buff_size_1600', 'buff_size_600_800', 'center', 'cmap_RdBu',
    'cmap_kamae', 'font', 'log', 'window_size', 'zlim_upper', 'zlim_unbound', 
    'zoom', 'toggle_rh']
phase_ids = ['annotate_no_color', 'annotate_color', 'title', 'log', 'unit', 'xlim', 'ylim']
profile_ids = ['annotate_Hello_YT', 'annotate_Hello_YTb', 'set_title', 'set_log', 'set_unit',
    'set_xlim', 'set_ylim']
callback_ids = ['contour', 'velocity']

_raw_field_names =  [('raw', 'Bx'),
                     ('raw', 'By'),
                     ('raw', 'Bz'),
                     ('raw', 'Ex'),
                     ('raw', 'Ey'),
                     ('raw', 'Ez'),
                     ('raw', 'jx'),
                     ('raw', 'jy'),
                     ('raw', 'jz')]


def pytest_generate_tests(metafunc):
    if metafunc.function.__name__ == "test_attributes":
        metafunc.parametrize('axis', ['x', 'y', 'z'], ids=['x', 'y', 'z'])
        metafunc.parametrize('attr_name, attr_args', attr_pairs, ids=attr_ids)
        metafunc.parametrize('callback', CALLBACK_TESTS, ids=callback_ids)
    if metafunc.function.__name__ == "test_attributes_wt":
        metafunc.parametrize('attr_name, attr_args', attr_pairs, ids=attr_ids)
        metafunc.parametrize('callback', CALLBACK_TESTS, ids=callback_ids)
    if metafunc.function.__name__ == "test_particle_projection_answers":
        metafunc.parametrize('axis', ['x', 'y', 'z'], ids=['x', 'y', 'z'])
        metafunc.parametrize('attr_name, attr_args', part_attr_pairs, ids=attr_ids)
    if metafunc.function.__name__ == "test_particle_projection_filter":
        metafunc.parametrize('axis', ['x', 'y', 'z'], ids=['x', 'y', 'z'])
        attr_name = "set_log"
        metafunc.parametrize('attr_args', [v for v in PROJ_ATTR_ARGS[attr_name]], ids=["set_log"])
    if metafunc.function.__name__ == "test_particle_phase_answers":
        metafunc.parametrize('attr_name, attr_args', phase_attr_pairs, ids=phase_ids)
    if metafunc.function.__name__ == "test_raw_field_slices":
        metafunc.parametrize('field', _raw_field_names, ids=['Bx', 'By', 'Bz',
            'Ex', 'Ey', 'Ez', 'jx', 'jy', 'jz'])
    if metafunc.function.__name__ == 'test_geo_slices_amr':
        ds = fake_amr_ds(geometry="geographic")
        metafunc.parametrize('transform', transform_list, ids=[t.__repr__() for t in transform_list])
        metafunc.parametrize('field, ds', [(f, ds) for f in ds.field_list], ids=[f.__repr__() for f in ds.field_list])
    # The parameterizations for test_mesh_slices.py's functions doesn't work b/c
    # the ds_amr fixture cannot be accesed from here. Creating the ds object here
    # gives issues with weak references. I'm not sure how to implement the indirect
    # method used in test_profile_plot.py's parameterizations below since the
    # fields are unknown. These tests are skipped, however, so this doesn't
    # currently break anything
    if metafunc.function.__name__ == 'test_mesh_slices_amr':
        metafunc.parametrize('field', ds_amr.field_list, ids=[f.__repr__() for f in fields])
    if metafunc.function.__name__ == 'test_mesh_slices_tetrahedral':
        metafunc.parametrize('field', ds_tetra.field_list, ids=[f.__repr__() for f in fields])
        metafunc.parametrize('idir', [0, 1, 2], ids=['0', '1', '2'])
    if metafunc.function.__name__ == 'test_mesh_slices_hexahedral':
        metafunc.parametrize('field', ds_hex.field_list, ids=[f.__repr__() for f in fields])
        metafunc.parametrize('idir', [0, 1, 2], ids=['0', '1', '2'])
    if metafunc.function.__name__ == 'test_phase_plot_attributes':
        metafunc.parametrize('attr_name, args', profile_pairs, ids=profile_ids)
    if metafunc.function.__name__ == 'test_profile_plot':
        pr_fields = [('density', 'temperature'), ('density', 'velocity_x'),
                     ('temperature', 'cell_mass'), ('density', 'radius'),
                     ('velocity_magnitude', 'cell_mass')]
        metafunc.parametrize('region', ['region1', 'region_all_data'], indirect=True)
        metafunc.parametrize('x_field, y_field', [(f[0], f[1]) for f in pr_fields],
            ids=['dens-temp', 'dens-vx', 'temp-cell_mass', 'dens-rad', 'v-cell_mass'])
    if metafunc.function.__name__ == 'test_phase_plot':
        ph_fields = [('density', 'temperature', 'cell_mass'),
                     ('density', 'velocity_x', 'cell_mass'),
                     ('radius', 'temperature', 'velocity_magnitude')]
        metafunc.parametrize('region', ['region1', 'region_all_data'], indirect=True)
        metafunc.parametrize('x_field, y_field, z_field', [(f[0], f[1], f[2]) for f in ph_fields],
            ids=['dens-temp-cell_mass', 'dens-vx-cell_mass', 'rad-temp-v'])


@pytest.fixture(scope='class')
def ds_amr():
    ds = fake_amr_ds()
    return ds

@pytest.fixture(scope='class')
def ds_tetra():
    ds = fake_tetrahedral_ds()
    return ds

@pytest.fixture(scope='class')
def ds_hex():
    ds = fake_hexahedral_ds()
    return ds

@pytest.fixture(scope='class')
def ds_random():
    ds = fake_random_ds(16, fields=('density', 'temperature'))
    return ds

@pytest.fixture(scope='class')
def ds_test():
    fields = ('density', 'temperature', 'velocity_x', 'velocity_y', 'velocity_z')
    units = ('g/cm**3', 'K', 'cm/s', 'cm/s', 'cm/s')
    ds = fake_random_ds(16, fields=fields, units=units)
    return ds

@pytest.fixture(scope='class')
def ds_mult():
    fields = ('density', 'temperature', 'dark_matter_density')
    ds = fake_random_ds(16, fields=fields)
    return ds

@pytest.fixture(scope='class')
def region(request, ds_test):
    if request.param == 'region1':
        return ds_test.region([0.5] * 3, [0.4] * 3, [0.6] * 3)
    if request.param == 'region_all_data':
        return ds_test.all_data()
