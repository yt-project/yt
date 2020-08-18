from collections import OrderedDict

import pytest

from yt.testing import fake_amr_ds, fake_random_ds, \
    fake_hexahedral_ds, fake_tetrahedral_ds
from yt.visualization.geo_plot_utils import transform_list


def simple_contour(plot_field, plot):
    plot.annotate_contour(plot_field)

def simple_velocity(plot_field, plot):
    plot.annotate_velocity()

CALLBACK_TESTS = (
    ("simple_contour", (simple_contour,)),
    ("simple_velocity", (simple_velocity,)),
)


def get_attr_pairs(attr_args):
    pairs = []
    for k in attr_args:
        for v in attr_args[k]:
            pairs.append((k, v))
    return pairs


FPROPS = {
        'family': 'sans-serif',
        'style': 'italic',
        'weight': 'bold',
        'size': 24
}

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
             "toggle_right_handed": [((),{})]
}


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


PROFILE_ATTR_ARGS = {"annotate_text": [(((5e-29, 5e7), "Hello YT"), {}),
                               (((5e-29, 5e7), "Hello YT"), {'color': 'b'})],

             "set_title": [(('cell_mass', 'A phase plot.'), {})],
             "set_log": [(('cell_mass', False), {})],
             "set_unit": [(('cell_mass', 'Msun'), {})],
             "set_xlim": [((1e-27, 1e-24), {})],
             "set_ylim": [((1e2, 1e6), {})]
}


attr_pairs = get_attr_pairs(ATTR_ARGS)
part_attr_pairs = get_attr_pairs(PROJ_ATTR_ARGS)
phase_attr_pairs = get_attr_pairs(PHASE_ATTR_ARGS)
profile_pairs = get_attr_pairs(PROFILE_ATTR_ARGS)

attr_ids = [
    'pan',
    'pan_rel',
    'axis_unit_kpc',
    'axis_unit_mpc',
    'axis_unit_kpc_kpc',
    'axis_unit_kpc_mpc',
    'buff_size_1600',
    'buff_size_600_800',
    'center',
    'cmap_RdBu',
    'cmap_kamae',
    'font',
    'log',
    'window_size',
    'zlim_upper',
    'zlim_unbound', 
    'zoom',
    'toggle_rh'
]

phase_ids = [
    'annotate_no_color',
    'annotate_color',
    'title',
    'log',
    'unit',
    'xlim',
    'ylim'
]

profile_ids = [
    'annotate_Hello_YT',
    'annotate_Hello_YTb',
    'set_title',
    'set_log',
    'set_unit',
    'set_xlim',
    'set_ylim'
]

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

pr_fields = [
    ('density', 'temperature'),
    ('density', 'velocity_x'),
    ('temperature', 'cell_mass'),
    ('density', 'radius'),
    ('velocity_magnitude', 'cell_mass')
]

ph_fields = [
    ('density', 'temperature', 'cell_mass'),
    ('density', 'velocity_x', 'cell_mass'),
    ('radius', 'temperature', 'velocity_magnitude')
]


# Test parameters. Format:
# {test1: {param1 : [(val1, val2,...), (id1, id2,...)], param2 : ...}, test2: ...}
test_params = {
    'test_attributes' : {
        'axis' : [('x', 'y', 'z'), ('x', 'y', 'z')],
        'attr_name, attr_args' : [attr_pairs, attr_ids],
        'callback' : [CALLBACK_TESTS, callback_ids]
    },
    'test_attributes_wt' : {
        'attr_name, attr_args' : [attr_pairs, attr_ids],
        'callback' : [CALLBACK_TESTS, callback_ids]
    },
    'test_particle_projection_answers' : {
        'axis' : [('x', 'y', 'z'), ('x', 'y', 'z')],
        'attr_name, attr_args' : [part_attr_pairs, attr_ids]
    },
    'test_particle_projection_filter' : {
        'axis' : [('x', 'y', 'z'), ('x', 'y', 'z')],
        'attr_args' : [(v for v in PROJ_ATTR_ARGS['set_log']), ('set_log',)]
    },
    'test_particle_phase_answers' : {
        'attr_name, attr_args' : [phase_attr_pairs, phase_ids]
    },
    'test_raw_field_slices' : {
        'field' : [_raw_field_names, ('Bx', 'By', 'Bz', 'Ex', 'Ey', 'Ez', 'jx', 'jy', 'jz')]
    },
    'test_mesh_slices_amr' : {
        'field' : [fake_amr_ds().field_list, (f.__repr__() for f in fake_amr_ds().field_list)]
    },
    'test_mesh_slices_tetrahedral' : {
        'field' : [fake_tetrahedral_ds().field_list, (f.__repr__() for f in fake_amr_ds().field_list)],
        'idir' : [(0, 1, 2), ('0', '1', '2')]
    },
    'test_mesh_slices_hexahedral' : {
        'field' : [fake_hexahedral_ds().field_list, (f.__repr__() for f in fake_amr_ds().field_list)],
        'idir' : [(0, 1, 2,), ('0', '1', '2')]
    },
    'test_phase_plot_attributes' : {
        'attr_name, args' : [profile_pairs, profile_ids]
    },
}


def pytest_generate_tests(metafunc):
    if metafunc.function.__name__ == 'test_geo_slices_amr':
        ds = fake_amr_ds(geometry="geographic")
        metafunc.parametrize('transform',
            transform_list,
            ids=[t.__repr__() for t in transform_list]
        )
        metafunc.parametrize('field, ds',
            [(f, ds) for f in ds.field_list],
            ids=[f.__repr__() for f in ds.field_list]
        )
    elif metafunc.function.__name__ == 'test_profile_plot':
        metafunc.parametrize('region',
            ['region1', 'region_all_data'],
            indirect=True
        )
        metafunc.parametrize('x_field, y_field',
            [(f[0], f[1]) for f in pr_fields],
            ids=['dens-temp', 'dens-vx', 'temp-cell_mass', 'dens-rad', 'v-cell_mass']
        )
    elif metafunc.function.__name__ == 'test_phase_plot':
        metafunc.parametrize('region',
            ['region1', 'region_all_data'],
            indirect=True
        )
        metafunc.parametrize('x_field, y_field, z_field',
            [(f[0], f[1], f[2]) for f in ph_fields],
            ids=['dens-temp-cell_mass', 'dens-vx-cell_mass', 'rad-temp-v']
        )
    else:
        # Loop over each test in test_params
        for test_name, params in test_params.items():
            if metafunc.function.__name__ == test_name:
                # Parametrize
                for param_name, param_vals in params.items():
                    metafunc.parametrize(param_name, param_vals[0], ids=param_vals[1])

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
