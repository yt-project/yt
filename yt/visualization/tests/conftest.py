import pytest


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


def get_attr_pairs(attr_args):
    pairs = []
    for k in attr_args:
        for v in args[k]:
            pairs.append((k, v))
    return pairs

attr_pairs = get_attr_pairs(ATTR_ARGS)
part_attr_pairs = get_attr_pairs(PROJ_ATTR_ARGS)
phase_attr_pairs = get_attr_pairs(PHASE_ATTR_ARGS)

attr_ids = ['pan', 'pan_rel', 'axis_unit_kpc', 'axis_unit_mpc', 'axis_unit_kpc_kpc',
    'axis_unit_kpc_mpc', 'buff_size_1600', 'buff_size_600_800', 'center', 'cmap_RdBu',
    'cmap_kamae', 'font', 'log', 'window_size', 'zlim_upper', 'zlim_unbound', 'zlim_dynamic',
    'zoom', 'toggle_rh']
phase_ids = ['annotate_no_color', 'annotate_color', 'title', 'log', 'unit', 'xlim', 'ylim']
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
        metafunc.parametrize('attr_name, attr_args', part_attr_pairs, ids=part_attr_ids)
    if metafunc.function.__name__ == "test_particle_projection_filter":
        metafunc.parametrize('axis', ['x', 'y', 'z'], ids=['x', 'y', 'z'])
        attr_name = "set_log"
        metafunc.parametrize('attr_args', [v for v in PROJ_ATTR_ARGS[attr_name]], ids=part_attr_ids)
    if metafunc.function.__name__ == "test_particle_phase_answers":
        metafunc.parametrize('attr_name, attr_args', phase_attr_pairs, ids=phase_ids)
    if metafunc.function.__name__ == "test_raw_field_slices":
        metafunc.parametrize('field', _raw_field_names, ids=['Bx', 'By', 'Bz',
            'Ex', 'Ey', 'Ez', 'jx', 'jy', 'jz'])
