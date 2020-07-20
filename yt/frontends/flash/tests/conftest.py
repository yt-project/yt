"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest


# Test parameters. Format:
# {test1: {param1 : [(val1, val2,...), (id1, id2,...)], param2 : ...}, test2: ...}
test_params = {
    'test_sloshing' : {
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('max', (0.1, 'unitary')))), ('None' , 'sphere')],
        'w' : [(None, 'density'), ('None' , 'density')],
        'f' : [('temperature', 'density', 'velocity_magnitude'),
                ('temperature', 'density', 'velocity_magnitude')]
    },
    'test_wind_tunnel' : {
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('max', (0.1, 'unitary')))), ('None' , 'sphere')],
        'w' : [(None, 'density'), ('None' , 'density')],
        'f' : [('temperature', 'density'), ('temperature', 'density')]
    },
    'test_fid_1to3_b1' : {
        'f, w' : [((("all", "particle_mass") , None),
                    (("all", "particle_ones") , None),
                    (("all", "particle_velocity_x") , ("all", "particle_mass")),
                    (("all", "particle_velocity_y") , ("all", "particle_mass")),
                    (("all", "particle_velocity_z") , ("all", "particle_mass"))),
                    ('all_mass', 'all_ones', 'all_vx', 'all_vy', 'all_vz')],
        'd' : [(None, ('sphere', ('c', (0.1, 'unitary')))), ('None' ,'sphere')],
        'a' : [(0, 1, 2), ('0', '1', '2')]
    }
}


def pytest_generate_tests(metafunc):
    for test_name, params in test_params.items():
        if metafunc.function.__name__ == test_name:
            for param_name, param_vals in params.items():
                metafunc.parametrize(param_name, param_vals[0], ids=param_vals[1])
