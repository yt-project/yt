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
        'f, w' : [((("deposit", "all_density") , None),
                    (("deposit", "all_count") , None),
                    (("deposit", "all_cic") , None),
                    (("deposit", "all_cic_velocity_x") , ("deposit", "all_cic")),
                    (("deposit", "all_cic_velocity_y") , ("deposit", "all_cic")),
                    (("deposit", "all_cic_velocity_z") , ("deposit", "all_cic"))),
                    ('all_dens', 'all_count', 'all_cic', 'all_cic_vx', 'all_cic_vy', 'all_cic_vz')],
        'd' : [(None, ('sphere', ('c', (0.1, 'unitary')))), ('None' ,'sphere')],
        'a' : [(0, 1, 2), ('0', '1', '2')]
    }
}


def pytest_generate_tests(metafunc):
    for test_name, params in test_params.items():
        if metafunc.function.__name__ == test_name:
            for param_name, param_vals in params.items():
                metafunc.parametrize(param_name, param_vals[0], ids=param_vals[1])
