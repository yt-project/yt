"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest

from yt.testing import assert_equal
from yt.utilities.answer_testing import utils


_orion_fields = ("temperature", "density", "velocity_magnitude")
_nyx_fields = ("Ne", "Temp", "particle_mass_density")
_castro_fields = ("Temp", "density", "particle_count")
_warpx_fields = ("Ex", "By", "jz")
_raw_fields = [('raw', 'Bx'), ('raw', 'Ey'), ('raw', 'jz')]


# Test parameters. Format:
# {test1: {param1 : [(val1, val2,...), (id1, id2,...)], param2 : ...}, test2: ...}
test_params = {
    'test_radadvect' : {
        'f' : [_orion_fields, ('temperature', 'density', 'velocity_magnitude')],
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('max', (0.1, 'unitary')))), ('None', 'sphere')],
        'w' : [(None, 'density'), ('None', 'density')]
    },
    'test_radtube' : {
        'f' : [_orion_fields, ('temperature', 'density', 'velocity_magnitude')],
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('max', (0.1, 'unitary')))), ('None', 'sphere')],
        'w' : [(None, 'density'), ('None', 'density')]
    },
    'test_star' : {
        'f' : [_orion_fields, ('temperature', 'density', 'velocity_magnitude')],
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('max', (0.1, 'unitary')))), ('None', 'sphere')],
        'w' : [(None, 'density'), ('None', 'density')]
    },
    'test_LyA' : {
        'f' : [_nyx_fields, ('Ne', 'Temp', 'particle_mass_density')],
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('c', (0.1, 'unitary')))), ('None', 'sphere')],
        'w' : [(None, 'Ne'), ('None', 'Ne')]
    },
    'test_RT_particles' : {
        'f' : [_castro_fields, ('Temp', 'density', 'count')],
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('max', (0.1, 'unitary')))), ('None', 'sphere')],
        'w' : [(None, 'density'), ('None', 'density')]
    },
    'test_langmuir' : {
        'f' : [_warpx_fields, ('Ex', 'By', 'jz')],
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('c', (0.1, 'unitary')))), ('None', 'sphere')],
        'w' : [(None, 'Ex'), ('None', 'Ex')]
    },
    'test_plasma' : {
        'f' : [_warpx_fields, ('Ex', 'By', 'jz')],
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('c', (0.1, 'unitary')))), ('None', 'sphere')],
        'w' : [(None, 'Ex'), ('None', 'Ex')]
    },
    'test_beam' : {
        'f' : [_warpx_fields, ('Ex', 'By', 'jz')],
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('c', (0.1, 'unitary')))), ('None', 'sphere')],
        'w' : [(None, 'Ex'), ('None', 'Ex')]
    },
    'test_raw_fields' : {
        'f' : [_raw_fields, ('Bx', 'Ey', 'jz')]
    }
}


def pytest_generate_tests(metafunc):
    # Loop over each test in test_params
    for test_name, params in test_params.items():
        if metafunc.function.__name__ == test_name:
            # Parametrize
            for param_name, param_vals in params.items():
                metafunc.parametrize(param_name, param_vals[0], ids=param_vals[1])
