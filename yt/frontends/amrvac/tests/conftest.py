"""
Title: conftest.py
Purpose: Generates parameters and loads data for tests.
"""
import pytest

from yt.utilities.answer_testing.utils import get_parameterization


bw_polar_2d = get_parameterization("amrvac/bw_polar_2D0000.dat")
bw_cart_3d = get_parameterization("amrvac/bw_3d0000.dat")
bw_sph_2d = get_parameterization("amrvac/bw_2d0000.dat")
bw_cyl_3d = get_parameterization("amrvac/bw_cylindrical_3D0000.dat")
khi_cart_2d = get_parameterization("amrvac/kh_2d0000.dat")
khi_cart_3d = get_parameterization("amrvac/kh_3D0000.dat")
jet_cyl_25d = get_parameterization("amrvac/Jet0003.dat")
rie_cart_175d = get_parameterization("amrvac/R_1d0005.dat")

test_params = {
    'test_bw_polar_2d' : {
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('max', (0.1, 'unitary')))), ('None', 'sphere')],
        'w' : [(None, 'density'), ('None', 'density')],
        'f' : [bw_polar_2d[0], bw_polar_2d[1]]
    },
    'test_blastwave_cartesian_3D' : {
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('max', (0.1, 'unitary')))), ('None', 'sphere')],
        'w' : [(None, 'density'), ('None', 'density')],
        'f' : [bw_cart_3d[0], bw_cart_3d[1]]
    },
    'test_blastwave_spherical_2D' : {
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('max', (0.1, 'unitary')))), ('None', 'sphere')],
        'w' : [(None, 'density'), ('None', 'density')],
        'f' : [bw_sph_2d[0], bw_sph_2d[1]]
    },
    'test_blastwave_cylindrical_3D' : {
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('max', (0.1, 'unitary')))), ('None', 'sphere')],
        'w' : [(None, 'density'), ('None', 'density')],
        'f' : [bw_cyl_3d[0], bw_cyl_3d[1]]
    },
    'test_khi_cartesian_2D' : {
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('max', (0.1, 'unitary')))), ('None', 'sphere')],
        'w' : [(None, 'density'), ('None', 'density')],
        'f' : [khi_cart_2d[0], khi_cart_2d[1]]
    },
    'test_khi_cartesian_3D' : {
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('max', (0.1, 'unitary')))), ('None', 'sphere')],
        'w' : [(None, 'density'), ('None', 'density')],
        'f' : [khi_cart_3d[0], khi_cart_3d[1]]
    },
    'test_jet_cylindrical_25D' : {
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('max', (0.1, 'unitary')))), ('None', 'sphere')],
        'w' : [(None, 'density'), ('None', 'density')],
        'f' : [jet_cyl_25d[0], jet_cyl_25d[1]]
    },
    'test_riemann_cartesian_175D' : {
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('max', (0.1, 'unitary')))), ('None', 'sphere')],
        'w' : [(None, 'density'), ('None', 'density')],
        'f' : [rie_cart_175d[0], rie_cart_175d[1]]
    },
}


def pytest_generate_tests(metafunc):
    for test_name, params in test_params.items():
        if metafunc.function.__name__ == test_name:
            for param_name, param_vals in params.items():
                metafunc.parametrize(param_name, param_vals[0], ids=param_vals[1])
