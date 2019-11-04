import pytest


_fields = (("halos", "particle_position_x"),
           ("halos", "particle_position_y"),
           ("halos", "particle_position_z"),
           ("halos", "particle_mass"))

rockstar_files = ["rockstar_halos/halos_0.0.bin", "rockstar_halos/halos_1.0.bin"]


def pytest_generate_tests(metafunc):
    if metafunc.function.__name__ == 'test_rockstar':
        metafunc.parametrize('h', rockstar_files, ids=['halos_0.0.bin', 'halos_1.0.bin'])
        metafunc.parametrize('field', _fields, ids=['x', 'y', 'z', 'mass'])
