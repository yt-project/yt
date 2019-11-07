_fields = (("halos", "particle_position_x"),
           ("halos", "particle_position_y"),
           ("halos", "particle_position_z"),
           ("halos", "particle_mass"))

methods = {"fof": 2, "hop": 2, "rockstar": 3}
           

def pytest_generate_tests(metafunc):
    if metafunc.function.__name__ == 'test_halo_finders':
        metafunc.parametrize('method, method_val', [(k, v) for k, v in methods.items()],
            ids=['fof', 'hop', 'rockstar'])
        metafunc.parametrize('field', _fields, ids=['x', 'y', 'z', 'mass'])
