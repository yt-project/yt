from yt.utilities.answer_testing.framework import \
    requires_ds, \
    data_dir_load

enzotiny = "enzo_tiny_cosmology/DD0046/DD0046"
@requires_ds(enzotiny)
def test_simulated_halo_mass_function():
    ds = data_dir_load(enzotiny)
    for field_type in dir(ds.fields):
        ft = getattr(ds.fields, field_type)
        for field_name in dir(ft):
            f = getattr(ft, field_name)
            assert ((field_type, field_name) == f.name)
