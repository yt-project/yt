from yt import \
    load
from yt.testing import \
    requires_file

def do_field_type(ft):
    for field_name in dir(ft):
        f = getattr(ft, field_name)
        assert ((ft.field_type, field_name) == f.name)
    for field in ft:
        f = getattr(ft, field.name[1])
        assert (f == field)
        assert (f in ft)
        assert (f.name in ft)
        assert (f.name[1] in ft)


enzotiny = "enzo_tiny_cosmology/DD0046/DD0046"
@requires_file(enzotiny)
def test_field_name_container():
    ds = load(enzotiny)
    for field_type in dir(ds.fields):
        assert (field_type in ds.fields)
        ft = getattr(ds.fields, field_type)
        do_field_type(ft)
    for field_type in ds.fields:
        assert (field_type in ds.fields)
        do_field_type(field_type)
