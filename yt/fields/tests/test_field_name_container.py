from yt import load
from yt.testing import fake_amr_ds, fake_hexahedral_ds, requires_file


def do_field_type(ft):
    assert dir(ft) == sorted(dir(ft))
    assert sorted(dir(ft)) == sorted(f.name[1] for f in ft)
    for field_name in dir(ft):
        f = getattr(ft, field_name)
        assert (ft.field_type, field_name) == f.name
    for field in ft:
        f = getattr(ft, field.name[1])
        assert f == field
        assert f in ft
        assert f.name in ft
        assert f.name[1] in ft


enzotiny = "enzo_tiny_cosmology/DD0046/DD0046"


@requires_file(enzotiny)
def test_field_name_container():
    ds = load(enzotiny)
    assert dir(ds.fields) == sorted(dir(ds.fields))
    assert sorted(ft.field_type for ft in ds.fields) == sorted(dir(ds.fields))
    for field_type in dir(ds.fields):
        assert field_type in ds.fields
        ft = getattr(ds.fields, field_type)
        do_field_type(ft)
    for field_type in ds.fields:
        assert field_type in ds.fields
        do_field_type(field_type)


def test_vertex_fields_only_in_unstructured_ds():
    def get_vertex_fields(ds):
        return [(ft, fn) for ft, fn in ds.derived_field_list if "vertex" in fn]

    ds = fake_amr_ds()
    vertex_fields = get_vertex_fields(ds)
    assert not vertex_fields

    ds = fake_hexahedral_ds()
    actual = get_vertex_fields(ds)
    expected = [
        ("all", "vertex_x"),
        ("all", "vertex_y"),
        ("all", "vertex_z"),
        ("connect1", "vertex_x"),
        ("connect1", "vertex_y"),
        ("connect1", "vertex_z"),
        ("index", "vertex_x"),
        ("index", "vertex_y"),
        ("index", "vertex_z"),
    ]
    assert actual == expected
