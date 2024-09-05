from copy import deepcopy

import pytest

import yt
from yt.config import ytcfg
from yt.testing import requires_file

output_00080 = "output_00080/info_00080.txt"
ramses_new_format = "ramses_new_format/output_00002/info_00002.txt"

custom_hydro = [
    "my-Density",
    "my-x-velocity",
    "my-y-velocity",
    "my-z-velocity",
    "my-Pressure",
    "my-Metallicity",
]
custom_grav = [
    "my-x-acceleration",
    "my-y-acceleration",
    "my-z-acceleration",
]


@pytest.fixture()
def custom_ramses_fields_conf():
    old_config = deepcopy(ytcfg.config_root.as_dict())
    ytcfg.add_section("ramses-hydro")
    ytcfg.add_section("ramses-grav")
    ytcfg.set("ramses-hydro", "fields", custom_hydro)
    ytcfg.set("ramses-grav", "fields", custom_grav)
    yield
    ytcfg.remove_section("ramses-hydro")
    ytcfg.remove_section("ramses-grav")
    ytcfg.update(old_config)


@requires_file(output_00080)
def test_field_config_1(custom_ramses_fields_conf):
    ds = yt.load(output_00080)

    for f in custom_hydro:
        assert ("ramses", f) in ds.field_list
    for f in custom_grav:
        assert ("gravity", f) in ds.field_list


@requires_file(ramses_new_format)
def test_field_config_2(custom_ramses_fields_conf):
    ds = yt.load(ramses_new_format)

    for f in custom_hydro:
        assert ("ramses", f) in ds.field_list
    for f in custom_grav:
        assert ("gravity", f) in ds.field_list
