import os.path

from yt.testing import assert_equal
from yt.utilities.answer_testing.framework import (
    FieldValuesTest,
    data_dir_load,
    requires_ds,
)

# from yt.frontends.owls_subfind.api import OWLSSubfindDataset

_fields = (
    ("all", "particle_position_x"),
    ("all", "particle_position_y"),
    ("all", "particle_position_z"),
    ("all", "particle_mass"),
)

# a dataset with empty files
g1 = "owls_fof_halos/groups_001/group_001.0.hdf5"
g8 = "owls_fof_halos/groups_008/group_008.0.hdf5"


@requires_ds(g8)
def test_fields_g8():
    ds = data_dir_load(g8)
    assert_equal(str(ds), os.path.basename(g8))
    for field in _fields:
        yield FieldValuesTest(g8, field, particle_type=True)


@requires_ds(g1)
def test_fields_g1():
    ds = data_dir_load(g1)
    assert_equal(str(ds), os.path.basename(g1))
    for field in _fields:
        yield FieldValuesTest(g1, field, particle_type=True)
