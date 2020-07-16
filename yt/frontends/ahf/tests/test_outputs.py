import os.path

from yt.frontends.ahf.api import AHFHalosDataset
from yt.testing import ParticleSelectionComparison, assert_equal, requires_file
from yt.utilities.answer_testing.framework import (
    FieldValuesTest,
    data_dir_load,
    requires_ds,
)

_fields = (
    "particle_position_x",
    "particle_position_y",
    "particle_position_z",
    "particle_mass",
)

ahf_halos = "ahf_halos/snap_N64L16_135.parameter"


def load(filename):
    return data_dir_load(filename, kwargs={"hubble_constant": 0.7})


@requires_ds(ahf_halos)
def test_fields_ahf_halos():
    ds = load(ahf_halos)
    assert_equal(str(ds), os.path.basename(ahf_halos))
    for field in _fields:
        yield FieldValuesTest(ds, field, particle_type=True)


@requires_file(ahf_halos)
def test_AHFHalosDataset():
    ds = load(ahf_halos)
    assert isinstance(ds, AHFHalosDataset)
    ad = ds.all_data()
    ad["particle_mass"]
    psc = ParticleSelectionComparison(ds)
    psc.run_defaults()
