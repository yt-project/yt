import os.path

from numpy.testing import assert_equal

from yt.frontends.rockstar.api import RockstarDataset
from yt.testing import ParticleSelectionComparison, requires_file
from yt.utilities.answer_testing.framework import (
    FieldValuesTest,
    data_dir_load,
    requires_ds,
)

_fields = (
    ("all", "particle_position_x"),
    ("all", "particle_position_y"),
    ("all", "particle_position_z"),
    ("all", "particle_mass"),
)

r1 = "rockstar_halos/halos_0.0.bin"


@requires_ds(r1)
def test_fields_r1():
    ds = data_dir_load(r1)
    assert_equal(str(ds), os.path.basename(r1))
    for field in _fields:
        yield FieldValuesTest(r1, field, particle_type=True)


@requires_file(r1)
def test_RockstarDataset():
    assert isinstance(data_dir_load(r1), RockstarDataset)


@requires_file(r1)
def test_particle_selection():
    ds = data_dir_load(r1)
    psc = ParticleSelectionComparison(ds)
    psc.run_defaults()


@requires_file(r1)
def test_halo_loading():
    ds = data_dir_load(r1)

    for halo_id, Npart in zip(
        ds.r["halos", "particle_identifier"],
        ds.r["halos", "num_p"],
        strict=False,
    ):
        halo = ds.halo("halos", halo_id)
        assert halo is not None

        # Try accessing properties
        halo.position
        halo.velocity
        halo.mass

        # Make sure we can access the member particles
        assert_equal(len(halo.member_ids), Npart)
