import tempfile

from yt.testing import fake_particle_ds, requires_module


@requires_module("firefly_api")
def test_firefly_JSON_object():
    tmpdir = tempfile.mkdtemp()

    ds = fake_particle_ds()
    ad = ds.all_data()
    reader = ad.create_firefly_object(
        path_to_firefly=tmpdir,
        velocity_units="cm/s",
        coordinate_units="cm",
        dataset_name="test",
    )
    reader.dumpToJSON()
