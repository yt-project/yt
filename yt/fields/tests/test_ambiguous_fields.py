import pytest

from yt.testing import fake_amr_ds


def test_ambiguous_fails():
    ds = fake_amr_ds(particles=10)
    msg = "The requested field name '{}' is ambiguous"

    def _mock_field(field, data):
        return data["ones"]

    # create a pair of ambiguous field
    ds.add_field(("io", "mock_field"), _mock_field, sampling_type="particle")
    ds.add_field(("gas", "mock_field"), _mock_field, sampling_type="cell")

    # Test errors are raised for ambiguous fields
    with pytest.raises(ValueError, match=msg.format("mock_field")):
        ds.r["mock_field"]

    # check that explicit name tuples don't raise a warning
    ds.r["io", "mock_field"]
    ds.r["gas", "mock_field"]


def test_nameonly_field_with_all_aliases_candidates():
    # see https://github.com/yt-project/yt/issues/3839
    ds = fake_amr_ds(fields=["density"], units=["g/cm**3"])

    # here we rely on implementations details from fake_amr_ds,
    # so we verify that it provides the appropriate conditions
    # for the actual test.
    candidates = [f for f in ds.derived_field_list if f[1] == "density"]
    assert len(candidates) == 2
    fi = ds.field_info
    assert fi[candidates[0]].is_alias_to(fi[candidates[1]])

    # this is the actual test (check that no error or warning is raised)
    ds.all_data()["density"]
