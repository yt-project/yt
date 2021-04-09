import pytest

from yt.data_objects.static_output import Dataset
from yt.loaders import load, load_simulation
from yt.utilities.exceptions import (
    YTAmbiguousDataType,
    YTSimulationNotIdentified,
    YTUnidentifiedDataType,
)
from yt.utilities.object_registries import output_type_registry


def test_load_not_a_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        load(tmp_path / "not_a_file")


@pytest.mark.parametrize("simtype", ["Enzo", "unregistered_simulation_type"])
def test_load_simulation_not_a_file(simtype, tmp_path):
    # it is preferable to report the most important problem in an error message
    # (missing data is worse than a typo insimulation_type)
    # so we make sure the error raised is not YTSimulationNotIdentified,
    # even with an absurd simulation type
    with pytest.raises(FileNotFoundError):
        load_simulation(tmp_path / "not_a_file", simtype)


@pytest.fixture()
def tmp_path_with_empty_file(tmp_path):
    empty_file_path = tmp_path / "empty_file"
    empty_file_path.touch()
    return tmp_path, empty_file_path


def test_load_unidentified_data_dir(tmp_path_with_empty_file):
    tmp_path, empty_file_path = tmp_path_with_empty_file
    with pytest.raises(YTUnidentifiedDataType):
        load(tmp_path)


def test_load_unidentified_data_file(tmp_path_with_empty_file):
    tmp_path, empty_file_path = tmp_path_with_empty_file
    with pytest.raises(YTUnidentifiedDataType):
        load(empty_file_path)


def test_load_simulation_unidentified_data_dir(tmp_path_with_empty_file):
    tmp_path, empty_file_path = tmp_path_with_empty_file
    with pytest.raises(YTSimulationNotIdentified):
        load_simulation(tmp_path, "unregistered_simulation_type")


def test_load_simulation_unidentified_data_file(tmp_path_with_empty_file):
    tmp_path, empty_file_path = tmp_path_with_empty_file
    with pytest.raises(YTSimulationNotIdentified):
        load_simulation(
            empty_file_path,
            "unregistered_simulation_type",
        )


def test_load_ambiguous_data(tmp_path):
    # we deliberately setup a situation where two Dataset subclasses
    # that aren't parents are consisdered valid
    class FakeDataset(Dataset):
        @classmethod
        def _is_valid(cls, *args, **kwargs):
            return True

    class FakeDataset2(Dataset):
        @classmethod
        def _is_valid(cls, *args, **kwargs):
            return True

    try:
        with pytest.raises(YTAmbiguousDataType):
            load(tmp_path)
    finally:
        # tear down to avoid possible breakage in following tests
        output_type_registry.pop("FakeDataset")
        output_type_registry.pop("FakeDataset2")
