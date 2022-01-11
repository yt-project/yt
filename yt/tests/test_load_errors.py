import pytest

from yt.data_objects.static_output import Dataset
from yt.geometry.grid_geometry_handler import GridIndex
from yt.loaders import load, load_simulation
from yt.utilities.exceptions import (
    YTAmbiguousDataType,
    YTSimulationNotIdentified,
    YTUnidentifiedDataType,
)
from yt.utilities.object_registries import output_type_registry


@pytest.fixture
def tmp_data_dir(tmp_path):
    from yt.config import ytcfg

    pre_test_data_dir = ytcfg["yt", "test_data_dir"]
    ytcfg.set("yt", "test_data_dir", str(tmp_path))

    yield tmp_path

    ytcfg.set("yt", "test_data_dir", pre_test_data_dir)


@pytest.mark.usefixtures("tmp_data_dir")
def test_load_not_a_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        load(tmp_path / "not_a_file")


@pytest.mark.parametrize("simtype", ["Enzo", "unregistered_simulation_type"])
@pytest.mark.usefixtures("tmp_data_dir")
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


@pytest.fixture()
def ambiguous_dataset_classes():
    # We deliberately setup a situation where two Dataset subclasses
    # that aren't parents are consisdered valid.
    # We implement the bare minimum for these classes to be actually
    # loadable in order to test hints.
    class MockHierarchy(GridIndex):
        pass

    class MockDataset(Dataset):
        _index_class = MockHierarchy

        def _parse_parameter_file(self, *args, **kwargs):
            self.current_time = -1.0
            self.cosmological_simulation = 0

        def _set_code_unit_attributes(self, *args, **kwargs):
            self.length_unit = self.quan(1, "m")
            self.mass_unit = self.quan(1, "kg")
            self.time_unit = self.quan(1, "s")

    class AlphaDataset(MockDataset):
        @classmethod
        def _is_valid(cls, *args, **kwargs):
            return True

    class BetaDataset(MockDataset):
        @classmethod
        def _is_valid(cls, *args, **kwargs):
            return True

    yield

    # teardown to avoid possible breakage in following tests
    output_type_registry.pop("MockDataset")
    output_type_registry.pop("AlphaDataset")
    output_type_registry.pop("BetaDataset")


@pytest.mark.usefixtures("ambiguous_dataset_classes")
def test_load_ambiguous_data(tmp_path):
    with pytest.raises(YTAmbiguousDataType):
        load(tmp_path)

    file = tmp_path / "fake_datafile0011.dump"
    file.touch()

    pattern = str(tmp_path / "fake_datafile00??.dump")

    # loading a DatasetSeries should not crash until an item is retrieved
    ts = load(pattern)
    with pytest.raises(YTAmbiguousDataType):
        ts[0]


@pytest.mark.parametrize(
    "hint, expected_type",
    [
        ("alpha", "AlphaDataset"),
        ("al", "AlphaDataset"),
        ("ph", "AlphaDataset"),
        ("beta", "BetaDataset"),
        ("BeTA", "BetaDataset"),
        ("b", "BetaDataset"),
    ],
)
@pytest.mark.usefixtures("ambiguous_dataset_classes")
def test_load_ambiguous_data_with_hint(hint, expected_type, tmp_path):
    ds = load(tmp_path, hint=hint)
    assert type(ds).__name__ == expected_type

    file1 = tmp_path / "fake_datafile0011.dump"
    file2 = tmp_path / "fake_datafile0022.dump"
    file1.touch()
    file2.touch()

    pattern = str(tmp_path / "fake_datafile00??.dump")

    ts = load(pattern, hint=hint)
    ds = ts[0]
    assert type(ds).__name__ == expected_type

    ds = ts[1]
    assert type(ds).__name__ == expected_type
