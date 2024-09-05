import importlib.metadata
import sys

import pytest

import yt
from yt.data_objects.static_output import Dataset
from yt.geometry.grid_geometry_handler import GridIndex
from yt.utilities.object_registries import output_type_registry


class MockEntryPoint:
    @classmethod
    def load(cls):
        class MockHierarchy(GridIndex):
            grid = None

        class ExtDataset(Dataset):
            _index_class = MockHierarchy

            def _parse_parameter_file(self):
                self.current_time = 1.0
                self.cosmological_simulation = 0

            def _set_code_unit_attributes(self):
                self.length_unit = self.quan(1.0, "code_length")
                self.mass_unit = self.quan(1.0, "code_mass")
                self.time_unit = self.quan(1.0, "code_time")

            @classmethod
            def _is_valid(cls, filename, *args, **kwargs):
                return filename.endswith("mock")


@pytest.fixture()
def mock_external_frontend(monkeypatch):
    def mock_entry_points(group=None):
        if sys.version_info >= (3, 10):
            return [MockEntryPoint]
        else:
            return {"yt.frontends": [MockEntryPoint]}

    monkeypatch.setattr(importlib.metadata, "entry_points", mock_entry_points)
    assert "ExtDataset" not in output_type_registry

    yield

    assert "ExtDataset" in output_type_registry
    # teardown to avoid test pollution
    output_type_registry.pop("ExtDataset")


@pytest.mark.usefixtures("mock_external_frontend")
def test_external_frontend(tmp_path):
    test_file = tmp_path / "tmp.mock"
    test_file.write_text("")  # create the file
    assert test_file.is_file()

    ds = yt.load(test_file)
    assert "ExtDataset" in ds.__class__.__name__
