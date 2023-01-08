import sys
import tempfile
from unittest import mock

from yt.data_objects.static_output import Dataset
from yt.geometry.grid_geometry_handler import GridIndex
from yt.loaders import load
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

            def _set_code_unit_attributes(self):
                self.length_unit = self.quan(1.0, "code_length")
                self.mass_unit = self.quan(1.0, "code_mass")
                self.time_unit = self.quan(1.0, "code_time")

            @classmethod
            def _is_valid(cls, filename, *args, **kwargs):
                prefix = tempfile.gettempdir()
                return filename.startswith(prefix) and filename.endswith("mock")


def test_external_frontend():
    assert "ExtDataset" not in output_type_registry

    with tempfile.NamedTemporaryFile(suffix="mock") as test_file:
        with mock.patch("yt.loaders.entry_points") as ep:
            if sys.version_info >= (3, 10):
                ep.return_value = [MockEntryPoint]
            else:
                ep.return_value = {"yt.frontends": [MockEntryPoint]}
            ds = load(test_file.name)
            assert "ExtDataset" in str(ds.__class__)

    assert "ExtDataset" in output_type_registry
    output_type_registry.pop("ExtDataset")
