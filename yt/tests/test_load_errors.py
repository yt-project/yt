import os
import tempfile
from pathlib import Path

from yt.convenience import load, simulation
from yt.data_objects.static_output import Dataset
from yt.testing import assert_raises
from yt.utilities.exceptions import (
    YTAmbiguousDataType,
    YTOutputNotIdentified,
    YTSimulationNotIdentified,
)
from yt.utilities.object_registries import output_type_registry


def test_load_nonexistent_data():
    with tempfile.TemporaryDirectory() as tmpdir:
        assert_raises(FileNotFoundError, load, os.path.join(tmpdir, "not_a_file"))
        assert_raises(
            FileNotFoundError, simulation, os.path.join(tmpdir, "not_a_file"), "Enzo"
        )

        # this one is a design choice:
        # it is preferable to report the most important problem in an error message
        # (missing data is worse than a typo insimulation_type)
        # so we make sure the error raised is not YTSimulationNotIdentified
        assert_raises(
            FileNotFoundError,
            simulation,
            os.path.join(tmpdir, "not_a_file"),
            "unregistered_simulation_type",
        )


def test_load_unidentified_data():
    with tempfile.TemporaryDirectory() as tmpdir:
        empty_file_path = Path(tmpdir) / "empty_file"
        empty_file_path.touch()
        assert_raises(YTOutputNotIdentified, load, tmpdir)
        assert_raises(YTOutputNotIdentified, load, empty_file_path)
        assert_raises(
            YTSimulationNotIdentified,
            simulation,
            tmpdir,
            "unregistered_simulation_type",
        )
        assert_raises(
            YTSimulationNotIdentified,
            simulation,
            empty_file_path,
            "unregistered_simulation_type",
        )


def test_load_ambiguous_data():
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
        with tempfile.TemporaryDirectory() as tmpdir:
            assert_raises(YTAmbiguousDataType, load, tmpdir)
    except Exception:
        raise
    finally:
        # tear down to avoid possible breakage in following tests
        output_type_registry.pop("FakeDataset")
        output_type_registry.pop("FakeDataset2")
