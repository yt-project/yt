import os
import tempfile
from pathlib import Path

from yt.convenience import load, simulation
from yt.testing import assert_raises
from yt.utilities.exceptions import YTOutputNotIdentified, YTSimulationNotIdentified


def test_load_nonexistent_data():
    with tempfile.TemporaryDirectory() as tmpdir:
        assert_raises(FileNotFoundError, load, os.path.join(tmpdir, "not_a_file"))
        assert_raises(
            FileNotFoundError, simulation, os.path.join(tmpdir, "not_a_file"), "Enzo"
        )

        # this one is a design choice: it is preferable to report the most important
        # problem in an error message (missing data is worse than a typo in
        # simulation_type), so we make sure the error raised is not YTSimulationNotIdentified
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
