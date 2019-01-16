import sys
if sys.version.startswith("2"):
    import pathlib2 as pathlib
else:
    import pathlib
import yt

from yt.utilities.parameter_file_storage import output_type_registry

test_dir = pathlib.Path(__file__).parent.resolve()
example_data_2D = str(test_dir / "data/blast_wave_Cartesian_2D/bw_2d0004.dat")

def test_frontends_contents():
    """Check that our frontend is correctly being registered."""
    assert "amrvac" in dir(yt.frontends)

def test_output_type_registry():
    """Check that the Dataset child class is part of classes tried for validation."""
    keys = list(output_type_registry.keys())
    assert "AMRVACDataset" in keys

def test_is_amrvac_valid():
    """Check that our format validation fonction works at least on example file"""
    assert pathlib.Path(example_data_2D).is_file()
    assert output_type_registry["AMRVACDataset"]._is_valid(example_data_2D)

def test_load_amrvac():
    """Check that yt.load() doesn't crash AND correctly guesses that data is amrvac-formated"""
    ds = yt.load(example_data_2D)
    assert ds.__class__ == yt.frontends.amrvac.AMRVACDataset
