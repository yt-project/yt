import numpy as np

import yt.utilities.metadata as metadata_mod


class _DummyDS:
    dimensionality = 3
    domain_dimensions = np.array([16, 8, 4])
    dataset_type = "dummy"
    omega_lambda = None

    def __init__(self):
        self.parameters = {
            "grid_shape": np.array([16, 8, 4]),
            "refine_by": 2,
        }
        self.closed = False

    def close(self):
        self.closed = True


def test_get_metadata_handles_lists_missing_values_and_full_output():
    dummy_ds = _DummyDS()

    original_load = metadata_mod.load
    try:
        metadata_mod.load = lambda path: dummy_ds

        metadata = metadata_mod.get_metadata(
            "dummy/output",
            full_output=True,
            attrs=("dimensionality", "domain_dimensions", "omega_lambda", "dataset_type"),
        )

        assert metadata == {
            "filename": "dummy/output",
            "dimensionality": 3,
            "domain_dimensions": [16, 8, 4],
            "dataset_type": "dummy",
            "params": {"grid_shape": [16, 8, 4], "refine_by": 2},
        }
        assert dummy_ds.closed is True
    finally:
        metadata_mod.load = original_load


def test_get_metadata_omits_params_when_full_output_is_false():
    dummy_ds = _DummyDS()

    original_load = metadata_mod.load
    try:
        metadata_mod.load = lambda path: dummy_ds

        metadata = metadata_mod.get_metadata("dummy/output", full_output=False)

        assert metadata["filename"] == "dummy/output"
        assert "params" not in metadata
        assert dummy_ds.closed is True
    finally:
        metadata_mod.load = original_load
