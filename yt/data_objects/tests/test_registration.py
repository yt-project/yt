import pytest

from yt.data_objects.static_output import Dataset
from yt.utilities.object_registries import output_type_registry


def test_reregistration_warning():
    true_EnzoDataset = output_type_registry["EnzoDataset"]
    try:
        with pytest.warns(
            UserWarning,
            match=(
                "Overwritting EnzoDataset, which was previously registered. "
                "This is expected if you're importing a yt extension with a "
                "frontend that was already migrated to the main code base."
            ),
        ):

            class EnzoDataset(Dataset):
                pass

    finally:
        output_type_registry["EnzoDataset"] = true_EnzoDataset
