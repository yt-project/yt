import pytest

from yt.sample_data.api import get_data_registry_table


@pytest.mark.parametrize("column", ["load_name", "hash", "load_kwargs"])
def test_registy_integrity(column):
    reg = get_data_registry_table()
    missing_values = reg[reg[column].isna()]
    assert missing_values.empty
