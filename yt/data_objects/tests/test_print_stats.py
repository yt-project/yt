import pytest

from yt.testing import (
    fake_amr_ds,
    fake_particle_ds,
)


def test_no_print_stats():
    ds = fake_particle_ds()
    with pytest.raises(NotImplementedError):
        ds.print_stats()


def test_print_stats():
    ds = fake_amr_ds()
    ds.print_stats()
