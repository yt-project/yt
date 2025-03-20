import matplotlib
import pytest

from yt.visualization.base_plot_types import BACKEND_SPECS


@pytest.mark.skipif(
    matplotlib.__version_info__ < (3, 9),
    reason="This test requires matplotlib 3.9 or higher.",
)
@pytest.mark.parametrize("backend_key", BACKEND_SPECS.keys())
def test_backend_specs(backend_key):
    # note: while this test itself requires matplotlib 3.9, it is
    # testing functionality in yt that can be removed once yt has
    # a minimal matplotlib version of 3.9,
    # see https://github.com/yt-project/yt/issues/5138
    from matplotlib.backends import backend_registry

    assert backend_registry.is_valid_backend(backend_key)
