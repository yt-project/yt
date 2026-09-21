import warnings

import numpy as np
import pytest
from numpy.testing import assert_equal

from yt.testing import fake_random_ds
from yt.visualization.api import SlicePlot
from yt.visualization.image_writer import apply_colormap


@pytest.mark.parametrize(
    "image, color_bounds",
    [
        # explicit degenerate bounds, as built by the grid annotation callbacks
        # ([0, max_level]) for any dataset whose only grid level is 0
        (np.zeros(4), [0, 0]),
        (np.array([0.0, 1.0, 2.0]), [0, 0]),
        (np.array([2.0, 2.0]), [2, 2]),
        # inferred bounds collapse the same way for a constant image
        (np.full((8, 8), 3.0), None),
        (np.zeros((8, 8)), None),
    ],
)
def test_apply_colormap_degenerate_bounds(image, color_bounds):
    # a zero-width normalization range must not be divided by
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        to_plot = apply_colormap(
            image, color_bounds=color_bounds, cmap_name="B-W LINEAR_r"
        )

    assert not np.isnan(to_plot).any()
    # everything lands at the bottom of the colormap, which is what
    # matplotlib.colors.Normalize does when vmin == vmax
    expected = apply_colormap(
        np.zeros(np.shape(image), dtype="float64"),
        color_bounds=[0, 1],
        cmap_name="B-W LINEAR_r",
    )
    assert_equal(to_plot, expected)


def test_annotate_grids_single_level_dataset():
    # a dataset with a single grid level makes the callback's color bounds
    # degenerate ([0, 0]), which used to divide by zero and render the grid
    # edges as transparent black
    ds = fake_random_ds(16)
    assert ds.index.max_level == 0

    p = SlicePlot(ds, "z", ("gas", "density"))
    p.annotate_grids()
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        p.render()

    collections = [
        c for c in p.plots["gas", "density"].axes.collections if len(c.get_paths())
    ]
    assert collections
    edgecolors = np.concatenate([c.get_edgecolors() for c in collections])
    assert len(edgecolors)
    assert not np.isnan(edgecolors).any()
    # the grids have to actually be visible
    assert not np.allclose(edgecolors[:, :3], 0.0)
