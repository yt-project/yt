"""
Tests for Dyablo frontend.
"""

import os

import h5py
import numpy as np
import numpy.typing as npt
import pytest

import yt
from yt.utilities.lib.geometry_utils import compute_morton


def create_test_file(
    tmpdir,
    n_blocks: npt.NDArray[np.int_],
    block_size: npt.NDArray[np.int_],
    left_edge: npt.NDArray[np.float64],
    right_edge: npt.NDArray[np.float64],
    seed: int = 42,
):
    """Create a small test HDF5 file for Dyablo."""
    n_blocks = np.asarray(n_blocks)
    block_size = np.asarray(block_size)
    left_edge = np.asarray(left_edge)
    right_edge = np.asarray(right_edge)

    iteration = 1234567

    fluid_fname = os.path.join(tmpdir, f"test_dyablo_iter{iteration}.h5")

    n_cells = np.prod(n_blocks) * np.prod(block_size)

    generator = np.random.default_rng(seed=seed)

    # Compute vertices of the blocks
    vertex_coordinates = np.mgrid[
        left_edge[0] : right_edge[0] : (block_size[0] * n_blocks[0] + 1) * 1j,
        left_edge[1] : right_edge[1] : (block_size[1] * n_blocks[1] + 1) * 1j,
        left_edge[2] : right_edge[2] : (block_size[2] * n_blocks[2] + 1) * 1j,
    ]

    # Block positions
    dx = (right_edge - left_edge) / n_blocks
    block_inds = np.mgrid[0 : n_blocks[0], 0 : n_blocks[1], 0 : n_blocks[2]]
    block_positions = (
        left_edge[:, None, None, None] + (block_inds + 0.5) * dx[:, None, None, None]
    )

    # Compute block order **along Morton curve**
    # Dyablo's Morton order uses x as the least-significant bit; yt's
    # compute_morton uses z as LSB, so we swap x and z.
    morton_index = compute_morton(
        block_positions[2].flatten(),
        block_positions[1].flatten(),
        block_positions[0].flatten(),
        left_edge[::-1],
        right_edge[::-1],
    )
    block_order = np.argsort(morton_index)

    # Compute cell ordering within each block **in column-major (Fortran) order**.
    # For each F-order output position, find the corresponding C-order column index
    # in cell_inds.reshape(3,-1). This is the F→C permutation (not C→F).
    ix_f, iy_f, iz_f = np.unravel_index(
        np.arange(np.prod(block_size)), block_size, order="F"
    )

    # Offset w.r.t. the bottom-left-front corner of each cell
    vertex_offset = np.asarray(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ]
    ).T

    vertex_inds = (
        block_inds.reshape(3, -1)[:, block_order][:, :, None, None]
        * block_size[:, None, None, None]
        + np.array([ix_f, iy_f, iz_f])[:, None, :, None]
        + vertex_offset[:, None, None, :]
    )
    vertex_order = np.ravel_multi_index(
        vertex_inds, vertex_coordinates.shape[1:], order="C"
    )

    with h5py.File(fluid_fname, "w") as f:
        # Create fluid fields
        f["rho"] = generator.random(n_cells)
        f["rho_vx"] = generator.random(n_cells)
        f["rho_vy"] = generator.random(n_cells)
        f["rho_vz"] = generator.random(n_cells)
        f["e_tot"] = generator.random(n_cells)
        f["test_field"] = generator.random(n_cells)

        # Create vertices of the blocks
        f["coordinates"] = vertex_coordinates.reshape(3, -1).T
        f["connectivity"] = vertex_order.reshape(-1, 8)

        # Create scalar_data
        f.create_group("scalar_data")
        f["scalar_data"].attrs.update({"aexp": 1.0, "iter": iteration, "time": 123.0})


# Expected fluid fields written by create_test_file
_FLUID_FIELDS = ["rho", "rho_vx", "rho_vy", "rho_vz", "e_tot", "test_field"]


@pytest.mark.parametrize(
    "n_blocks, block_size",
    [
        # cubic, minimal 2×2×2 (block_size=1 is ambiguous to inference, start from 2)
        ([2, 2, 2], [2, 2, 2]),
        ([2, 2, 2], [4, 4, 4]),
        # non-cubic (N≠M≠L), uniform n_blocks
        ([4, 4, 4], [3, 4, 5]),
        # asymmetric n_blocks (all dims ≥ 2; avoid z-neighbor-as-block-1 orderings)
        ([3, 4, 6], [4, 4, 4]),
        ([2, 4, 6], [3, 5, 7]),
    ],
)
def test_dyablo_mock_dataset(tmp_path, n_blocks, block_size):
    n_blocks = np.asarray(n_blocks)
    block_size = np.asarray(block_size)

    create_test_file(
        str(tmp_path),
        n_blocks=n_blocks,
        block_size=block_size,
        left_edge=np.zeros(3),
        right_edge=np.ones(3),
    )

    fpath = next(tmp_path.glob("*.h5"))
    yt.mylog.setLevel(50)
    ds = yt.load(str(fpath))

    # Check inferred block_size and n_blocks
    np.testing.assert_array_equal(
        ds.block_size,
        block_size,
        err_msg=f"block_size mismatch: got {ds.block_size}, expected {block_size}",
    )
    np.testing.assert_array_equal(
        ds.domain_dimensions,
        n_blocks * block_size,
        err_msg=f"domain_dimensions mismatch: got {ds.domain_dimensions}, expected {n_blocks * block_size}",
    )

    ad = ds.all_data()
    expected_n_cells = int(np.prod(n_blocks) * np.prod(block_size))

    # Check all fluid fields are readable and have correct size
    for fname in _FLUID_FIELDS:
        data = ad["dyablo", fname]
        assert data.shape == (expected_n_cells,), (
            f"Field {fname}: expected {expected_n_cells} cells, got {data.shape}"
        )
        assert not np.any(np.isnan(data)), f"Field {fname} contains NaN values"

    # Check derived density field
    rho = ad["gas", "density"]
    assert rho.shape == (expected_n_cells,)
    assert not np.any(np.isnan(rho))


def test_dyablo_mock_dataset_edges(tmp_path):
    """Test that arbitrary left_edge and right_edge are preserved correctly."""
    n_blocks = np.asarray([2, 2, 2])
    block_size = np.asarray([4, 4, 4])
    left_edge = np.array([-3.0, 0.5, 1.0])
    right_edge = np.array([5.0, 4.5, 7.0])

    create_test_file(
        str(tmp_path),
        n_blocks=n_blocks,
        block_size=block_size,
        left_edge=left_edge,
        right_edge=right_edge,
    )

    fpath = next(tmp_path.glob("*.h5"))
    yt.mylog.setLevel(50)
    ds = yt.load(str(fpath))

    np.testing.assert_allclose(ds.domain_left_edge.d, left_edge)
    np.testing.assert_allclose(ds.domain_right_edge.d, right_edge)

    ad = ds.all_data()
    rho = ad["gas", "density"]
    assert rho.shape == (int(np.prod(n_blocks) * np.prod(block_size)),)
    assert not np.any(np.isnan(rho))
