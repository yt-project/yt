"""Regression tests for fake-octree mesh-sampling locality and perf behavior.

These tests are intended to exercise the current mesh-sampling implementation in
`yt/geometry/particle_deposit.pyx`, specifically
`ParticleDepositOperation.process_octree()`. They validate correctness for
clustered/ordered/random/boundary particle access patterns and provide a local
benchmark for perf regressions.

These tests are marked `local` / `perf` so they can be run manually with:

    pytest -m local yt/geometry/tests/test_process_octree_locality.py
    pytest -m perf yt/geometry/tests/test_process_octree_locality.py

They are not excluded from the regular suite via a global pytest marker filter,
and they do not require any special global CI flag to keep them separate.
"""

import time

import numpy as np
import pytest
from numpy.testing import assert_allclose

from yt.testing import fake_octree_ds


def _create_fake_octree_dataset():
    """Return a small synthetic octree dataset for mesh-sampling tests."""
    return fake_octree_ds(num_zones=8, partial_coverage=0)


def _timeit(label, func, *args, **kwargs):
    """Measure wall-clock time for a callable.

    This helper is only used by the perf-marked timing regression below.
    """
    start = time.perf_counter()
    result = func(*args, **kwargs)
    elapsed = time.perf_counter() - start
    print(f"{label}: {elapsed:.6f}s")
    return result, elapsed


@pytest.mark.local
def test_fake_octree_mesh_sampling_clustered_points():
    """Clustered particle positions should still return correct mesh-sampled values."""
    ds = _create_fake_octree_dataset()
    # This test is intended to exercise the bulk octree lookup path inside
    # ParticleDepositOperation.process_octree() for points that are spatially
    # clustered and therefore likely to traverse only a small set of nearby cells.
    ds.add_mesh_sampling_particle_field(("gas", "density"), ptype="all")
    ad = ds.all_data()

    positions = ds.r["all", "particle_position"].to_value("code_length")
    center = positions[0]
    d = np.linalg.norm(positions - center, axis=1)
    cluster_idx = np.argsort(d)[:64]

    actual = ad["all", "cell_gas_density"][cluster_idx].to_value("code_density")
    baseline = ds.find_field_values_at_points(("gas", "density"), positions[cluster_idx]).to_value(
        "code_density"
    )

    assert_allclose(actual, baseline)


@pytest.mark.local
def test_fake_octree_mesh_sampling_ordered_and_random_positions():
    """Ordered and random position sets should match point-sampled reference values."""
    ds = _create_fake_octree_dataset()
    # Compare ordered particle lookup against random access to ensure the updated
    # octree implementation behaves correctly for both locality-friendly and
    # locality-unfriendly access patterns.
    ds.add_mesh_sampling_particle_field(("gas", "density"), ptype="all")
    ad = ds.all_data()

    positions = ds.r["all", "particle_position"].to_value("code_length")
    sorted_idx = np.lexsort((positions[:, 2], positions[:, 1], positions[:, 0]))[:128]
    random_idx = np.random.RandomState(1).choice(len(positions), size=128, replace=False)

    actual_sorted = ad["all", "cell_gas_density"][sorted_idx].to_value("code_density")
    baseline_sorted = ds.find_field_values_at_points(("gas", "density"), positions[sorted_idx]).to_value(
        "code_density"
    )
    actual_random = ad["all", "cell_gas_density"][random_idx].to_value("code_density")
    baseline_random = ds.find_field_values_at_points(("gas", "density"), positions[random_idx]).to_value(
        "code_density"
    )

    assert_allclose(actual_sorted, baseline_sorted)
    assert_allclose(actual_random, baseline_random)


@pytest.mark.local
def test_fake_octree_mesh_sampling_boundary_positions():
    """Boundary particles should still sample the correct octree cells."""
    ds = _create_fake_octree_dataset()
    # This test exercises the edge handling logic in process_octree() by using
    # particles near domain boundaries where octree cell membership can be subtle.
    ds.add_mesh_sampling_particle_field(("gas", "density"), ptype="all")
    ad = ds.all_data()

    positions = ds.r["all", "particle_position"].to_value("code_length")
    mask = np.any((positions < 0.1) | (positions > 0.9), axis=1)
    boundary_idx = np.nonzero(mask)[0][:64]
    if boundary_idx.size == 0:
        boundary_idx = np.arange(min(64, len(positions)))

    actual = ad["all", "cell_gas_density"][boundary_idx].to_value("code_density")
    baseline = ds.find_field_values_at_points(("gas", "density"), positions[boundary_idx]).to_value(
        "code_density"
    )

    assert_allclose(actual, baseline)


@pytest.mark.perf
def test_fake_octree_mesh_sampling_timing():
    """Measure fake-octree mesh sampling performance for a small dataset."""
    ds = _create_fake_octree_dataset()
    # This timing test is intended for local benchmarking of the updated octree
    # lookup path, not for strict CI timing thresholds.
    ds.add_mesh_sampling_particle_field(("gas", "density"), ptype="all")
    ad = ds.all_data()

    positions = ds.r["all", "particle_position"]
    actual, elapsed = _timeit("fake_octree_mesh_sampling", lambda: ad["all", "cell_gas_density"])
    assert actual.shape[0] == positions.shape[0]
    assert elapsed >= 0.0
