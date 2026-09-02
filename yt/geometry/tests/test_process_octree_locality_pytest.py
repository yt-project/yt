"""
Fake-octree mesh-sampling tests for
ParticleDepositOperation.process_octree() (yt/geometry/particle_deposit.pyx).

Run with:

    pytest yt/geometry/tests/test_process_octree_locality.py -s

To sweep different particle counts, set YT_PERF_N_PARTICLES (comma-separated)
before that same command -- it replaces the default sweep in
_DEFAULT_N_PARTICLES for both test_fake_octree_mesh_sampling and
test_deep_octree_locality:

    YT_PERF_N_PARTICLES=1000,100000,5000000 pytest yt/geometry/tests/test_process_octree_locality.py -s
"""

import contextlib
import functools
import io
import os
import time

import numpy as np
from numpy.random import RandomState
from numpy.testing import assert_allclose
from yt.testing import fake_octree_ds

_CORRECTNESS_SAMPLE = 200
_DEFAULT_N_PARTICLES = [128, 1_024, 8_192, 65_536, 262_144]


def _create_fake_octree_dataset():
    # mesh_sampling_particle_field() only handles the octree's native 2x2x2
    # branching, so num_zones must stay at 2.
    #prng is shared, so initialize it randomly every time
    return fake_octree_ds(num_zones=2, partial_coverage=0, prng=RandomState(0x4D3D3D3))


def _clustered(ds, n, seed):
    prng = np.random.RandomState(seed)
    center = prng.random_sample(3)
    pos = np.clip(center + prng.normal(scale=0.02, size=(n, 3)), 0.0, 1.0)
    return ds.arr(pos, "code_length")


def _random(ds, n, seed):
    prng = np.random.RandomState(seed)
    return ds.arr(prng.random_sample((n, 3)), "code_length")


def _ordered(ds, n, seed):
    positions = _random(ds, n, seed)
    pos = positions.to_value("code_length")
    idx = np.lexsort((pos[:, 2], pos[:, 1], pos[:, 0]))
    return positions[idx]


def _boundary(ds, n, seed):
    prng = np.random.RandomState(seed)
    pos = prng.random_sample((n, 3))
    on_low_edge = prng.randint(0, 2, size=n).astype(bool)
    pos[:, 0] = np.where(on_low_edge, pos[:, 0] * 0.01, 1.0 - pos[:, 0] * 0.01)
    return ds.arr(pos, "code_length")


_PATTERNS = {
    "clustered": _clustered,
    "ordered": _ordered,
    "random": _random,
    "boundary": _boundary,
}


def _reference_density_at(ds, positions):
    # Not using ds.find_field_values_at_points() -- it breaks under NumPy 2.x.
    pos = positions.to_value("code_length")
    values = [ds.point(p)["gas", "density"][0] for p in pos]
    return ds.arr(values, "code_density")


def _timeit(label, func, *args, **kwargs):
    start = time.perf_counter()
    result = func(*args, **kwargs)
    elapsed = time.perf_counter() - start
    print(f"{label}: {elapsed:.6f}s")
    return result, elapsed


def pytest_generate_tests(metafunc):
    if "pattern" in metafunc.fixturenames:
        metafunc.parametrize("pattern", list(_PATTERNS))
    if "n_particles" in metafunc.fixturenames:
        raw = os.environ.get("YT_PERF_N_PARTICLES")
        values = [int(v) for v in raw.split(",")] if raw else _DEFAULT_N_PARTICLES
        metafunc.parametrize("n_particles", values)


def test_fake_octree_mesh_sampling(n_particles, pattern):
    ds = _create_fake_octree_dataset()
    positions = _PATTERNS[pattern](ds, n_particles, seed=1)

    ad = ds.all_data()
    density = ad["gas", "density"]
    obj = ad._current_chunk.objs[0]
    mesh_field = np.asarray(density.T.reshape(-1))

    # Only the actual octree lookup is inside the clock -- setup above and
    # the correctness check below both happen outside it.
    sampled, elapsed = _timeit(
        f"fake_octree_mesh_sampling ({pattern}, n_particles={n_particles})",
        obj.mesh_sampling_particle_field,
        positions,
        mesh_field,
    )
    assert sampled.shape[0] == n_particles
    assert elapsed >= 0.0

    # The reference lookup is a slow per-point Python loop, so only check a
    # capped sample instead of all n_particles.
    check_n = min(n_particles, _CORRECTNESS_SAMPLE)
    actual = ds.arr(sampled[:check_n], density.units).to_value("code_density")
    baseline = _reference_density_at(ds, positions[:check_n]).to_value("code_density")
    assert_allclose(actual, baseline)


# The fixture above is a shallow octree (depth ~2) -- there's nothing for a
# same-oct cache or a climb-instead-of-restart lookup to save, no matter how
# good the implementation is. This section builds a genuinely deep octree
# directly (bypassing the fake_octree_ds/load_octree path, which only goes a
# couple levels deep for this test's mask) so the perf numbers below actually
# have tree depth to exploit.

_DEEP_MAX_LEVEL = 24
_DEEP_MAX_NOCT = 4000
_DEEP_FSUBDIVIDE = 0.05


@functools.lru_cache(maxsize=1)
def _build_deep_octree():
    from yt.geometry.fake_octree import create_fake_octree
    from yt.geometry.oct_container import RAMSESOctreeContainer

    # create_fake_octree() only ever populates a single root oct (at index
    # [0, 0, 0]); domain_dimensions must be [1, 1, 1] so that root's cell
    # spans the whole [0, 1)^3 domain. With [2, 2, 2] (8 root cells),
    # particles in the other 7 root cells would map to a NULL oct and get
    # silently dropped by process_octree.
    dd = np.array([1, 1, 1], dtype="i4")
    dle = np.array([0.0, 0.0, 0.0], dtype="f8")
    dre = np.array([1.0, 1.0, 1.0], dtype="f8")
    oct_handler = RAMSESOctreeContainer(dd, dle, dre)
    # create_fake_octree prints a line per oct it creates -- silence it.
    with contextlib.redirect_stdout(io.StringIO()):
        create_fake_octree(
            oct_handler, _DEEP_MAX_NOCT, _DEEP_MAX_LEVEL, dd, dle, dre, _DEEP_FSUBDIVIDE
        )
    return oct_handler


def _random_raw(n, seed):
    return RandomState(seed).random_sample((n, 3))


def _ordered_raw(n, seed):
    pos = _random_raw(n, seed)
    idx = np.lexsort((pos[:, 2], pos[:, 1], pos[:, 0]))
    return pos[idx]


def _clustered_raw(n, seed):
    prng = RandomState(seed)
    center = prng.random_sample(3)
    return np.clip(center + prng.normal(scale=0.02, size=(n, 3)), 0.0, 1.0)


def _boundary_raw(n, seed):
    prng = RandomState(seed)
    pos = prng.random_sample((n, 3))
    on_low_edge = prng.randint(0, 2, size=n).astype(bool)
    pos[:, 0] = np.where(on_low_edge, pos[:, 0] * 0.01, 1.0 - pos[:, 0] * 0.01)
    return pos


_RAW_PATTERNS = {
    "clustered": _clustered_raw,
    "ordered": _ordered_raw,
    "random": _random_raw,
    "boundary": _boundary_raw,
}


def test_deep_octree_locality(n_particles, pattern):
    from yt.geometry.particle_deposit import CountParticles

    oct_handler = _build_deep_octree()
    positions = _RAW_PATTERNS[pattern](n_particles, seed=1)

    dom_ind = np.arange(oct_handler.nocts, dtype=np.int64)
    # RAMSESOctreeContainer was built with the default num_zones=2 -- see
    # _build_deep_octree, which never overrides it.
    nz = (2, 2, 2, oct_handler.nocts)
    op = CountParticles(nz, "cubic")
    op.initialize()

    _, elapsed = _timeit(
        f"deep_octree_process_octree ({pattern}, n_particles={n_particles}, "
        f"max_level={_DEEP_MAX_LEVEL}, nocts={oct_handler.nocts})",
        op.process_octree,
        oct_handler,
        dom_ind,
        positions,
        None,
        1,
        0,
    )
    counted = op.finalize()
    assert counted.sum() == n_particles
    assert elapsed >= 0.0
