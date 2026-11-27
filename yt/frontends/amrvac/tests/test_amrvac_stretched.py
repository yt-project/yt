from types import SimpleNamespace

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from yt.frontends.amrvac.data_structures import (
    AMRVACHierarchy,
    _parse_qstretch_baselevel,
    _parse_stretch_dim,
    _stretched_cell_edges,
)
from yt.frontends.amrvac.io import read_amrvac_namelist
from yt.testing import requires_module

# --- parfile parsing -----------------------------------------------------------


@pytest.mark.parametrize(
    "raw, expected",
    [
        (None, ["none", "none"]),
        ("uni", ["uni", "none"]),  # scalar: first dimension (Fortran semantics)
        (["uni", ""], ["uni", "none"]),
        ([None, "uniform"], ["none", "uni"]),  # indexed assignment stretch_dim(2)
        (["UNI", "none", None], ["uni", "none"]),  # trailing unused entries are ok
    ],
)
def test_parse_stretch_dim(raw, expected):
    meshlist = {} if raw is None else {"stretch_dim": raw}
    assert _parse_stretch_dim(meshlist, ndim=2) == expected


@pytest.mark.parametrize(
    "raw, error",
    [
        ("symm", NotImplementedError),
        ("foo", ValueError),
        ([1.0], ValueError),
        (["none", "uni"], ValueError),  # stretched dimension beyond ndim
    ],
)
def test_parse_stretch_dim_errors(raw, error):
    with pytest.raises(error):
        _parse_stretch_dim({"stretch_dim": raw}, ndim=1)


@pytest.mark.parametrize(
    "raw, expected",
    [
        (None, [None, None]),
        (1.02, [1.02, None]),  # scalar: first dimension (Fortran semantics)
        (2, [2.0, None]),
        ([None, 1.05], [None, 1.05]),
        ([1.01, 1.02, 1.03], [1.01, 1.02]),
    ],
)
def test_parse_qstretch_baselevel(raw, expected):
    meshlist = {} if raw is None else {"qstretch_baselevel": raw}
    assert _parse_qstretch_baselevel(meshlist, ndim=2) == expected


def test_parse_qstretch_baselevel_errors():
    with pytest.raises(ValueError):
        _parse_qstretch_baselevel({"qstretch_baselevel": "1.02"}, ndim=1)


@requires_module("f90nml")
@pytest.mark.parametrize(
    "statement, expected",
    [
        ("stretch_dim(2) = 'uni'", ["none", "uni"]),
        ("stretch_dim = 'uni'", ["uni", "none"]),
        ("stretch_dim = 'uni', 'uni'", ["uni", "uni"]),
    ],
)
def test_parse_stretch_dim_from_parfile(tmp_path, statement, expected):
    parfile = tmp_path / "stretch.par"
    parfile.write_text(f"&meshlist\n  {statement}\n/\n")
    meshlist = read_amrvac_namelist(str(parfile))["meshlist"]
    assert _parse_stretch_dim(meshlist, ndim=2) == expected


# --- geometry ------------------------------------------------------------------


class _StubHierarchy(AMRVACHierarchy):
    # __init__ is bypassed, so there are no grid arrays for GridIndex.__del__ to drop
    def __del__(self):
        pass


def _make_index(stretch_dim, qstretch_baselevel, xmin, xmax, domain_nx, max_level):
    """A minimal AMRVACHierarchy, enough to compute cell edges."""
    index = _StubHierarchy.__new__(_StubHierarchy)
    index.dataset = SimpleNamespace(
        dimensionality=len(stretch_dim),
        refine_by=2,
        parameters={
            "xmin": np.asarray(xmin, dtype="float64"),
            "xmax": np.asarray(xmax, dtype="float64"),
            "domain_nx": np.asarray(domain_nx),
        },
    )
    index.stretch_dim = stretch_dim
    index._qstretch_baselevel = qstretch_baselevel
    index.is_stretched = any(s != "none" for s in stretch_dim)
    index.max_level = max_level
    index._setup_stretching()
    return index


def test_stretched_cell_edges_closed_form():
    xmin, dxfirst, q = 2.0, 0.1, 1.1
    k = np.arange(10)
    widths = np.diff(_stretched_cell_edges(xmin, dxfirst, q, k))
    assert_allclose(widths, dxfirst * q ** k[:-1])
    assert _stretched_cell_edges(xmin, dxfirst, q, 0) == xmin


@pytest.mark.parametrize("q", [None, 1.2, 1.0001])
def test_base_level_spans_domain(q):
    nx = 32
    index = _make_index(["uni"], [q], [1.0], [10.0], [nx], max_level=0)
    edges = index._cell_edges(0, 0, np.arange(nx + 1))
    assert_allclose(edges[[0, -1]], [1.0, 10.0])
    assert np.all(np.diff(edges) > 0)


def test_default_qstretch_is_logarithmic():
    nx = 16
    index = _make_index(["uni"], [None], [1.0], [100.0], [nx], max_level=0)
    edges = index._cell_edges(0, 0, np.arange(nx + 1))
    assert_allclose(edges, np.logspace(0, 2, nx + 1))


def test_refined_levels_nest_in_parents():
    """Every coarse edge is also an edge on finer levels (2 children per cell)."""
    nx = 8
    index = _make_index(["uni"], [1.1], [0.0], [5.0], [nx], max_level=3)
    coarse = index._cell_edges(0, 0, np.arange(nx + 1))
    for ytlevel in range(1, 4):
        fine = index._cell_edges(0, ytlevel, np.arange(nx * 2**ytlevel + 1))
        assert_allclose(fine[:: 2**ytlevel], coarse, rtol=1e-12)


def _make_2d_index():
    # x uniform on [0, 2], y stretched on [1, 3], 4x4 cells, two levels
    return _make_index(
        ["none", "uni"], [None, 1.05], [0.0, 1.0], [2.0, 3.0], [4, 4], max_level=1
    )


def test_unstretched_dimension_is_uniform():
    index = _make_2d_index()
    assert_allclose(index._cell_edges(0, 1, np.arange(9)), np.linspace(0, 2, 9))


def test_icoords_to_fcoords_matches_edges():
    index = _make_2d_index()
    index.dataset.domain_left_edge = SimpleNamespace(d=np.array([0.0, 1.0, 0.0]))
    index.dataset.domain_right_edge = SimpleNamespace(d=np.array([2.0, 3.0, 1.0]))
    icoords = np.array([[0, 0, 0], [3, 5, 0], [1, 7, 0]])
    ires = np.array([0, 1, 1])
    coords, widths = index._icoords_to_fcoords(icoords, ires)
    for (i, j, _), lev, c, w in zip(icoords, ires, coords, widths, strict=True):
        y_edges = index._cell_edges(1, lev, np.array([j, j + 1]))
        assert_allclose(c[1], y_edges.mean())
        assert_allclose(w[1], np.diff(y_edges)[0])
        assert_allclose(w[0], 2.0 / (4 * 2**lev))
        assert_allclose(c[0], (i + 0.5) * w[0])
    assert_array_equal(coords[:, 2], 0.5)
    assert_array_equal(widths[:, 2], 1.0)


def test_qstretch_errors():
    with pytest.raises(ValueError, match="must be set"):
        _make_index(["uni"], [None], [0.0], [1.0], [8], max_level=0)
    with pytest.raises(ValueError, match="no stretching"):
        _make_index(["uni"], [1.0], [0.0], [1.0], [8], max_level=0)
