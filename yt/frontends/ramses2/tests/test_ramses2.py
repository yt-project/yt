"""Tests for the RAMSES II frontend, using small synthetic outputs."""

import numpy as np
import pytest

from yt.frontends.ramses2.api import RAMSES2Dataset
from yt.frontends.ramses2.data_structures import RAMSES2FileSanitizer
from yt.loaders import load
from yt.testing import assert_allclose_units

LEVELMIN = 2
LEVELMAX = 3
NVAR = 5
DM_FIELDS = "pos vel mass level birth_id"
STAR_FIELDS = "pos vel mass metallicity birth_date level birth_id"

# Octs stored in each file, as {level: [(ckey, refinement bitmask), ...]}.
# The root oct (0, 0, 0) is stored in the first file and has its first cell
# refined, but the corresponding child oct is stored in the second file.
ROOT_OCTS = [(i & 1, (i >> 1) & 1, (i >> 2) & 1) for i in range(8)]
FILES = [
    {2: [(ROOT_OCTS[0], 1)] + [(ck, 0) for ck in ROOT_OCTS[1:4]], 3: []},
    {2: [(ck, 0) for ck in ROOT_OCTS[4:]], 3: [((0, 0, 0), 0)]},
]


def _write_info(folder, nfile, aexp=1.0, H0=1.0, time=0.5):
    params = [
        ("nfile", nfile),
        ("ncpu", nfile),
        ("ndim", 3),
        ("levelmin", LEVELMIN),
        ("levelmax", LEVELMAX),
        ("ngridmax", 1000),
        ("nstep_coarse", 10),
    ]
    floats = [
        ("boxlen", 2.0),
        ("time", time),
        ("texp", time),
        ("aexp", aexp),
        ("H0", H0),
        ("omega_m", 1.0 if H0 == 1.0 else 0.3),
        ("omega_l", 0.0 if H0 == 1.0 else 0.7),
        ("omega_k", 0.0),
        ("omega_b", 0.0),
        ("gamma", 1.4),
        ("unit_l", 3.0856776e21),
        ("unit_d", 1.6726e-24),
        ("unit_t", 3.1557e13),
    ]
    with open(folder / "info.txt", "w") as f:
        for key, val in params:
            f.write(f"{key:<12s}={val:11d}\n")
        f.write("\n")
        for key, val in floats:
            f.write(f"{key:<12s}={val:23.15E}\n")
        f.write("\n")


def _cell_value(ckey, level, icell, ivar):
    """A unique value for each cell and variable."""
    x, y, z = ckey
    return 1.0 + ivar + 0.1 * level + 0.01 * icell + 1e-3 * (x + 2 * y + 4 * z)


def _write_particles(
    folder, prefix, fields, pos_per_file, id_dtype="<i4", legacy_header=False
):
    names = fields.split()
    with open(folder / f"{prefix}_header.txt", "w") as f:
        f.write(" Total number of particles\n")
        f.write(f"{sum(len(p) for p in pos_per_file):21d}\n")
        f.write(" Total number of files\n")
        f.write(f"{len(pos_per_file):12d}\n")
        f.write(" Particle fields\n")
        if legacy_header:
            f.write(fields + " ")
        else:
            for name in names:
                dtype = {"level": "i4", "birth_id": id_dtype[1:]}.get(name, "f4")
                ncomp = 3 if name in ("pos", "vel") else 1
                f.write(f" {name} {dtype} {ncomp}\n")

    first_id = 1
    for ifile, pos in enumerate(pos_per_file, start=1):
        pos = np.asarray(pos, dtype="float64").reshape(-1, 3)
        npart = len(pos)
        with open(folder / f"{prefix}.{ifile:05d}", "wb") as f:
            f.write(np.array([3, npart], dtype="<i4").tobytes())
            for name in names:
                if name == "pos":
                    data = [pos[:, i] for i in range(3)]
                elif name == "vel":
                    data = [np.full(npart, 0.1 * (i + 1)) for i in range(3)]
                elif name == "mass":
                    data = [np.full(npart, 0.5)]
                elif name in ("metallicity", "birth_date"):
                    data = [np.full(npart, 0.25)]
                if name == "level":
                    data = [np.full(npart, LEVELMIN, dtype="<i4")]
                elif name == "birth_id":
                    data = [np.arange(first_id, first_id + npart).astype(id_dtype)]
                else:
                    data = [d.astype("<f4") for d in data]
                for d in data:
                    f.write(d.tobytes())
        first_id += npart


@pytest.fixture
def ramses2_output(tmp_path):
    """Build a 2-file output and return its folder and the expected leaf cells."""
    return _make_output(tmp_path)


def _make_output(
    tmp_path,
    logical_refinement=False,
    id_dtype="<i4",
    legacy_header=False,
    **info_kwa,
):
    folder = tmp_path / "output_00001"
    folder.mkdir()
    nlevels = LEVELMAX - LEVELMIN + 1
    _write_info(folder, len(FILES), **info_kwa)

    with open(folder / "hydro_header.txt", "w") as f:
        f.write(f"nvar        ={NVAR:11d}\n")
        for i, name in enumerate(
            ["density", "velocity_x", "velocity_y", "velocity_z", "thermal_pressure"]
        ):
            f.write(f"variable #{i + 1:2d}: {name}\n")

    leaves = []
    for ifile, octs in enumerate(FILES, start=1):
        noct = np.array([len(octs[lvl]) for lvl in range(LEVELMIN, LEVELMAX + 1)])
        with (
            open(folder / f"amr.{ifile:05d}", "wb") as famr,
            open(folder / f"hydro.{ifile:05d}", "wb") as fhydro,
        ):
            famr.write(np.array([3, LEVELMIN, LEVELMAX], dtype="<i4").tobytes())
            famr.write(noct.astype("<i4").tobytes())
            fhydro.write(np.array([3, NVAR, LEVELMIN, LEVELMAX], dtype="<i4").tobytes())
            fhydro.write(noct.astype("<i4").tobytes())
            for level in range(LEVELMIN, LEVELMAX + 1):
                dx = 1.0 / 2**level
                for ckey, refined in octs[level]:
                    famr.write(np.array(ckey, dtype="<i4").tobytes())
                    if logical_refinement:
                        bits = [(refined >> i) & 1 for i in range(8)]
                        famr.write(np.array(bits, dtype="<i4").tobytes())
                    else:
                        famr.write(np.array([refined], dtype="<i4").tobytes())
                    block = np.array(
                        [
                            [
                                _cell_value(ckey, level, icell, ivar)
                                for icell in range(8)
                            ]
                            for ivar in range(NVAR)
                        ],
                        dtype="<f4",
                    )
                    fhydro.write(block.tobytes())
                    for icell in range(8):
                        if (refined >> icell) & 1:
                            continue
                        offset = np.array([(icell >> i) & 1 for i in range(3)])
                        center = (2 * np.array(ckey) + offset + 0.5) * dx
                        leaves.append((*center, block[0, icell]))
        assert len(noct) == nlevels

    # DM particles are spread over both files, stars only in the second one
    _write_particles(
        folder,
        "part",
        DM_FIELDS,
        [[(0.2, 0.2, 0.2), (1.5, 0.5, 0.5)], [(1.8, 1.8, 1.8)]],
        id_dtype=id_dtype,
        legacy_header=legacy_header,
    )
    _write_particles(
        folder,
        "star",
        STAR_FIELDS,
        [[], [(1.0, 1.0, 1.0)]],
        legacy_header=legacy_header,
    )

    return folder, np.array(leaves)


def test_sanitizer(ramses2_output):
    folder, _ = ramses2_output
    for fname in (folder, folder / "info.txt", folder / "amr.00002"):
        sanitizer = RAMSES2FileSanitizer(fname)
        assert sanitizer.is_valid
        assert sanitizer.root_folder == folder
    assert not RAMSES2FileSanitizer(folder.parent).is_valid
    assert RAMSES2Dataset._is_valid(str(folder))


def test_load(ramses2_output):
    folder, _ = ramses2_output
    ds = load(folder)
    assert isinstance(ds, RAMSES2Dataset)
    assert ds.nfile == 2
    assert ds.dimensionality == 3
    assert ds.min_level == LEVELMIN - 1
    assert ds.max_level == LEVELMAX - LEVELMIN
    np.testing.assert_equal(ds.domain_dimensions, 2**LEVELMIN)
    assert not ds.cosmological_simulation
    assert_allclose_units(ds.length_unit, ds.quan(2 * 3.0856776e21, "cm"))
    assert_allclose_units(ds.current_time, ds.quan(0.5 * 3.1557e13, "s"))
    assert set(ds.particle_types_raw) == {"part", "star"}
    assert ("ramses", "density") in ds.field_list
    assert ("star", "particle_birth_time") in ds.field_list


@pytest.mark.parametrize("logical_refinement", [False, True])
def test_cells(tmp_path, logical_refinement):
    folder, leaves = _make_output(tmp_path, logical_refinement=logical_refinement)
    ds = load(folder)
    ad = ds.all_data()

    # The refined cell of the first file is covered by an oct of the second
    # file: it must not be counted as a leaf.
    assert ad["index", "ones"].size == len(leaves) == 8 * 8 - 1 + 8
    assert_allclose_units(
        ad["index", "cell_volume"].sum(), ds.quan(1, "code_length**3")
    )

    xyz = np.stack([ad["index", ax].d for ax in "xyz"], axis=1)
    order_yt = np.lexsort(xyz.T)
    order_ref = np.lexsort(leaves[:, :3].T)
    np.testing.assert_allclose(xyz[order_yt], leaves[order_ref, :3])
    np.testing.assert_allclose(
        ad["ramses", "density"].d[order_yt], leaves[order_ref, 3], rtol=1e-6
    )
    assert ad["index", "grid_level"].max() == 1


def test_region_selection(ramses2_output):
    folder, leaves = ramses2_output
    ds = load(folder)
    reg = ds.region([0.5] * 3, [0, 0, 0], [0.5] * 3)
    centers = leaves[:, :3]
    expected = np.all(centers < 0.5, axis=1).sum()
    assert reg["ramses", "density"].size == expected
    # Only the root octs with x < 0.5 are selected
    assert reg["part", "particle_mass"].size == 1


def test_particles(ramses2_output):
    folder, _ = ramses2_output
    ds = load(folder)
    assert ds.particle_type_counts == {"part": 3, "star": 1}

    ad = ds.all_data()
    np.testing.assert_allclose(
        np.sort(ad["part", "particle_position_x"].to("code_length").d),
        [0.1, 0.75, 0.9],
        rtol=1e-6,
    )
    np.testing.assert_equal(np.sort(ad["part", "particle_identity"].d), [1, 2, 3])
    np.testing.assert_allclose(
        ad["part", "particle_velocity_z"].to("code_velocity").d, 0.3, rtol=1e-6
    )
    assert ad["star", "particle_mass"].size == 1
    assert_allclose_units(
        ad["star", "age"], ds.current_time - ds.quan(0.25, "code_time"), rtol=1e-6
    )


@pytest.mark.parametrize("legacy_header", [False, True])
def test_particle_long_ids(tmp_path, legacy_header):
    folder, _ = _make_output(tmp_path, id_dtype="<i8", legacy_header=legacy_header)
    ds = load(folder)
    ad = ds.all_data()
    np.testing.assert_equal(np.sort(ad["part", "particle_identity"].d), [1, 2, 3])
    np.testing.assert_allclose(ad["part", "particle_mass"].to("code_mass").d, 0.5)


def test_cosmology(tmp_path):
    folder, _ = _make_output(tmp_path, aexp=0.5, H0=70.0, time=-1.0)
    ds = load(folder)
    assert ds.cosmological_simulation
    assert ds.current_redshift == pytest.approx(1.0)
    assert ds.hubble_constant == pytest.approx(0.7)
    assert ds.current_time.to("Gyr") == pytest.approx(5.75, rel=0.02)


def test_empty_file(tmp_path):
    # A file can hold no octs at all, e.g. when nfile is large
    folder, leaves = _make_output(tmp_path)
    nlevels = LEVELMAX - LEVELMIN + 1
    header = np.array([3, LEVELMIN, LEVELMAX] + [0] * nlevels, dtype="<i4")
    (folder / "amr.00003").write_bytes(header.tobytes())
    header = np.array([3, NVAR, LEVELMIN, LEVELMAX] + [0] * nlevels, dtype="<i4")
    (folder / "hydro.00003").write_bytes(header.tobytes())
    info = (folder / "info.txt").read_text()
    (folder / "info.txt").write_text(
        info.replace("nfile       =          2", "nfile       =          3")
    )

    ds = load(folder)
    assert ds.nfile == 3
    assert ds.all_data()["ramses", "density"].size == len(leaves)


def test_legacy_header_three_fields(tmp_path):
    # A legacy header listing three fields looks like a new-style line
    folder, _ = _make_output(tmp_path)
    _write_particles(
        folder,
        "trac",
        "pos level birth_id",
        [[(0.2, 0.2, 0.2)], []],
        legacy_header=True,
    )
    ds = load(folder)
    ad = ds.all_data()
    np.testing.assert_equal(ad["trac", "particle_identity"].d, [1])
    np.testing.assert_equal(ad["trac", "particle_level"].d, [LEVELMIN])
