"""Tests for the RAMSES II frontend, using small synthetic outputs."""

import logging
from types import SimpleNamespace

import numpy as np
import pytest

from yt.frontends.ramses2.api import IOHandlerRAMSES2, RAMSES2Dataset
from yt.frontends.ramses2.data_structures import RAMSES2FileSanitizer
from yt.loaders import load
from yt.testing import assert_allclose_units
from yt.utilities.logger import ytLogger

LEVELMIN = 3
LEVELMAX = 6
NFILE = 3
BOXLEN = 2.0
# Probability for a cell coarser than LEVELMAX to be refined
REFINE_PROBABILITY = 0.15
HYDRO_FIELDS = ["density", "velocity_x", "velocity_y", "velocity_z", "thermal_pressure"]
NVAR = len(HYDRO_FIELDS)
DM_FIELDS = "pos vel mass level birth_id"
STAR_FIELDS = "pos vel mass metallicity birth_date level birth_id"
SINK_FIELDS = (
    "pos vel mass angmom birth_date internal_energy level accretion_step birth_id"
)

# Offsets of the 8 cells of an oct, in the RAMSES ordering (x varies fastest)
CELL_OFFSETS = np.array([[(icell >> i) & 1 for i in range(3)] for icell in range(8)])


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
        ("boxlen", BOXLEN),
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


def _morton_key(ckey, level):
    """Position of octs along the Morton curve, at the resolution of LEVELMAX."""
    key = np.zeros(len(ckey), dtype="int64")
    for bit in range(level - 1):
        for dim in range(3):
            key |= ((ckey[:, dim] >> bit) & 1) << (3 * bit + dim)
    return key << 3 * (LEVELMAX - level)


def _make_octree():
    """Build a random octree and split it into NFILE files.

    Return, for each level, the coordinates, Morton keys, refinement maps
    and files of the octs. Octs are sorted along the Morton curve, which is
    split into NFILE segments holding the same number of octs, in the same
    way as RAMSES II splits the Hilbert curve. Octs that straddle the
    boundary of two segments have some of their children in the next file.
    """
    rng = np.random.default_rng(5611)
    n = 2 ** (LEVELMIN - 1)
    ckey = np.indices((n, n, n)).reshape(3, -1).T
    octs = []
    for level in range(LEVELMIN, LEVELMAX + 1):
        key = _morton_key(ckey, level)
        order = np.argsort(key)
        ckey, key = ckey[order], key[order]
        if level < LEVELMAX:
            refined = rng.random((len(ckey), 8)) < REFINE_PROBABILITY
        else:
            refined = np.zeros((len(ckey), 8), dtype=bool)
        octs.append({"level": level, "ckey": ckey, "key": key, "refined": refined})
        iparent, icell = np.nonzero(refined)
        ckey = 2 * ckey[iparent] + CELL_OFFSETS[icell]

    keys = np.sort(np.concatenate([o["key"] for o in octs]))
    bounds = keys[len(keys) * np.arange(1, NFILE) // NFILE]
    for o in octs:
        o["file"] = 1 + np.searchsorted(bounds, o["key"], side="right")
    return octs


def _pick_particles(leaves, rng, per_level, where=None):
    """Place per_level particles at each level in distinct random leaf cells
    (among those selected by where), and return their positions (in code
    units), levels and files."""
    candidates = np.arange(len(leaves["level"]))
    if where is not None:
        candidates = candidates[where]
    icell = np.concatenate(
        [
            rng.choice(
                candidates[leaves["level"][candidates] == level],
                per_level,
                replace=False,
            )
            for level in range(LEVELMIN, LEVELMAX + 1)
        ]
    )
    dx = 0.5 ** leaves["level"][icell]
    pos = leaves["pos"][icell] + rng.uniform(-0.25, 0.25, (icell.size, 3)) * dx[:, None]
    # Positions are stored as float32
    pos = (pos * BOXLEN).astype("<f4").astype("float64") / BOXLEN
    return {"pos": pos, "level": leaves["level"][icell], "file": leaves["file"][icell]}


def _write_particles(
    folder,
    prefix,
    fields,
    particles,
    id_dtype="<i4",
    legacy_header=False,
    birth_date=0.25,
):
    """Write the files of a particle type. The id of each particle is its
    index (starting at 1) in particles."""
    names = fields.split()
    dtypes = {
        name: {
            "level": "<i4",
            "accretion_step": "<i4",
            "birth_id": id_dtype,
            "id": id_dtype,
        }.get(name, "<f4")
        for name in names
    }
    with open(folder / f"{prefix}_header.txt", "w") as f:
        f.write(" Total number of particles\n")
        f.write(f"{len(particles['level']):21d}\n")
        f.write(" Total number of files\n")
        f.write(f"{NFILE:12d}\n")
        f.write(" Particle fields\n")
        if legacy_header:
            f.write(fields + " ")
        else:
            for name in names:
                ncomp = 3 if name in ("pos", "vel", "angmom") else 1
                f.write(f" {name} {dtypes[name][1:]} {ncomp}\n")

    ids = np.arange(1, len(particles["level"]) + 1)
    for ifile in range(1, NFILE + 1):
        in_file = particles["file"] == ifile
        npart = in_file.sum()
        with open(folder / f"{prefix}.{ifile:05d}", "wb") as f:
            f.write(np.array([3, npart], dtype="<i4").tobytes())
            for name in names:
                if name == "pos":
                    data = (particles["pos"][in_file] * BOXLEN).T
                elif name in ("vel", "angmom"):
                    data = [np.full(npart, 0.1 * (i + 1)) for i in range(3)]
                elif name == "mass":
                    data = [np.full(npart, 0.5)]
                elif name in ("metallicity", "internal_energy"):
                    data = [np.full(npart, 0.25)]
                elif name == "birth_date":
                    data = [np.full(npart, birth_date)]
                elif name == "level":
                    data = [particles["level"][in_file]]
                elif name == "accretion_step":
                    data = [np.full(npart, 10)]
                elif name in ("birth_id", "id"):
                    data = [ids[in_file]]
                for d in data:
                    f.write(np.asarray(d).astype(dtypes[name]).tobytes())


@pytest.fixture
def ramses2_output(tmp_path):
    return _make_output(tmp_path)


def _make_output(
    tmp_path,
    logical_refinement=False,
    id_dtype="<i4",
    legacy_header=False,
    hydro_fields=HYDRO_FIELDS,
    **info_kwa,
):
    """Write a synthetic output, and return its folder, its octs, the
    expected leaf cells and its particles."""
    folder = tmp_path / "output_00001"
    folder.mkdir()
    nvar = len(hydro_fields)
    _write_info(folder, NFILE, **info_kwa)

    with open(folder / "hydro_header.txt", "w") as f:
        f.write(f"nvar        ={nvar:11d}\n")
        for i, name in enumerate(hydro_fields):
            f.write(f"variable #{i + 1:2d}: {name}\n")

    octs = _make_octree()
    # Fluid values are unique for each cell and variable
    first = 1
    for o in octs:
        size = len(o["ckey"]) * nvar * 8
        o["values"] = np.arange(first, first + size).astype("<f4").reshape(-1, nvar, 8)
        first += size

    nref = 8 if logical_refinement else 1
    record = np.dtype([("ckey", "<i4", (3,)), ("refined", "<i4", (nref,))])
    for ifile in range(1, NFILE + 1):
        in_file = [o["file"] == ifile for o in octs]
        noct = np.array([mask.sum() for mask in in_file], dtype="<i4")
        with (
            open(folder / f"amr.{ifile:05d}", "wb") as famr,
            open(folder / f"hydro.{ifile:05d}", "wb") as fhydro,
        ):
            famr.write(np.array([3, LEVELMIN, LEVELMAX], dtype="<i4").tobytes())
            famr.write(noct.tobytes())
            fhydro.write(np.array([3, nvar, LEVELMIN, LEVELMAX], dtype="<i4").tobytes())
            fhydro.write(noct.tobytes())
            for o, mask in zip(octs, in_file, strict=True):
                amr = np.empty(mask.sum(), dtype=record)
                amr["ckey"] = o["ckey"][mask]
                refined = o["refined"][mask]
                if logical_refinement:
                    amr["refined"] = refined
                else:
                    amr["refined"] = (refined << np.arange(8)).sum(
                        axis=1, keepdims=True
                    )
                famr.write(amr.tobytes())
                fhydro.write(o["values"][mask].tobytes())

    leaves = {"pos": [], "level": [], "file": [], "values": []}
    for o in octs:
        ioct, icell = np.nonzero(~o["refined"])
        dx = 0.5 ** o["level"]
        leaves["pos"].append((2 * o["ckey"][ioct] + CELL_OFFSETS[icell] + 0.5) * dx)
        leaves["level"].append(np.full(ioct.size, o["level"]))
        leaves["file"].append(o["file"][ioct])
        leaves["values"].append(o["values"][ioct, :, icell])
    leaves = {k: np.concatenate(v) for k, v in leaves.items()}

    # DM particles are spread over all files, stars are not in the first one
    rng = np.random.default_rng(0)
    particles = {
        "part": _pick_particles(leaves, rng, 4),
        "star": _pick_particles(leaves, rng, 1, where=leaves["file"] > 1),
    }
    _write_particles(
        folder,
        "part",
        DM_FIELDS,
        particles["part"],
        id_dtype=id_dtype,
        legacy_header=legacy_header,
    )
    _write_particles(
        folder, "star", STAR_FIELDS, particles["star"], legacy_header=legacy_header
    )

    return SimpleNamespace(folder=folder, octs=octs, leaves=leaves, particles=particles)


def test_synthetic_output(ramses2_output):
    # The random octree covers the cases the frontend has to handle
    octs = ramses2_output.octs
    assert all(len(o["ckey"]) > 0 for o in octs)
    # Cells are refined at every position within their oct
    assert np.concatenate([o["refined"] for o in octs]).any(axis=0).all()
    # Some octs are stored in another file than their parent
    cross_file = 0
    for parent, child in zip(octs[:-1], octs[1:], strict=True):
        shift = 3 * (LEVELMAX - parent["level"])
        iparent = np.searchsorted(parent["key"], child["key"] >> shift << shift)
        np.testing.assert_equal(parent["key"][iparent], child["key"] >> shift << shift)
        cross_file += np.count_nonzero(parent["file"][iparent] != child["file"])
    assert cross_file > 0
    assert np.unique(ramses2_output.particles["part"]["file"]).size == NFILE


def test_sanitizer(ramses2_output):
    folder = ramses2_output.folder
    for fname in (folder, folder / "info.txt", folder / "amr.00002"):
        sanitizer = RAMSES2FileSanitizer(fname)
        assert sanitizer.is_valid
        assert sanitizer.root_folder == folder
    assert not RAMSES2FileSanitizer(folder.parent).is_valid
    assert RAMSES2Dataset._is_valid(str(folder))


def test_load(ramses2_output):
    ds = load(ramses2_output.folder)
    assert isinstance(ds, RAMSES2Dataset)
    assert ds.nfile == NFILE
    assert ds.dimensionality == 3
    assert ds.min_level == LEVELMIN - 1
    assert ds.max_level == LEVELMAX - LEVELMIN
    assert ds.index.max_level == LEVELMAX - LEVELMIN
    assert ds.index.num_grids == sum(len(o["ckey"]) for o in ramses2_output.octs)
    np.testing.assert_equal(ds.domain_dimensions, 2**LEVELMIN)
    assert not ds.cosmological_simulation
    assert_allclose_units(ds.length_unit, ds.quan(BOXLEN * 3.0856776e21, "cm"))
    assert_allclose_units(ds.current_time, ds.quan(0.5 * 3.1557e13, "s"))
    assert set(ds.particle_types_raw) == {"part", "star"}
    assert ("ramses", "density") in ds.field_list
    assert ("star", "particle_birth_time") in ds.field_list


@pytest.mark.parametrize("logical_refinement", [False, True])
def test_cells(tmp_path, logical_refinement):
    out = _make_output(tmp_path, logical_refinement=logical_refinement)
    leaves = out.leaves
    ds = load(out.folder)
    ad = ds.all_data()

    # Refined cells, including those whose child oct is stored in another
    # file, must not be counted as leaves
    assert ad["index", "ones"].size == len(leaves["level"])
    assert_allclose_units(
        ad["index", "cell_volume"].sum(), ds.quan(1, "code_length**3")
    )

    xyz = np.stack([ad["index", ax].d for ax in "xyz"], axis=1)
    order_yt = np.lexsort(xyz.T)
    order_ref = np.lexsort(leaves["pos"].T)
    np.testing.assert_allclose(xyz[order_yt], leaves["pos"][order_ref])
    np.testing.assert_equal(
        ad["index", "grid_level"].d[order_yt], leaves["level"][order_ref] - LEVELMIN
    )
    for ivar, fname in enumerate(HYDRO_FIELDS):
        np.testing.assert_equal(
            ad["ramses", fname].d[order_yt], leaves["values"][order_ref, ivar]
        )


def test_region_selection(ramses2_output):
    out = ramses2_output
    ds = load(out.folder)
    left = np.array([0.1, 0.1, 0.05])
    right = np.array([0.6, 0.7, 0.6])
    reg = ds.region((left + right) / 2, left, right)

    def _inside(pos):
        return np.all((pos >= left) & (pos < right), axis=1)

    inside = _inside(out.leaves["pos"])
    assert np.unique(out.leaves["file"][inside]).size > 1
    np.testing.assert_equal(
        np.sort(reg["ramses", "density"].d), np.sort(out.leaves["values"][inside, 0])
    )

    for ptype, particles in out.particles.items():
        inside = _inside(particles["pos"])
        np.testing.assert_equal(
            np.sort(reg[ptype, "particle_identity"].d), np.flatnonzero(inside) + 1
        )
    inside = _inside(out.particles["part"]["pos"])
    assert 0 < inside.sum() < inside.size


def test_particles(ramses2_output):
    out = ramses2_output
    ds = load(out.folder)
    assert ds.particle_type_counts == {
        ptype: len(p["level"]) for ptype, p in out.particles.items()
    }

    ad = ds.all_data()
    for ptype, p in out.particles.items():
        order = np.argsort(ad[ptype, "particle_identity"].d)
        np.testing.assert_equal(
            ad[ptype, "particle_identity"].d[order], np.arange(1, len(p["level"]) + 1)
        )
        for i, ax in enumerate("xyz"):
            np.testing.assert_allclose(
                ad[ptype, f"particle_position_{ax}"].to("code_length").d[order],
                p["pos"][:, i],
            )
        np.testing.assert_equal(ad[ptype, "particle_level"].d[order], p["level"])
    np.testing.assert_allclose(
        ad["part", "particle_velocity_z"].to("code_velocity").d, 0.3, rtol=1e-6
    )
    assert_allclose_units(
        ad["star", "age"], ds.current_time - ds.quan(0.25, "code_time"), rtol=1e-6
    )


def test_all_data_particles(ramses2_output, monkeypatch):
    # Selecting all particles does not require reading their positions
    ds = load(ramses2_output.folder)
    fields_read = []
    read_particle_file = IOHandlerRAMSES2._read_particle_file

    def _read_particle_file(self, fn, ptype, fields):
        fields_read.extend(fields)
        return read_particle_file(self, fn, ptype, fields)

    monkeypatch.setattr(IOHandlerRAMSES2, "_read_particle_file", _read_particle_file)
    mass = ds.all_data()["part", "particle_mass"]
    assert mass.size == len(ramses2_output.particles["part"]["level"])
    assert set(fields_read) == {"particle_mass"}


@pytest.mark.parametrize("legacy_header", [False, True])
def test_particle_long_ids(tmp_path, legacy_header):
    out = _make_output(tmp_path, id_dtype="<i8", legacy_header=legacy_header)
    ds = load(out.folder)
    ad = ds.all_data()
    npart = len(out.particles["part"]["level"])
    np.testing.assert_equal(
        np.sort(ad["part", "particle_identity"].d), np.arange(1, npart + 1)
    )
    np.testing.assert_allclose(ad["part", "particle_mass"].to("code_mass").d, 0.5)


def test_cosmology(tmp_path):
    out = _make_output(tmp_path, aexp=0.5, H0=70.0, time=-1.0)
    # Birth dates are proper times since today, in units of unit_t / aexp**2,
    # like texp (-1)
    stars = out.particles["star"]
    _write_particles(out.folder, "star", STAR_FIELDS, stars, birth_date=-1.25)
    ds = load(out.folder)
    assert ds.cosmological_simulation
    assert ds.current_redshift == pytest.approx(1.0)
    assert ds.hubble_constant == pytest.approx(0.7)
    assert ds.current_time.to("Gyr") == pytest.approx(5.75, rel=0.02)
    ad = ds.all_data()
    assert ad["star", "age"].size == len(stars["level"])
    assert_allclose_units(
        ad["star", "age"], ds.quan((-1.0 + 1.25) / 0.5**2, "code_time"), rtol=1e-6
    )


def test_empty_file(tmp_path):
    # A file can hold no octs at all, e.g. when nfile is large
    out = _make_output(tmp_path)
    folder = out.folder
    nlevels = LEVELMAX - LEVELMIN + 1
    header = np.array([3, LEVELMIN, LEVELMAX] + [0] * nlevels, dtype="<i4")
    (folder / f"amr.{NFILE + 1:05d}").write_bytes(header.tobytes())
    header = np.array([3, NVAR, LEVELMIN, LEVELMAX] + [0] * nlevels, dtype="<i4")
    (folder / f"hydro.{NFILE + 1:05d}").write_bytes(header.tobytes())
    _write_info(folder, NFILE + 1)

    ds = load(folder)
    assert ds.nfile == NFILE + 1
    assert ds.all_data()["ramses", "density"].size == len(out.leaves["level"])


def test_legacy_header_three_fields(tmp_path):
    # A legacy header listing three fields looks like a new-style line
    out = _make_output(tmp_path)
    trac = _pick_particles(out.leaves, np.random.default_rng(1), 1)
    _write_particles(out.folder, "trac", "pos level birth_id", trac, legacy_header=True)
    ds = load(out.folder)
    ad = ds.all_data()
    order = np.argsort(ad["trac", "particle_identity"].d)
    np.testing.assert_equal(
        ad["trac", "particle_identity"].d[order], np.arange(1, len(trac["level"]) + 1)
    )
    np.testing.assert_equal(ad["trac", "particle_level"].d[order], trac["level"])


def test_mhd(tmp_path):
    bfields = [f"magnetic_field_{ax}" for ax in "xyz"]
    out = _make_output(tmp_path, hydro_fields=HYDRO_FIELDS + bfields)
    ds = load(out.folder)
    ad = ds.all_data()
    b2 = sum(ad["ramses", f].to("code_magnetic").d ** 2 for f in bfields)
    for f in bfields:
        assert_allclose_units(ad["gas", f], ad["ramses", f])
    np.testing.assert_allclose(
        ad["gas", "magnetic_field_strength"].to("code_magnetic").d, np.sqrt(b2)
    )
    # B**2 / 2 in code units, with the RAMSES (Lorentz-Heaviside) convention
    np.testing.assert_allclose(
        ad["gas", "magnetic_pressure"].to("code_pressure").d, b2 / 2, rtol=1e-6
    )
    assert ("gas", "alfven_speed") in ds.derived_field_list


def test_sink_fields(tmp_path):
    out = _make_output(tmp_path)
    center = np.array([0.25, 0.25, 0.25])
    near = np.linalg.norm(out.leaves["pos"] - center, axis=1) < 0.4
    sink = _pick_particles(out.leaves, np.random.default_rng(1), 1, where=near)
    _write_particles(out.folder, "sink", SINK_FIELDS, sink)
    ds = load(out.folder)
    ad = ds.all_data()
    # angmom is the specific angular momentum of the accreted gas, about the
    # sink, with lengths in unit_l, whereas code_length is unit_l * boxlen
    for i, ax in enumerate("xyz"):
        spin = ad["sink", f"particle_spin_{ax}"]
        assert spin.units.dimensions == ds.quan(1, "cm**2/s").units.dimensions
        np.testing.assert_allclose(
            spin.to("code_length * code_velocity").d,
            0.1 * (i + 1) / BOXLEN,
            rtol=1e-6,
        )
    # The spin does not replace the orbital angular momentum derived by yt
    # about the center of the data object
    sp = ds.sphere(center, (0.5, "code_length"))
    sp.set_field_parameter("bulk_velocity", ds.arr([0, 0, 0], "code_velocity"))
    assert sp["sink", "particle_mass"].size == len(sink["level"])
    r = (
        sp["sink", "particle_position"].to("code_length").d
        - sp.center.to("code_length").d
    )
    v = sp["sink", "particle_velocity"].to("code_velocity").d
    np.testing.assert_allclose(
        sp["sink", "particle_specific_angular_momentum"]
        .to("code_length * code_velocity")
        .d,
        np.cross(r, v),
        rtol=1e-6,
    )
    np.testing.assert_allclose(
        ad["sink", "particle_internal_energy"].to("code_velocity**2").d, 0.25
    )
    np.testing.assert_equal(ad["sink", "particle_accretion_step"].d, 10)


def test_legacy_header_unknown_field(tmp_path, caplog, monkeypatch):
    # The size of unknown fields is not known: their particles cannot be read
    out = _make_output(tmp_path)
    (out.folder / "trac_header.txt").write_text(
        " Total number of particles\n"
        f"{1:21d}\n"
        " Total number of files\n"
        f"{NFILE:12d}\n"
        " Particle fields\n"
        "pos level unknown_field \n"
    )
    monkeypatch.setattr(ytLogger, "propagate", True)
    with caplog.at_level(logging.WARNING, "yt"):
        ds = load(out.folder)
    assert set(ds.particle_types_raw) == {"part", "star"}
    assert "unknown fields ['unknown_field']" in caplog.text


def test_legacy_header_id(tmp_path):
    # Some old outputs name the particle ids "id" instead of "birth_id"
    out = _make_output(tmp_path)
    part = out.particles["part"]
    _write_particles(
        out.folder, "part", "pos vel mass level id", part, legacy_header=True
    )
    ds = load(out.folder)
    assert ds.particle_type_counts["part"] == len(part["level"])
    np.testing.assert_equal(
        np.sort(ds.all_data()["part", "particle_identity"].d),
        np.arange(1, len(part["level"]) + 1),
    )
