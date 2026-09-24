"""Tests for the mini-ramses frontend."""

import os
import struct
import tempfile

import numpy as np
import pytest

from yt.frontends.mini_ramses.data_structures import (
    MiniRAMSESDataset,
    MiniRAMSESFileSanitizer,
)


def _create_mini_ramses_output(
    tmpdir, ndim=3, levelmin=3, nlevelmax=5, nvar=5, npart=10
):
    """Create a minimal synthetic mini-ramses output directory."""
    outdir = os.path.join(tmpdir, "output_00001")
    os.makedirs(outdir, exist_ok=True)

    boxlen = 1.0
    gamma = 1.4

    # Write info.txt (mini-ramses format: nfile first!)
    with open(os.path.join(outdir, "info.txt"), "w") as f:
        f.write(f"nfile       ={'1':>11s}\n")
        f.write(f"ncpu        ={'1':>11s}\n")
        f.write(f"ndim        ={ndim:>11d}\n")
        f.write(f"levelmin    ={levelmin:>11d}\n")
        f.write(f"levelmax    ={nlevelmax:>11d}\n")
        f.write(f"ngridmax    ={'100':>11s}\n")
        f.write(f"nstep_coarse={'0':>11s}\n")
        f.write("\n")
        f.write(f"boxlen      ={boxlen:>23.15E}\n")
        f.write(f"time        ={0.0:>23.15E}\n")
        f.write(f"texp        ={0.0:>23.15E}\n")
        f.write(f"aexp        ={1.0:>23.15E}\n")
        f.write(f"H0          ={1.0:>23.15E}\n")
        f.write(f"omega_m     ={0.0:>23.15E}\n")
        f.write(f"omega_l     ={0.0:>23.15E}\n")
        f.write(f"omega_k     ={0.0:>23.15E}\n")
        f.write(f"omega_b     ={0.0:>23.15E}\n")
        f.write(f"gamma       ={gamma:>23.15E}\n")
        f.write(f"unit_l      ={3.08568e+21:>23.15E}\n")
        f.write(f"unit_d      ={1.6726e-24:>23.15E}\n")
        f.write(f"unit_t      ={3.1557e+13:>23.15E}\n")
        f.write("\n")

    # Write amr.00001 (stream binary)
    twotondim = 2**ndim
    # Place a few octs at levelmin
    nocts_per_level = np.zeros(nlevelmax, dtype="int32")
    n_base_octs = 2  # number of octs at the base level
    nocts_per_level[levelmin - 1] = n_base_octs

    with open(os.path.join(outdir, "amr.00001"), "wb") as f:
        f.write(struct.pack("i", ndim))
        f.write(struct.pack("i", levelmin))
        f.write(struct.pack("i", nlevelmax))
        for ilevel in range(levelmin - 1, nlevelmax):
            f.write(struct.pack("i", nocts_per_level[ilevel]))

        # Write grid data
        for ilevel in range(levelmin - 1, nlevelmax):
            for igrid in range(nocts_per_level[ilevel]):
                # Cartesian key (ndim int32 values)
                ckey = [igrid] * ndim
                f.write(struct.pack(f"{ndim}i", *ckey))
                # Refinement map (twotondim int32 values) - not refined
                refined = [0] * twotondim
                f.write(struct.pack(f"{twotondim}i", *refined))

    # Write hydro.00001 (stream binary)
    with open(os.path.join(outdir, "hydro.00001"), "wb") as f:
        f.write(struct.pack("i", ndim))
        f.write(struct.pack("i", nvar))
        f.write(struct.pack("i", levelmin))
        f.write(struct.pack("i", nlevelmax))
        for ilevel in range(levelmin - 1, nlevelmax):
            f.write(struct.pack("i", nocts_per_level[ilevel]))

        # Write cell data per level
        for ilevel in range(levelmin - 1, nlevelmax):
            ncache = nocts_per_level[ilevel]
            for igrid in range(ncache):
                # qout(twotondim, nvar) as float32
                qout = np.random.rand(twotondim, nvar).astype("float32")
                # Set density to be positive
                qout[:, 0] = np.abs(qout[:, 0]) + 0.1
                # Set pressure to be positive
                if nvar >= 5:
                    qout[:, 4] = np.abs(qout[:, 4]) + 0.01
                f.write(qout.tobytes())

    # Write hydro_header.txt
    with open(os.path.join(outdir, "hydro_header.txt"), "w") as f:
        f.write(f"nvar        ={nvar:>11d}\n")
        f.write("variable # 1: density\n")
        f.write("variable # 2: velocity_x\n")
        f.write("variable # 3: velocity_y\n")
        f.write("variable # 4: velocity_z\n")
        f.write("variable # 5: thermal_pressure\n")

    # Write part.00001 (stream binary) - particle data
    with open(os.path.join(outdir, "part.00001"), "wb") as f:
        f.write(struct.pack("i", ndim))
        f.write(struct.pack("i", npart))

        # Positions (float32)
        for ax in range(ndim):
            pos = np.random.rand(npart).astype("float32") * boxlen
            f.write(pos.tobytes())

        # Velocities (float32)
        for ax in range(ndim):
            vel = (np.random.rand(npart).astype("float32") - 0.5) * 100
            f.write(vel.tobytes())

        # Mass (float32)
        mass = np.random.rand(npart).astype("float32") * 1e-3
        f.write(mass.tobytes())

        # Level (int32)
        levels = np.ones(npart, dtype="int32") * levelmin
        f.write(levels.tobytes())

        # Birth ID (int32)
        ids = np.arange(1, npart + 1, dtype="int32")
        f.write(ids.tobytes())

    # Write part_header.txt
    with open(os.path.join(outdir, "part_header.txt"), "w") as f:
        f.write("Total number of particles\n")
        f.write(f"{npart}\n")
        f.write("Total number of files\n")
        f.write("1\n")
        f.write("Particle fields\n")
        f.write("pos vel mass level birth_id \n")

    return outdir


class TestMiniRAMSESFileSanitizer:
    def test_valid_info_file(self, tmp_path):
        outdir = _create_mini_ramses_output(str(tmp_path))
        info_path = os.path.join(outdir, "info.txt")
        sanitizer = MiniRAMSESFileSanitizer(info_path)
        assert sanitizer.is_valid

    def test_valid_directory(self, tmp_path):
        outdir = _create_mini_ramses_output(str(tmp_path))
        sanitizer = MiniRAMSESFileSanitizer(outdir)
        assert sanitizer.is_valid

    def test_valid_amr_file(self, tmp_path):
        outdir = _create_mini_ramses_output(str(tmp_path))
        amr_path = os.path.join(outdir, "amr.00001")
        sanitizer = MiniRAMSESFileSanitizer(amr_path)
        assert sanitizer.is_valid

    def test_invalid_file(self, tmp_path):
        sanitizer = MiniRAMSESFileSanitizer(str(tmp_path / "nonexistent"))
        assert not sanitizer.is_valid

    def test_ramses_info_not_valid(self, tmp_path):
        """Ensure RAMSES info files (ncpu first, not nfile) are rejected."""
        outdir = tmp_path / "output_00001"
        outdir.mkdir()
        info = outdir / "info.txt"
        info.write_text("ncpu        =          1\nndim        =          3\n")
        sanitizer = MiniRAMSESFileSanitizer(str(info))
        assert not sanitizer.is_valid


class TestMiniRAMSESDataset:
    def test_is_valid(self, tmp_path):
        outdir = _create_mini_ramses_output(str(tmp_path))
        assert MiniRAMSESDataset._is_valid(outdir)

    def test_is_valid_info_file(self, tmp_path):
        outdir = _create_mini_ramses_output(str(tmp_path))
        info_path = os.path.join(outdir, "info.txt")
        assert MiniRAMSESDataset._is_valid(info_path)

    def test_not_valid_for_ramses(self, tmp_path):
        """Ensure standard RAMSES outputs are NOT detected as mini-ramses."""
        outdir = tmp_path / "output_00001"
        outdir.mkdir()
        info = outdir / "info.txt"
        info.write_text(
            "ncpu        =          1\n"
            "ndim        =          3\n"
            "levelmin    =          3\n"
        )
        assert not MiniRAMSESDataset._is_valid(str(outdir))

    def test_load_dataset(self, tmp_path):
        outdir = _create_mini_ramses_output(str(tmp_path))
        ds = MiniRAMSESDataset(outdir)
        assert ds is not None
        assert ds.dimensionality == 3
        assert ds.min_level == 3
        assert ds.max_level == 5
        assert ds.gamma == 1.4
        assert ds.cosmological_simulation is False

    def test_domain_dimensions(self, tmp_path):
        outdir = _create_mini_ramses_output(str(tmp_path), levelmin=3)
        ds = MiniRAMSESDataset(outdir)
        # domain_dimensions should be 2^levelmin in each direction
        assert np.all(ds.domain_dimensions == 8)

    def test_unit_attributes(self, tmp_path):
        outdir = _create_mini_ramses_output(str(tmp_path))
        ds = MiniRAMSESDataset(outdir)
        assert ds.length_unit is not None
        assert ds.mass_unit is not None
        assert ds.time_unit is not None

    def test_fluid_fields_detected(self, tmp_path):
        outdir = _create_mini_ramses_output(str(tmp_path))
        ds = MiniRAMSESDataset(outdir)
        assert len(ds._fields_in_file) > 0
        field_names = [f[1] for f in ds._fields_in_file]
        assert "Density" in field_names

    def test_particle_types_detected(self, tmp_path):
        outdir = _create_mini_ramses_output(str(tmp_path))
        ds = MiniRAMSESDataset(outdir)
        assert "io" in ds.particle_types

    def test_str_representation(self, tmp_path):
        outdir = _create_mini_ramses_output(str(tmp_path))
        ds = MiniRAMSESDataset(outdir)
        assert str(ds) == "output_00001"
