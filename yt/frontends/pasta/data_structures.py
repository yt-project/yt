import glob
import os
import re
import struct
import weakref

import numpy as np

from yt.data_objects.index_subobjects.grid_patch import AMRGridPatch
from yt.data_objects.static_output import Dataset
from yt.fields.magnetic_field import get_magnetic_normalization
from yt.funcs import setdefaultattr
from yt.geometry.api import Geometry
from yt.geometry.grid_geometry_handler import GridIndex
from yt.utilities.chemical_formulas import compute_mu

from .definitions import geometry_map
from .fields import PastaFieldInfo


class PastaGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, id, index, level=-1):
        super().__init__(id, filename=index.index_filename, index=index)
        self.Parent = None
        self.Children = []
        self.Level = level


def parse_binary_header(filename):
    """Parse Pasta binary output header."""
    with open(filename, "rb") as f:
        endian_bytes = f.read(4)
        if len(endian_bytes) < 4:
            raise OSError(f"{filename} is too short to be a Pasta binary file")
        endian = struct.unpack("=i", endian_bytes)[0]
        bo = "<" if endian != 0 else ">"

        header_lines = []
        peek_pos = f.tell()
        line = f.readline()
        if b"# Pasta" not in line:
            raise TypeError(f"'# Pasta' marker not found in {filename}")
        header_lines.append(line)

        nfre = 0
        while True:
            peek_pos = f.tell()
            line = f.readline()
            if not line:
                break
            if line.startswith(b"# number of frequencies"):
                header_lines.append(line)
                if b"=" in line:
                    try:
                        nfre = int(line.split(b"=")[1].strip())
                    except ValueError:
                        nfre = 0
                if nfre > 1:
                    fre_line = f.readline()
                    header_lines.append(fre_line)
                continue
            if not line.startswith(b"#"):
                f.seek(peek_pos)
                break
            header_lines.append(line)
            if line.startswith(b"# per-block time table"):
                break

        def _val(prefix):
            for hl in header_lines:
                if hl.startswith(prefix) and b"=" in hl:
                    return hl.split(b"=", 1)[1].strip().decode("utf-8")
            return None

        # Parse header attributes
        time_str = None
        for hl in header_lines:
            if b"# Pasta" in hl and b"=" in hl:
                time_str = hl.split(b"=", 1)[1].strip().decode("utf-8")
                break
        current_time = float(time_str) if time_str is not None else 0.0

        coord_type = _val(b"# coord_type") or "Cartesian"
        cycle_val = _val(b"# cycle")
        cycle = int(cycle_val) if cycle_val is not None else 0

        grid_type_1 = int(_val(b"# grid_type_1") or 1)
        grid_type_2 = int(_val(b"# grid_type_2") or 1)
        grid_type_3 = int(_val(b"# grid_type_3") or 1)

        cgs_units = {
            "density": float(_val(b"# density unit/cgs") or 1.0),
            "length": float(_val(b"# length unit/cgs") or 1.0),
            "time": float(_val(b"# time unit/cgs") or 1.0),
            "velocity": float(_val(b"# velocity unit/cgs") or 1.0),
            "pressure": float(_val(b"# pressure unit/cgs") or 1.0),
            "temperature": float(_val(b"# temperature unit/cgs") or 1.0),
            "Er": float(_val(b"# radiation energy density unit/cgs") or 1.0),
            "Fr": float(_val(b"# radiation flux unit/cgs") or 1.0),
        }

        coord_data_size = int(_val(b"# size of coordinate data") or 8)
        data_size = int(_val(b"# size of variable data") or 8)
        num_vars = int(_val(b"# number of variables") or 0)
        num_ranks = int(_val(b"# number of MPI ranks") or 1)

        var_names = []
        for hl in header_lines:
            if hl.startswith(b"# variables:"):
                names = hl.split(b":", 1)[1].split()
                var_names = [n.decode("utf-8") for n in names]
                break

        num_face_field = int(_val(b"# face_field") or 0)
        num_edge_field = int(_val(b"# edge_field") or 0)

        redshift_str = _val(b"# redshift")
        scale_factor_str = _val(b"# scale_factor")
        redshift = float(redshift_str) if redshift_str is not None else None
        scale_factor = float(scale_factor_str) if scale_factor_str is not None else None

        # 6 integers for block info: {io_nbd, nblocks, nvar, ncells3, ncells2, ncells1}
        block_info_bytes = f.read(24)
        if len(block_info_bytes) < 24:
            raise OSError(f"Failed to read block info in {filename}")
        binfo = struct.unpack(bo + "6i", block_info_bytes)

        block_dic = {
            "nbd": int(binfo[0]),
            "nblocks": int(binfo[1]),
            "nvar": int(binfo[2]),
            "ncells3": int(binfo[3]),
            "ncells2": int(binfo[4]),
            "ncells1": int(binfo[5]),
        }

        data_start = f.tell()

    head = {
        "time": current_time,
        "coord_type": coord_type.strip(),
        "cycle": cycle,
        "grid_type": (grid_type_1, grid_type_2, grid_type_3),
        "cgs_units": cgs_units,
        "coord_data_size": coord_data_size,
        "data_size": data_size,
        "num_vars": num_vars,
        "num_ranks": num_ranks,
        "var_names": var_names,
        "num_face_field": num_face_field,
        "num_edge_field": num_edge_field,
        "redshift": redshift,
        "scale_factor": scale_factor,
        "byte_order": bo,
        "data_start": data_start,
    }
    return head, block_dic


def parse_particle_header(filename):
    """Parse Pasta particle dump header."""
    if not os.path.isfile(filename):
        return None
    with open(filename, "rb") as f:
        endian_bytes = f.read(4)
        if len(endian_bytes) < 4:
            return None
        endian = struct.unpack("=i", endian_bytes)[0]
        if endian not in (0, 1):
            return None
        bo = "<" if endian == 1 else ">"

        header_lines = []
        int_names, real_names = [], []
        head = {"byte_order": bo}

        while True:
            pos = f.tell()
            line = f.readline()
            if not line.startswith(b"#"):
                f.seek(pos)
                break
            header_lines.append(line)
            text = line[1:].decode("utf-8", errors="replace").strip()
            if text.startswith("int properties:"):
                int_names = text.split(":", 1)[1].split()
            elif text.startswith("real properties:"):
                real_names = text.split(":", 1)[1].split()
            elif "=" in text:
                k, v = [s.strip() for s in text.split("=", 1)]
                try:
                    head[k] = float(v) if ("." in v or "e" in v or "E" in v) else int(v)
                except ValueError:
                    head[k] = v

        ints_bytes = f.read(16)
        if len(ints_bytes) < 16:
            return None
        ints = np.frombuffer(ints_bytes, dtype=bo + "i4")
        nint, nreal, isz, rsz = [int(v) for v in ints]

        npar_bytes = f.read(8)
        nblk_bytes = f.read(8)
        if len(npar_bytes) < 8 or len(nblk_bytes) < 8:
            return None
        npar = int(np.frombuffer(npar_bytes, dtype=bo + "i8")[0])
        nblk = int(np.frombuffer(nblk_bytes, dtype=bo + "i8")[0])

        block_count = np.frombuffer(f.read(8 * nblk), dtype=bo + "i8").copy()
        data_start = f.tell()

    prop_offsets = {}
    prop_dtypes = {}
    curr_offset = data_start
    int_dtype = bo + "i8"
    real_dtype = bo + ("f8" if rsz == 8 else "f4")

    for name in int_names:
        prop_offsets[name] = curr_offset
        prop_dtypes[name] = int_dtype
        curr_offset += isz * npar

    for name in real_names:
        prop_offsets[name] = curr_offset
        prop_dtypes[name] = real_dtype
        curr_offset += rsz * npar

    def _val(prefix):
        for hl in header_lines:
            if hl.startswith(prefix) and b"=" in hl:
                return hl.split(b"=", 1)[1].strip().decode("utf-8")
        return None

    cgs_units = {
        "density": float(_val(b"# density unit/cgs") or 1.0),
        "length": float(_val(b"# length unit/cgs") or 1.0),
        "time": float(_val(b"# time unit/cgs") or 1.0),
        "velocity": float(_val(b"# velocity unit/cgs") or 1.0),
        "pressure": float(_val(b"# pressure unit/cgs") or 1.0),
        "temperature": float(_val(b"# temperature unit/cgs") or 1.0),
    }

    time_str = None
    for hl in header_lines:
        if b"# Pasta" in hl and b"=" in hl:
            time_str = hl.split(b"=", 1)[1].strip().decode("utf-8")
            break
    current_time = float(time_str) if time_str is not None else head.get("time", 0.0)
    cycle = int(_val(b"# cycle") or head.get("cycle", 0))
    coord_type = _val(b"# coord_type") or head.get("coord_type", "Cartesian")

    return {
        "head": head,
        "nint": nint,
        "nreal": nreal,
        "isz": isz,
        "rsz": rsz,
        "npar": npar,
        "nblk": nblk,
        "block_count": block_count,
        "int_names": int_names,
        "real_names": real_names,
        "prop_offsets": prop_offsets,
        "prop_dtypes": prop_dtypes,
        "data_start": data_start,
        "byte_order": bo,
        "cgs_units": cgs_units,
        "time": current_time,
        "cycle": cycle,
        "coord_type": coord_type.strip(),
        "redshift": head.get("redshift"),
        "scale_factor": head.get("scale_factor"),
    }


def parse_filename_info(filename):
    """Extract label/prefix and cycle from a Pasta output filename."""
    basename = os.path.basename(filename)
    stem, ext = os.path.splitext(basename)
    parts = stem.split(".")
    cycle = None
    cycle_str = None
    for part in reversed(parts):
        if part.isdigit():
            cycle_str = part
            cycle = int(part)
            break
    prefix = parts[0]
    if prefix.endswith("_par"):
        prefix = prefix[:-4]
    return {
        "prefix": prefix,
        "cycle": cycle,
        "cycle_str": cycle_str,
        "ext": ext,
        "basename": basename,
        "stem": stem,
    }


def find_companion_particle_file(grid_filename, cycle=None):
    """Find companion particle dump for a grid file."""
    dirname = os.path.dirname(os.path.abspath(grid_filename))
    info = parse_filename_info(grid_filename)
    label = info["prefix"]
    cycle_strs = []
    if cycle is not None:
        cycle_strs.append(f"{cycle:05d}")
        cycle_strs.append(str(cycle))
    if info["cycle_str"] is not None:
        cycle_strs.append(info["cycle_str"])

    cands = []
    for c_str in set(cycle_strs):
        cands.extend(
            [
                os.path.join(dirname, "..", "par", f"{label}.{c_str}.par"),
                os.path.join(dirname, "..", "par", f"{label}_par.{c_str}.par"),
                os.path.join(dirname, f"{label}.{c_str}.par"),
                os.path.join(dirname, f"{label}_par.{c_str}.par"),
                os.path.join(dirname, f"{label}.par.{c_str}.binary"),
            ]
        )
        cands.extend(
            sorted(
                glob.glob(os.path.join(dirname, "..", "par", f"{label}*{c_str}*.par"))
            )
        )
        cands.extend(sorted(glob.glob(os.path.join(dirname, f"{label}*{c_str}*.par"))))
        cands.extend(
            sorted(glob.glob(os.path.join(dirname, "..", "par", f"*{c_str}*.par")))
        )
        cands.extend(sorted(glob.glob(os.path.join(dirname, f"*{c_str}*.par"))))

    for c in cands:
        if os.path.isfile(c):
            return os.path.abspath(c)
    return None


def find_companion_grid_file(par_filename, cycle=None):
    """Find companion grid file for a particle dump."""
    dirname = os.path.dirname(os.path.abspath(par_filename))
    info = parse_filename_info(par_filename)
    label = info["prefix"]
    cycle_strs = []
    if cycle is not None:
        cycle_strs.append(f"{cycle:05d}")
        cycle_strs.append(str(cycle))
    if info["cycle_str"] is not None:
        cycle_strs.append(info["cycle_str"])

    cands = []
    for c_str in set(cycle_strs):
        for var in [
            "fluid_prim",
            "prim",
            "fluid_cons",
            "cons",
            "particle_rho",
            "gravity_phi",
            "field_bcc",
            "field_bf",
        ]:
            cands.extend(
                [
                    os.path.join(dirname, "..", "bin", f"{label}.{var}.{c_str}.binary"),
                    os.path.join(dirname, f"{label}.{var}.{c_str}.binary"),
                ]
            )
        cands.extend(
            [
                os.path.join(dirname, "..", "bin", f"{label}.{c_str}.binary"),
                os.path.join(dirname, f"{label}.{c_str}.binary"),
            ]
        )
        cands.extend(
            sorted(
                glob.glob(
                    os.path.join(dirname, "..", "bin", f"{label}*{c_str}*.binary")
                )
            )
        )
        cands.extend(
            sorted(glob.glob(os.path.join(dirname, f"{label}*{c_str}*.binary")))
        )
        cands.extend(
            sorted(glob.glob(os.path.join(dirname, "..", "bin", f"*{c_str}*.binary")))
        )
        cands.extend(sorted(glob.glob(os.path.join(dirname, f"*{c_str}*.binary"))))

    for c in cands:
        if os.path.isfile(c):
            return os.path.abspath(c)
    return None


def find_companion_grid_files(grid_filename, cycle=None):
    """Find all matching binary files for the same snapshot mesh and cycle."""
    grid_filename = os.path.abspath(grid_filename)
    dirname = os.path.dirname(grid_filename)
    info = parse_filename_info(grid_filename)
    prefix = info["prefix"]
    cycle_strs = []
    if cycle is not None:
        cycle_strs.append(f"{cycle:05d}")
        cycle_strs.append(str(cycle))
    if info["cycle_str"] is not None:
        cycle_strs.append(info["cycle_str"])

    candidates = [grid_filename]
    for c_str in set(cycle_strs):
        candidates.extend(
            sorted(glob.glob(os.path.join(dirname, f"{prefix}.*.{c_str}.binary")))
        )
        candidates.extend(
            sorted(glob.glob(os.path.join(dirname, f"{prefix}.{c_str}.binary")))
        )
        candidates.extend(sorted(glob.glob(os.path.join(dirname, f"*{c_str}*.binary"))))

    seen = {grid_filename}
    ordered_cands = [grid_filename]
    for c in candidates:
        c_abs = os.path.abspath(c)
        if c_abs not in seen and os.path.isfile(c_abs):
            seen.add(c_abs)
            ordered_cands.append(c_abs)

    return ordered_cands


def parse_cosmology_metadata(filename, cycle=None):
    """Parse cosmology metadata from companion files or logs."""
    dirname = os.path.dirname(os.path.abspath(filename))
    info = parse_filename_info(filename)
    prefix = info["prefix"]
    cycle_str = info["cycle_str"]
    if cycle is not None:
        cycle_str = f"{cycle:05d}"

    cosmo_params = {}

    # 1. Check companion *.fof_halos.h5
    h5_cands = []
    if cycle_str:
        h5_cands.extend(
            [
                os.path.join(dirname, f"{prefix}.{cycle_str}.fof_halos.h5"),
                os.path.join(dirname, f"{prefix}_{cycle_str}.fof_halos.h5"),
                os.path.join(
                    dirname, "..", "bin", f"{prefix}.{cycle_str}.fof_halos.h5"
                ),
            ]
        )
        h5_cands.extend(
            sorted(glob.glob(os.path.join(dirname, f"*{cycle_str}*.fof_halos.h5")))
        )
        h5_cands.extend(
            sorted(
                glob.glob(
                    os.path.join(dirname, "..", "bin", f"*{cycle_str}*.fof_halos.h5")
                )
            )
        )

    for cand in h5_cands:
        if os.path.isfile(cand):
            try:
                import h5py

                with h5py.File(cand, "r") as hf:
                    if "Header" in hf:
                        hattrs = hf["Header"].attrs
                        if "BoxSize" in hattrs:
                            cosmo_params["comoving_box_size"] = float(hattrs["BoxSize"])
                        if "HubbleParam" in hattrs:
                            cosmo_params["hubble_constant"] = float(
                                hattrs["HubbleParam"]
                            )
                        if "Omega0" in hattrs:
                            cosmo_params["omega_matter"] = float(hattrs["Omega0"])
                        elif "OmegaM" in hattrs:
                            cosmo_params["omega_matter"] = float(hattrs["OmegaM"])
                        if "OmegaLambda" in hattrs:
                            cosmo_params["omega_lambda"] = float(hattrs["OmegaLambda"])
                        elif "OmegaL" in hattrs:
                            cosmo_params["omega_lambda"] = float(hattrs["OmegaL"])
                        if "Redshift" in hattrs:
                            cosmo_params["redshift"] = float(hattrs["Redshift"])
                        if "ScaleFactor" in hattrs:
                            cosmo_params["scale_factor"] = float(hattrs["ScaleFactor"])
                        cosmo_params["cosmological_simulation"] = 1
                        break
            except Exception:
                pass

    # 2. Check pasta_logs.txt
    log_cands = [
        os.path.join(dirname, "pasta_logs.txt"),
        os.path.join(dirname, "..", "pasta_logs.txt"),
        os.path.join(dirname, "..", "..", "pasta_logs.txt"),
    ]
    for cand in log_cands:
        if os.path.isfile(cand):
            try:
                with open(cand) as f:
                    content = f.read()
                if "Cosmology:" in content:
                    cosmo_sec = content.split("Cosmology:")[1].split(
                        "----------------"
                    )[0]
                    lines = [ln.strip() for ln in cosmo_sec.splitlines()]
                    for ln in lines:
                        if ":" in ln:
                            k, v = [s.strip() for s in ln.split(":", 1)]
                            if k == "comoving_coordinates":
                                if v.upper() in ("ON", "TRUE", "1"):
                                    cosmo_params["cosmological_simulation"] = 1
                                elif v.upper() in ("OFF", "FALSE", "0"):
                                    cosmo_params["cosmological_simulation"] = 0
                            elif (
                                k == "hubble_constant_now"
                                and "hubble_constant" not in cosmo_params
                            ):
                                cosmo_params["hubble_constant"] = float(v)
                            elif (
                                k == "omega_matter_now"
                                and "omega_matter" not in cosmo_params
                            ):
                                cosmo_params["omega_matter"] = float(v)
                            elif (
                                k == "omega_lambda_now"
                                and "omega_lambda" not in cosmo_params
                            ):
                                cosmo_params["omega_lambda"] = float(v)
                            elif (
                                k == "omega_radiation_now"
                                and "omega_radiation" not in cosmo_params
                            ):
                                cosmo_params["omega_radiation"] = float(v)
                            elif (
                                k == "comoving_box_size"
                                and "comoving_box_size" not in cosmo_params
                            ):
                                cosmo_params["comoving_box_size"] = float(v)
                            elif (
                                k == "initial_redshift"
                                and "initial_redshift" not in cosmo_params
                            ):
                                cosmo_params["initial_redshift"] = float(v)
            except Exception:
                pass

    # 3. Check *.conf or input.txt
    conf_cands = [
        os.path.join(dirname, "..", f"{prefix}.conf"),
        os.path.join(dirname, "..", "..", f"{prefix}.conf"),
        os.path.join(dirname, "input.txt"),
        os.path.join(dirname, "..", "input.txt"),
    ]
    for cand in conf_cands:
        if os.path.isfile(cand):
            try:
                with open(cand) as f:
                    for line in f:
                        line = line.strip()
                        if line.startswith("#") or not line:
                            continue
                        if "=" in line:
                            k, v = [s.strip() for s in line.split("=", 1)]
                            kl = k.lower()
                            if (
                                kl in ("boxlength", "box_size", "comoving_box_size")
                                and "comoving_box_size" not in cosmo_params
                            ):
                                cosmo_params["comoving_box_size"] = float(v)
                            elif (
                                kl in ("omega_m", "omegam")
                                and "omega_matter" not in cosmo_params
                            ):
                                cosmo_params["omega_matter"] = float(v)
                            elif (
                                kl in ("omega_l", "omegal", "omega_lambda")
                                and "omega_lambda" not in cosmo_params
                            ):
                                cosmo_params["omega_lambda"] = float(v)
                            elif (
                                kl in ("h0", "hubble0")
                                and "hubble_constant" not in cosmo_params
                            ):
                                h_val = float(v)
                                if h_val > 10.0:
                                    h_val /= 100.0
                                cosmo_params["hubble_constant"] = h_val
            except Exception:
                pass

    return cosmo_params


def parse_pasta_logs_bc(filename):
    """Attempt to parse boundary condition flags from pasta_logs.txt."""
    dirname = os.path.dirname(os.path.abspath(filename))
    candidates = [
        os.path.join(dirname, "pasta_logs.txt"),
        os.path.join(dirname, "..", "pasta_logs.txt"),
        os.path.join(dirname, "..", "..", "pasta_logs.txt"),
    ]
    for cand in candidates:
        if os.path.isfile(cand):
            try:
                with open(cand) as f:
                    content = f.read()
                # Find active boundary flags
                match = re.search(
                    r"Active boundary flags.*?\n\s*fluid\s*:\s*(.+)",
                    content,
                )
                if match:
                    line = match.group(1)
                    flags = {}
                    for item in line.split():
                        if "=" in item:
                            k, v = item.split("=", 1)
                            flags[k.lower()] = v.lower()
                    p1 = (
                        flags.get("leftx1") == "periodic"
                        and flags.get("rightx1") == "periodic"
                    )
                    p2 = (
                        flags.get("leftx2") == "periodic"
                        and flags.get("rightx2") == "periodic"
                    )
                    p3 = (
                        flags.get("leftx3") == "periodic"
                        and flags.get("rightx3") == "periodic"
                    )
                    return (p1, p2, p3)
            except Exception:
                pass
    return None


class PastaHierarchy(GridIndex):
    grid = PastaGrid
    _dataset_type = "pasta"

    def __init__(self, ds, dataset_type="pasta"):
        self.dataset_type = dataset_type
        self.dataset = weakref.proxy(ds)
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        self.float_type = np.float64
        super().__init__(ds, dataset_type)

    def _detect_output_fields(self):
        self.field_list = [("pasta", k) for k in self.dataset._var_names]
        if self.dataset._head.get("num_face_field", 0) > 0:
            for i in range(1, 4):
                fld = ("pasta", f"field_bf_{i}")
                if fld not in self.field_list:
                    self.field_list.append(fld)
        if self.dataset._head.get("num_edge_field", 0) > 0:
            for i in range(1, 4):
                fld = ("pasta", f"field_emf_{i}")
                if fld not in self.field_list:
                    self.field_list.append(fld)
        if self.dataset._has_particles:
            pinfo = self.dataset._particle_info
            for name in pinfo["int_names"] + pinfo["real_names"]:
                self.field_list.append(("io", name))
            self.field_list.extend(
                [
                    ("io", "particle_position_x"),
                    ("io", "particle_position_y"),
                    ("io", "particle_position_z"),
                    ("io", "particle_velocity_x"),
                    ("io", "particle_velocity_y"),
                    ("io", "particle_velocity_z"),
                ]
            )
            if "mass" in pinfo["real_names"]:
                self.field_list.append(("io", "particle_mass"))
            if "pid" in pinfo["int_names"]:
                self.field_list.append(("io", "particle_index"))
            if self.dataset.geometry is Geometry.CYLINDRICAL:
                self.field_list.extend(
                    [
                        ("io", "particle_position_r"),
                        ("io", "particle_position_theta"),
                        ("io", "particle_velocity_r"),
                        ("io", "particle_velocity_theta"),
                    ]
                )
            elif self.dataset.geometry is Geometry.SPHERICAL:
                self.field_list.extend(
                    [
                        ("io", "particle_position_r"),
                        ("io", "particle_position_theta"),
                        ("io", "particle_position_phi"),
                        ("io", "particle_velocity_r"),
                        ("io", "particle_velocity_theta"),
                        ("io", "particle_velocity_phi"),
                    ]
                )

    def _count_grids(self):
        if self.dataset.parameter_filename.endswith(".par"):
            self.num_grids = 1
        else:
            self.num_grids = self.dataset._block_info["nblocks"]

    def _parse_index(self):
        if self.dataset.parameter_filename.endswith(".par"):
            self.num_grids = 1
            self.grid_left_edge = self.ds.arr(
                [self.dataset.domain_left_edge], "code_length"
            )
            self.grid_right_edge = self.ds.arr(
                [self.dataset.domain_right_edge], "code_length"
            )
            self.grid_dimensions = np.array([[1, 1, 1]], dtype="int32")
            self.grid_levels.flat[:] = [0]
            self.max_level = 0
            self.grids = np.empty(1, dtype="object")
            self.grids[0] = self.grid(0, self, 0)
            npar = self.dataset._particle_info["npar"]
            self.grid_particle_count = np.array([[npar]], dtype="int64")
            self.particle_indices = np.array([0, npar], dtype="int64")
            self.grids[0].NumberOfParticles = npar
            return

        head = self.dataset._head
        info = self.dataset._block_info
        bo = head["byte_order"]
        csz = head["coord_data_size"]
        dsz = head["data_size"]
        cfmt = bo + ("f8" if csz == 8 else "f4")

        nblocks = info["nblocks"]
        nc1, nc2, nc3 = info["ncells1"], info["ncells2"], info["ncells3"]
        nvar = info["nvar"]
        tot_cells = nc1 * nc2 * nc3

        # Compute block record byte stride
        coord_count = (nc3 + 1) + nc3 + (nc2 + 1) + nc2 + (nc1 + 1) + nc1
        coord_bytes = coord_count * csz
        data_bytes = nvar * tot_cells * dsz

        face_bytes = 0
        if head["num_face_field"] > 0:
            face_count = (
                nc3 * nc2 * (nc1 + 1) + nc3 * (nc2 + 1) * nc1 + (nc3 + 1) * nc2 * nc1
            )
            face_bytes = face_count * dsz

        edge_bytes = 0
        if head["num_edge_field"] > 0:
            edge_count = (
                (nc3 + 1) * (nc2 + 1) * nc1
                + (nc3 + 1) * nc2 * (nc1 + 1)
                + nc3 * (nc2 + 1) * (nc1 + 1)
            )
            edge_bytes = edge_count * dsz

        block_bytes = coord_bytes + data_bytes + face_bytes + edge_bytes
        self.dataset._block_bytes = block_bytes
        self.dataset._coord_bytes = coord_bytes
        self.dataset._var_bytes = tot_cells * dsz

        self.grid_left_edge = np.zeros((nblocks, 3), dtype="float64")
        self.grid_right_edge = np.zeros((nblocks, 3), dtype="float64")
        self.grid_dimensions = np.zeros((nblocks, 3), dtype="int32")

        # Read face coordinates for all blocks
        with open(self.index_filename, "rb") as f:
            for n in range(nblocks):
                f.seek(head["data_start"] + n * block_bytes)
                x3f = np.fromfile(f, dtype=cfmt, count=nc3 + 1)
                f.seek(csz * nc3, os.SEEK_CUR)  # skip x3c
                x2f = np.fromfile(f, dtype=cfmt, count=nc2 + 1)
                f.seek(csz * nc2, os.SEEK_CUR)  # skip x2c
                x1f = np.fromfile(f, dtype=cfmt, count=nc1 + 1)

                self.grid_left_edge[n] = [x1f[0], x2f[0], x3f[0]]
                self.grid_right_edge[n] = [x1f[-1], x2f[-1], x3f[-1]]
                self.grid_dimensions[n] = [nc1, nc2, nc3]

        self.grid_left_edge = self.ds.arr(self.grid_left_edge, "code_length")
        self.grid_right_edge = self.ds.arr(self.grid_right_edge, "code_length")

        # Calculate cell widths and levels
        dx_blocks = (self.grid_right_edge[:, 0].d - self.grid_left_edge[:, 0].d) / nc1
        dx_root = np.max(dx_blocks)
        levels = np.rint(np.log2(np.maximum(dx_root / dx_blocks, 1.0))).astype("int32")

        # Check for companion .grid file
        companion_grid = self._find_companion_grid()
        if companion_grid is not None:
            grid_levels = self._read_grid_levels(companion_grid, nblocks)
            if grid_levels is not None:
                levels = grid_levels

        self.grid_levels.flat[:] = levels
        self.max_level = int(np.max(levels)) if len(levels) > 0 else 0

        self.grids = np.empty(self.num_grids, dtype="object")
        for i in range(self.num_grids):
            self.grids[i] = self.grid(i, self, levels[i])

        if self.dataset._has_particles:
            pinfo = self.dataset._particle_info
            block_count = pinfo["block_count"]
            if len(block_count) == self.num_grids:
                self.grid_particle_count = block_count.reshape(self.num_grids, 1)
            else:
                self.grid_particle_count = np.zeros((self.num_grids, 1), dtype="int64")
        else:
            self.grid_particle_count = np.zeros((self.num_grids, 1), dtype="int64")

        self.particle_indices = np.zeros(self.num_grids + 1, dtype="int64")
        np.cumsum(self.grid_particle_count.flatten(), out=self.particle_indices[1:])
        for i in range(self.num_grids):
            self.grids[i].NumberOfParticles = int(self.grid_particle_count[i, 0])

        # Setup 1D/2D boundaries if needed
        if self.dataset.dimensionality <= 2:
            self.grid_left_edge[:, 2] = self.dataset.domain_left_edge[2]
            self.grid_right_edge[:, 2] = self.dataset.domain_right_edge[2]
        if self.dataset.dimensionality == 1:
            self.grid_left_edge[:, 1:] = self.dataset.domain_left_edge[1:]
            self.grid_right_edge[:, 1:] = self.dataset.domain_right_edge[1:]

    def _find_companion_grid(self):
        dirname = os.path.dirname(os.path.abspath(self.index_filename))
        basename = os.path.basename(self.index_filename)
        # Check adjacent .grid or ../grid/<name>.grid
        stem = basename.rsplit(".", 1)[0]
        cands = [
            os.path.join(dirname, f"{stem}.grid"),
            os.path.join(dirname, "..", "grid", f"{stem}.grid"),
        ]
        for c in cands:
            if os.path.isfile(c):
                return c
        return None

    def _read_grid_levels(self, grid_file, expected_count):
        try:
            levels = []
            with open(grid_file) as f:
                for line in f:
                    if line.startswith("#") or "=" in line:
                        continue
                    parts = line.strip().split()
                    if len(parts) >= 4:
                        levels.append(int(parts[3]))
            if len(levels) == expected_count:
                return np.array(levels, dtype="int32")
        except Exception:
            pass
        return None

    def _populate_grid_objects(self):
        for g in self.grids:
            g._prepare_grid()
            g._setup_dx()


class PastaDataset(Dataset):
    _field_info_class = PastaFieldInfo
    _dataset_type = "pasta"
    _index_class = PastaHierarchy

    def __init__(
        self,
        filename,
        dataset_type="pasta",
        storage_filename=None,
        parameters=None,
        cosmology_parameters=None,
        units_override=None,
        unit_system="cgs",
        default_species_fields=None,
        magnetic_normalization="gaussian",
        particle_filename=None,
        *,
        axis_order=None,
    ):
        self.fluid_types += ("pasta",)
        if parameters is None:
            parameters = {}
        self.specified_parameters = parameters
        self._cosmology_parameters = cosmology_parameters
        if units_override is None:
            units_override = {}
        self._magnetic_factor = get_magnetic_normalization(magnetic_normalization)

        self.particle_filename = particle_filename
        self._particle_info = None
        self._has_particles = False
        self.primordial_chemistry = None

        if filename.endswith(".par"):
            self.particle_filename = os.path.abspath(filename)
            companion_grid = find_companion_grid_file(filename)
            if companion_grid is not None:
                filename = companion_grid
            else:
                self._particle_info = parse_particle_header(self.particle_filename)
                if self._particle_info is not None and self._particle_info["npar"] > 0:
                    self._has_particles = True

        if not filename.endswith(".par"):
            head, block_info = parse_binary_header(filename)
            self._head = head
            self._block_info = block_info

            companion_files = find_companion_grid_files(filename, cycle=head["cycle"])
            field_sources = {}
            valid_grid_files = []

            nc1, nc2, nc3 = (
                block_info["ncells1"],
                block_info["ncells2"],
                block_info["ncells3"],
            )
            tot_cells = nc1 * nc2 * nc3

            for fn in companion_files:
                try:
                    h, b = (
                        (head, block_info)
                        if fn == filename
                        else parse_binary_header(fn)
                    )
                    if (
                        b["nblocks"] != block_info["nblocks"]
                        or b["ncells1"] != nc1
                        or b["ncells2"] != nc2
                        or b["ncells3"] != nc3
                        or h["cycle"] != head["cycle"]
                    ):
                        continue

                    valid_grid_files.append(fn)
                    csz = h["coord_data_size"]
                    dsz = h["data_size"]
                    coord_count = (nc3 + 1) + nc3 + (nc2 + 1) + nc2 + (nc1 + 1) + nc1
                    coord_bytes = coord_count * csz
                    var_bytes = tot_cells * dsz
                    face_bytes = 0
                    if h["num_face_field"] > 0:
                        face_count = (
                            nc3 * nc2 * (nc1 + 1)
                            + nc3 * (nc2 + 1) * nc1
                            + (nc3 + 1) * nc2 * nc1
                        )
                        face_bytes = face_count * dsz
                    edge_bytes = 0
                    if h["num_edge_field"] > 0:
                        edge_count = (
                            (nc3 + 1) * (nc2 + 1) * nc1
                            + (nc3 + 1) * nc2 * (nc1 + 1)
                            + nc3 * (nc2 + 1) * (nc1 + 1)
                        )
                        edge_bytes = edge_count * dsz
                    block_bytes = (
                        coord_bytes + b["nvar"] * var_bytes + face_bytes + edge_bytes
                    )
                    dfmt = h["byte_order"] + ("f8" if dsz == 8 else "f4")

                    for v_idx, v_name in enumerate(h["var_names"]):
                        if v_name not in field_sources:
                            field_sources[v_name] = {
                                "filename": fn,
                                "head": h,
                                "block_info": b,
                                "v_idx": v_idx,
                                "is_face": False,
                                "is_edge": False,
                                "data_start": h["data_start"],
                                "coord_bytes": coord_bytes,
                                "var_bytes": var_bytes,
                                "block_bytes": block_bytes,
                                "dfmt": dfmt,
                                "dsz": dsz,
                            }

                    if h["num_face_field"] > 0:
                        for comp in (1, 2, 3):
                            fname = f"field_bf_{comp}"
                            if fname not in field_sources:
                                field_sources[fname] = {
                                    "filename": fn,
                                    "head": h,
                                    "block_info": b,
                                    "v_idx": comp - 1,
                                    "is_face": True,
                                    "is_edge": False,
                                    "data_start": h["data_start"],
                                    "coord_bytes": coord_bytes,
                                    "var_bytes": var_bytes,
                                    "block_bytes": block_bytes,
                                    "dfmt": dfmt,
                                    "dsz": dsz,
                                }

                    if h["num_edge_field"] > 0:
                        for comp in (1, 2, 3):
                            fname = f"field_emf_{comp}"
                            if fname not in field_sources:
                                field_sources[fname] = {
                                    "filename": fn,
                                    "head": h,
                                    "block_info": b,
                                    "v_idx": comp - 1,
                                    "is_face": False,
                                    "is_edge": True,
                                    "data_start": h["data_start"],
                                    "coord_bytes": coord_bytes,
                                    "var_bytes": var_bytes,
                                    "block_bytes": block_bytes,
                                    "dfmt": dfmt,
                                    "dsz": dsz,
                                }
                except Exception:
                    continue

            self._field_sources = field_sources
            self._companion_grid_files = valid_grid_files
            self._var_names = list(field_sources.keys())

            if self.particle_filename is None:
                self.particle_filename = find_companion_particle_file(
                    filename, cycle=head["cycle"]
                )
            if self.particle_filename is not None and os.path.isfile(
                self.particle_filename
            ):
                self._particle_info = parse_particle_header(self.particle_filename)
                if self._particle_info is not None and self._particle_info["npar"] > 0:
                    self._has_particles = True
            coord = head["coord_type"].strip().lower()

            cosmo_params = parse_cosmology_metadata(filename, cycle=head["cycle"])
            user_cosmo = self._cosmology_parameters or {}
            user_params = self.specified_parameters or {}
            for k in (
                "hubble_constant",
                "omega_matter",
                "omega_lambda",
                "omega_radiation",
                "comoving_box_size",
            ):
                if k in user_params:
                    cosmo_params[k] = float(user_params[k])
                elif k in user_cosmo:
                    cosmo_params[k] = float(user_cosmo[k])
            if "h0" in user_params:
                h0_val = float(user_params["h0"])
                cosmo_params["hubble_constant"] = (
                    h0_val / 100.0 if h0_val > 10.0 else h0_val
                )
            elif "h0" in user_cosmo:
                h0_val = float(user_cosmo["h0"])
                cosmo_params["hubble_constant"] = (
                    h0_val / 100.0 if h0_val > 10.0 else h0_val
                )

            if (
                user_params.get("cosmological_simulation") == 1
                or user_cosmo.get("cosmological_simulation") == 1
            ):
                cosmo_params["cosmological_simulation"] = 1
            elif (
                user_params.get("cosmological_simulation") == 0
                or user_cosmo.get("cosmological_simulation") == 0
            ):
                cosmo_params["cosmological_simulation"] = 0
            elif head.get("redshift") is not None:
                cosmo_params["cosmological_simulation"] = 1
                cosmo_params["redshift"] = head["redshift"]
                cosmo_params["scale_factor"] = head.get("scale_factor", 1.0)
            elif (
                self._particle_info and self._particle_info.get("redshift") is not None
            ):
                cosmo_params["cosmological_simulation"] = 1
                cosmo_params["redshift"] = self._particle_info["redshift"]
                cosmo_params["scale_factor"] = self._particle_info.get(
                    "scale_factor", 1.0
                )
            self._cosmo_params = cosmo_params
        else:
            self._head = {
                "time": self._particle_info["time"],
                "coord_type": self._particle_info["coord_type"],
                "cycle": self._particle_info["cycle"],
                "cgs_units": self._particle_info["cgs_units"],
                "byte_order": self._particle_info["byte_order"],
            }
            self._block_info = {
                "nblocks": 1,
                "ncells1": 1,
                "ncells2": 1,
                "ncells3": 1,
                "nvar": 0,
            }
            self._var_names = []
            self._field_sources = {}
            coord = self._particle_info["coord_type"].strip().lower()

            cosmo_params = parse_cosmology_metadata(
                self.particle_filename, cycle=self._particle_info.get("cycle")
            )
            user_cosmo = self._cosmology_parameters or {}
            user_params = self.specified_parameters or {}
            for k in (
                "hubble_constant",
                "omega_matter",
                "omega_lambda",
                "omega_radiation",
                "comoving_box_size",
            ):
                if k in user_params:
                    cosmo_params[k] = float(user_params[k])
                elif k in user_cosmo:
                    cosmo_params[k] = float(user_cosmo[k])
            if "h0" in user_params:
                h0_val = float(user_params["h0"])
                cosmo_params["hubble_constant"] = (
                    h0_val / 100.0 if h0_val > 10.0 else h0_val
                )
            elif "h0" in user_cosmo:
                h0_val = float(user_cosmo["h0"])
                cosmo_params["hubble_constant"] = (
                    h0_val / 100.0 if h0_val > 10.0 else h0_val
                )

            if (
                user_params.get("cosmological_simulation") == 1
                or user_cosmo.get("cosmological_simulation") == 1
            ):
                cosmo_params["cosmological_simulation"] = 1
            elif (
                user_params.get("cosmological_simulation") == 0
                or user_cosmo.get("cosmological_simulation") == 0
            ):
                cosmo_params["cosmological_simulation"] = 0
            elif self._particle_info.get("redshift") is not None:
                cosmo_params["cosmological_simulation"] = 1
                cosmo_params["redshift"] = self._particle_info["redshift"]
                cosmo_params["scale_factor"] = self._particle_info.get(
                    "scale_factor", 1.0
                )
            self._cosmo_params = cosmo_params

        self.geometry = Geometry(geometry_map.get(coord, "cartesian"))
        if self.geometry is Geometry.CYLINDRICAL:
            axis_order = ("r", "theta", "z")
        elif self.geometry is Geometry.SPHERICAL:
            axis_order = ("r", "theta", "phi")
        else:
            axis_order = ("x", "y", "z")

        super().__init__(
            filename,
            dataset_type,
            units_override=units_override,
            unit_system=unit_system,
            default_species_fields=default_species_fields,
            axis_order=axis_order,
        )
        self.storage_filename = storage_filename

    def _parse_parameter_file(self):
        if self.parameter_filename.endswith(".par"):
            pinfo = self._particle_info
            bo = pinfo["byte_order"]
            rdt = bo + ("f8" if pinfo["rsz"] == 8 else "f4")
            with open(self.parameter_filename, "rb") as f:
                f.seek(pinfo["prop_offsets"]["x1"])
                x1 = np.fromfile(f, dtype=rdt, count=pinfo["npar"])
                f.seek(pinfo["prop_offsets"]["x2"])
                x2 = np.fromfile(f, dtype=rdt, count=pinfo["npar"])
                f.seek(pinfo["prop_offsets"]["x3"])
                x3 = np.fromfile(f, dtype=rdt, count=pinfo["npar"])
            x1_min, x1_max = float(np.min(x1)), float(np.max(x1))
            x2_min, x2_max = float(np.min(x2)), float(np.max(x2))
            x3_min, x3_max = float(np.min(x3)), float(np.max(x3))
            pad1 = 0.01 * (x1_max - x1_min) if x1_max > x1_min else 0.1
            pad2 = 0.01 * (x2_max - x2_min) if x2_max > x2_min else 0.1
            pad3 = 0.01 * (x3_max - x3_min) if x3_max > x3_min else 0.1
            self.domain_left_edge = np.array(
                [x1_min - pad1, x2_min - pad2, x3_min - pad3], dtype="float64"
            )
            self.domain_right_edge = np.array(
                [x1_max + pad1, x2_max + pad2, x3_max + pad3], dtype="float64"
            )
            self.domain_width = self.domain_right_edge - self.domain_left_edge
            self.domain_dimensions = np.array([1, 1, 1], dtype="int64")
            self.dimensionality = 3
            self.current_time = pinfo["time"]
            self.current_cycle = pinfo["cycle"]
            cosmo_params = getattr(self, "_cosmo_params", {})
            if cosmo_params.get("cosmological_simulation", 0) == 1:
                self.cosmological_simulation = 1
                self.current_redshift = cosmo_params.get("redshift", 0.0)
                self.omega_matter = cosmo_params.get("omega_matter", 0.3)
                self.omega_lambda = cosmo_params.get("omega_lambda", 0.7)
                self.omega_radiation = cosmo_params.get("omega_radiation", 0.0)
                self.hubble_constant = cosmo_params.get("hubble_constant", 0.7)
                scale_factor = cosmo_params.get("scale_factor", 1.0)
                if "comoving_box_size" in cosmo_params:
                    box_size = cosmo_params["comoving_box_size"]
                else:
                    z_i = scale_factor * (1.0 + self.current_redshift) - 1.0
                    mpc_cgs = 3.085677581282e24
                    h0 = self.hubble_constant
                    box_size = (
                        self._head["cgs_units"]["length"] * h0 * (1.0 + z_i) / mpc_cgs
                    )
                self.parameters["comoving_box_size"] = float(box_size)
                self.parameters["scale_factor"] = float(scale_factor)
            else:
                self.cosmological_simulation = 0
                self.current_redshift = 0.0
                self.omega_matter = 0.0
                self.omega_lambda = 0.0
                self.omega_radiation = 0.0
                self.hubble_constant = 0.0
            self.refine_by = 2
            self._periodicity = (False, False, False)
            self.gamma = float(self.specified_parameters.get("gamma", 5.0 / 3.0))
            self.mu = float(
                self.specified_parameters.get(
                    "mu", compute_mu(self.default_species_fields)
                )
            )
            self._block_bytes = 0
            self._coord_bytes = 0
            self._var_bytes = 0
            return

        head = self._head
        info = self._block_info
        bo = head["byte_order"]
        csz = head["coord_data_size"]
        cfmt = bo + ("f8" if csz == 8 else "f4")
        nblocks = info["nblocks"]
        nc1, nc2, nc3 = info["ncells1"], info["ncells2"], info["ncells3"]

        tot_cells = nc1 * nc2 * nc3
        coord_bytes = ((nc3 + 1) + nc3 + (nc2 + 1) + nc2 + (nc1 + 1) + nc1) * csz
        data_bytes = info["nvar"] * tot_cells * head["data_size"]
        face_bytes = 0
        if head["num_face_field"] > 0:
            face_count = (
                nc3 * nc2 * (nc1 + 1) + nc3 * (nc2 + 1) * nc1 + (nc3 + 1) * nc2 * nc1
            )
            face_bytes = face_count * head["data_size"]
        edge_bytes = 0
        if head["num_edge_field"] > 0:
            edge_count = (
                (nc3 + 1) * (nc2 + 1) * nc1
                + (nc3 + 1) * nc2 * (nc1 + 1)
                + nc3 * (nc2 + 1) * (nc1 + 1)
            )
            edge_bytes = edge_count * head["data_size"]
        block_bytes = coord_bytes + data_bytes + face_bytes + edge_bytes

        # Read overall domain extents and root dx from block faces
        x1_min, x1_max = np.inf, -np.inf
        x2_min, x2_max = np.inf, -np.inf
        x3_min, x3_max = np.inf, -np.inf
        dx1_max, dx2_max, dx3_max = 0.0, 0.0, 0.0

        with open(self.parameter_filename, "rb") as f:
            for n in range(nblocks):
                f.seek(head["data_start"] + n * block_bytes)
                x3f = np.fromfile(f, dtype=cfmt, count=nc3 + 1)
                f.seek(csz * nc3, os.SEEK_CUR)
                x2f = np.fromfile(f, dtype=cfmt, count=nc2 + 1)
                f.seek(csz * nc2, os.SEEK_CUR)
                x1f = np.fromfile(f, dtype=cfmt, count=nc1 + 1)

                x1_min = min(x1_min, float(x1f[0]))
                x1_max = max(x1_max, float(x1f[-1]))
                x2_min = min(x2_min, float(x2f[0]))
                x2_max = max(x2_max, float(x2f[-1]))
                x3_min = min(x3_min, float(x3f[0]))
                x3_max = max(x3_max, float(x3f[-1]))

                dx1_max = max(dx1_max, float(x1f[-1] - x1f[0]) / nc1)
                dx2_max = max(dx2_max, float(x2f[-1] - x2f[0]) / nc2)
                dx3_max = max(dx3_max, float(x3f[-1] - x3f[0]) / nc3)

        # Dimensionality
        if nc3 == 1 and nc2 == 1:
            self.dimensionality = 1
            if x2_min >= x2_max:
                x2_min, x2_max = 0.0, 1.0
            if x3_min >= x3_max:
                x3_min, x3_max = 0.0, 1.0
        elif nc3 == 1:
            self.dimensionality = 2
            if x3_min >= x3_max:
                x3_min, x3_max = 0.0, 1.0
        else:
            self.dimensionality = 3

        self.domain_left_edge = np.array([x1_min, x2_min, x3_min], dtype="float64")
        self.domain_right_edge = np.array([x1_max, x2_max, x3_max], dtype="float64")
        self.domain_width = self.domain_right_edge - self.domain_left_edge

        # Compute root domain dimensions
        nx = int(np.rint(self.domain_width[0] / dx1_max)) if dx1_max > 0 else 1
        ny = (
            int(np.rint(self.domain_width[1] / dx2_max))
            if (self.dimensionality >= 2 and dx2_max > 0)
            else 1
        )
        nz = (
            int(np.rint(self.domain_width[2] / dx3_max))
            if (self.dimensionality >= 3 and dx3_max > 0)
            else 1
        )
        self.domain_dimensions = np.array([nx, ny, nz], dtype="int64")

        self.current_time = head["time"]
        self.current_cycle = head["cycle"]
        cosmo_params = getattr(self, "_cosmo_params", {})
        if cosmo_params.get("cosmological_simulation", 0) == 1:
            self.cosmological_simulation = 1
            self.current_redshift = cosmo_params.get("redshift", 0.0)
            self.omega_matter = cosmo_params.get("omega_matter", 0.3)
            self.omega_lambda = cosmo_params.get("omega_lambda", 0.7)
            self.omega_radiation = cosmo_params.get("omega_radiation", 0.0)
            self.hubble_constant = cosmo_params.get("hubble_constant", 0.7)
            scale_factor = cosmo_params.get("scale_factor", 1.0)
            if "comoving_box_size" in cosmo_params:
                box_size = cosmo_params["comoving_box_size"]
            else:
                z_i = scale_factor * (1.0 + self.current_redshift) - 1.0
                mpc_cgs = 3.085677581282e24
                h0 = self.hubble_constant
                box_size = (
                    self._head["cgs_units"]["length"] * h0 * (1.0 + z_i) / mpc_cgs
                )
            self.parameters["comoving_box_size"] = float(box_size)
            self.parameters["scale_factor"] = float(scale_factor)
        else:
            self.cosmological_simulation = 0
            self.current_redshift = 0.0
            self.omega_matter = 0.0
            self.omega_lambda = 0.0
            self.omega_radiation = 0.0
            self.hubble_constant = 0.0
        self.refine_by = 2

        # Boundary conditions
        parsed_bc = parse_pasta_logs_bc(self.parameter_filename)
        if "periodicity" in self.specified_parameters:
            self._periodicity = tuple(self.specified_parameters["periodicity"])
        elif parsed_bc is not None:
            self._periodicity = parsed_bc
        elif self.geometry is Geometry.SPHERICAL:
            self._periodicity = (False, False, True)
        elif self.geometry is Geometry.CYLINDRICAL:
            self._periodicity = (False, True, True)
        else:
            self._periodicity = (True, True, True)

        self.gamma = float(self.specified_parameters.get("gamma", 5.0 / 3.0))
        self.mu = float(
            self.specified_parameters.get("mu", compute_mu(self.default_species_fields))
        )

        prim_chem = self.specified_parameters.get("chemistry/primordial_chemistry")
        if prim_chem is None:
            prim_chem = self.specified_parameters.get("primordial_chemistry")
        if prim_chem is None:
            prim_chem = self.parameters.get("chemistry/primordial_chemistry")
        if prim_chem is None:
            prim_chem = self.parameters.get("primordial_chemistry")

        if prim_chem is None and hasattr(self, "_var_names"):
            scalar_indices = set()
            for fn in self._var_names:
                if fn.startswith("scalar_cons_"):
                    try:
                        scalar_indices.add(int(fn.split("scalar_cons_")[1]))
                    except ValueError:
                        pass
                elif fn.startswith("scalar_prim_"):
                    try:
                        scalar_indices.add(int(fn.split("scalar_prim_")[1]))
                    except ValueError:
                        pass
            if scalar_indices:
                max_idx = max(scalar_indices)
                if max_idx >= 11:
                    prim_chem = 3
                elif max_idx >= 8:
                    prim_chem = 2
                elif max_idx >= 5 or 0 in scalar_indices:
                    prim_chem = 1
                else:
                    prim_chem = 0

        self.primordial_chemistry = int(prim_chem) if prim_chem is not None else 0
        self.parameters["primordial_chemistry"] = self.primordial_chemistry

    def _set_code_unit_attributes(self):
        units = self._head["cgs_units"]
        if self.cosmological_simulation:
            box_size = self.parameters.get("comoving_box_size", 1.0)
            setdefaultattr(self, "length_unit", self.quan(box_size, "Mpccm/h"))
            setdefaultattr(
                self,
                "mass_unit",
                self.quan(units["density"], "g/cm**3")
                * (self.length_unit.in_cgs()) ** 3,
            )
            setdefaultattr(self, "time_unit", self.quan(units["time"], "s"))
            setdefaultattr(self, "velocity_unit", self.quan(units["velocity"], "cm/s"))
            setdefaultattr(
                self, "temperature_unit", self.quan(units["temperature"], "K")
            )
        else:
            setdefaultattr(self, "length_unit", self.quan(units["length"], "cm"))
            setdefaultattr(self, "time_unit", self.quan(units["time"], "s"))
            mass = units["density"] * (units["length"] ** 3)
            setdefaultattr(self, "mass_unit", self.quan(mass, "g"))
            setdefaultattr(self, "velocity_unit", self.quan(units["velocity"], "cm/s"))
            setdefaultattr(
                self, "temperature_unit", self.quan(units["temperature"], "K")
            )

        # Magnetic unit in Gaussian cgs
        mag_unit = np.sqrt(
            self._magnetic_factor
            * self.mass_unit
            / (self.time_unit**2 * self.length_unit)
        )
        setdefaultattr(
            self, "magnetic_unit", self.quan(float(mag_unit.in_cgs()), "gauss")
        )

    @classmethod
    def _is_valid(cls, filename: str, *args, **kwargs) -> bool:
        if not os.path.isfile(filename):
            return False
        try:
            with open(filename, "rb") as f:
                endian_bytes = f.read(4)
                if len(endian_bytes) < 4:
                    return False
                endian = struct.unpack("=i", endian_bytes)[0]
                if endian not in (0, 1):
                    return False
                header = f.readline()
                return b"# Pasta" in header
        except Exception:
            return False

    @property
    def species_names(self):
        _ = self.derived_field_list
        return getattr(self.field_info, "species_names", [])

    def __str__(self):
        return os.path.basename(self.parameter_filename)
