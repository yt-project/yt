import shutil

import numpy as np

import yt
from yt.frontends.pasta.api import PastaDataset
from yt.testing import assert_allclose_units


def _write_pasta_binary(
    filename,
    coord_type="Cartesian",
    nblocks=1,
    ncells=(8, 8, 8),
    vars_to_write=(
        "fluid_den_prim",
        "fluid_vel_1",
        "fluid_vel_2",
        "fluid_vel_3",
        "fluid_pressure",
    ),
    num_face_field=0,
    num_edge_field=0,
    cycle=0,
    time=0.0,
    redshift=None,
    scale_factor=None,
    cgs_units=None,
    x1_range=(0.0, 1.0),
    x2_range=(0.0, 1.0),
    x3_range=(0.0, 1.0),
):
    nc1, nc2, nc3 = ncells
    nvar = len(vars_to_write)
    var_str = "  ".join(vars_to_write)

    if cgs_units is None:
        cgs_units = {}

    header = (
        f"# Pasta output at global time = {time:.16e}\n"
        f"# coord_type = {coord_type}\n"
        f"# cycle = {cycle}\n"
        "# grid_type_1 = 1\n"
        "# grid_type_2 = 1\n"
        "# grid_type_3 = 1\n"
        f"# density unit/cgs = {cgs_units.get('density', 1.0):.16e}\n"
        f"# length unit/cgs = {cgs_units.get('length', 1.0):.16e}\n"
        f"# time unit/cgs = {cgs_units.get('time', 1.0):.16e}\n"
        f"# velocity unit/cgs = {cgs_units.get('velocity', 1.0):.16e}\n"
        f"# pressure unit/cgs = {cgs_units.get('pressure', 1.0):.16e}\n"
        f"# temperature unit/cgs = {cgs_units.get('temperature', 1.0):.16e}\n"
        "# radiation energy density unit/cgs = 1.0000000000000000e+00\n"
        "# radiation flux unit/cgs = 1.0000000000000000e+00\n"
        "# size of coordinate data = 8\n"
        "# size of variable data = 8\n"
        f"# number of variables = {nvar}\n"
        "# number of MPI ranks = 1\n"
        f"# variables:  {var_str}\n"
        "# number of frequencies = 0\n"
    )
    if num_face_field > 0:
        header += f"# face_field = {num_face_field}\n"
    if num_edge_field > 0:
        header += f"# edge_field = {num_edge_field}\n"
    if redshift is not None:
        header += f"# redshift = {redshift:.16e}\n"
    if scale_factor is not None:
        header += f"# scale_factor = {scale_factor:.16e}\n"

    with open(filename, "wb") as f:
        f.write(np.int32(1).tobytes())
        f.write(header.encode("ascii"))
        binfo = np.array([0, nblocks, nvar, nc3, nc2, nc1], dtype="<i4")
        f.write(binfo.tobytes())

        tot_cells = nc1 * nc2 * nc3
        for b in range(nblocks):
            # Coordinates
            x3f = np.linspace(x3_range[0], x3_range[1], nc3 + 1, dtype="<f8")
            x3c = 0.5 * (x3f[:-1] + x3f[1:])
            x2f = np.linspace(x2_range[0], x2_range[1], nc2 + 1, dtype="<f8")
            x2c = 0.5 * (x2f[:-1] + x2f[1:])
            dx1 = (x1_range[1] - x1_range[0]) / nblocks
            x1f = np.linspace(
                x1_range[0] + b * dx1, x1_range[0] + (b + 1) * dx1, nc1 + 1, dtype="<f8"
            )
            x1c = 0.5 * (x1f[:-1] + x1f[1:])

            f.write(x3f.tobytes())
            f.write(x3c.tobytes())
            f.write(x2f.tobytes())
            f.write(x2c.tobytes())
            f.write(x1f.tobytes())
            f.write(x1c.tobytes())

            # Fluid variables
            for v_idx, _v_name in enumerate(vars_to_write):
                val = (b + 1.0) * (v_idx + 1.0)
                arr = np.full(tot_cells, val, dtype="<f8")
                f.write(arr.tobytes())

            # Staggered face fields if any
            if num_face_field > 0:
                f1 = np.ones((nc3, nc2, nc1 + 1), dtype="<f8") * 0.1
                f2 = np.ones((nc3, nc2 + 1, nc1), dtype="<f8") * 0.2
                f3 = np.ones((nc3 + 1, nc2, nc1), dtype="<f8") * 0.3
                f.write(f1.tobytes())
                f.write(f2.tobytes())
                f.write(f3.tobytes())

            # Staggered edge fields if any
            if num_edge_field > 0:
                e1 = np.ones((nc3 + 1, nc2 + 1, nc1), dtype="<f8") * 0.01
                e2 = np.ones((nc3 + 1, nc2, nc1 + 1), dtype="<f8") * 0.02
                e3 = np.ones((nc3, nc2 + 1, nc1 + 1), dtype="<f8") * 0.03
                f.write(e1.tobytes())
                f.write(e2.tobytes())
                f.write(e3.tobytes())


def _write_pasta_particle(
    filename,
    coord_type="Cartesian",
    nblocks=1,
    npar_per_block=(4,),
    cycle=0,
    time=0.0,
    redshift=None,
    scale_factor=None,
    must_refine=False,
):
    total_npar = sum(npar_per_block)
    int_names = ["pid", "blk"]
    if must_refine:
        int_names.append("must_refine")
    int_str = " ".join(int_names)

    header = (
        f"# Pasta output at global time = {time:.16e}\n"
        f"# coord_type = {coord_type}\n"
        f"# cycle = {cycle}\n"
    )
    if redshift is not None:
        header += f"# redshift = {redshift:.16e}\n"
    if scale_factor is not None:
        header += f"# scale_factor = {scale_factor:.16e}\n"
    header += (
        f"# particle count = {total_npar}\n"
        f"# int properties:  {int_str}\n"
        "# real properties: x1 x2 x3 v1 v2 v3 mass\n"
        f"# number of blocks = {nblocks}\n"
    )

    with open(filename, "wb") as f:
        f.write(np.int32(1).tobytes())
        f.write(header.encode("ascii"))
        f.write(np.array([len(int_names), 7, 8, 8], dtype="<i4").tobytes())
        f.write(np.int64(total_npar).tobytes())
        f.write(np.int64(nblocks).tobytes())
        f.write(np.array(npar_per_block, dtype="<i8").tobytes())

        # int properties: pid, blk, (optional must_refine)
        pid = np.arange(total_npar, dtype="<i8")
        blk_list = []
        for b, count in enumerate(npar_per_block):
            blk_list.extend([b] * count)
        blk = np.array(blk_list, dtype="<i8")
        f.write(pid.tobytes())
        f.write(blk.tobytes())
        if must_refine:
            mr = np.ones(total_npar, dtype="<i8")
            f.write(mr.tobytes())

        # real properties: x1, x2, x3, v1, v2, v3, mass
        x1 = np.linspace(0.1, 0.9, total_npar, dtype="<f8")
        x2 = np.linspace(0.2, 0.8, total_npar, dtype="<f8")
        x3 = np.linspace(0.3, 0.7, total_npar, dtype="<f8")
        v1 = np.ones(total_npar, dtype="<f8") * 10.0
        v2 = np.ones(total_npar, dtype="<f8") * 20.0
        v3 = np.ones(total_npar, dtype="<f8") * 30.0
        mass = np.ones(total_npar, dtype="<f8") * 0.05
        for arr in [x1, x2, x3, v1, v2, v3, mass]:
            f.write(arr.tobytes())


def test_pasta_1d_synthetic(tmp_path):
    fn = str(tmp_path / "data.prim.00000.binary")
    _write_pasta_binary(
        fn,
        coord_type="Cartesian",
        nblocks=2,
        ncells=(8, 1, 1),
        vars_to_write=("fluid_den_prim", "fluid_vel_1", "fluid_pressure"),
    )
    ds = yt.load(fn)
    assert isinstance(ds, PastaDataset)
    assert ds.dimensionality == 1
    assert ds.index.num_grids == 2
    ad = ds.all_data()
    assert ad["gas", "density"].size == 16
    assert ("gas", "velocity_x") in ds.derived_field_list
    assert ("gas", "pressure") in ds.derived_field_list


def test_pasta_3d_synthetic(tmp_path):
    fn = str(tmp_path / "data.prim.00000.binary")
    _write_pasta_binary(
        fn,
        coord_type="Cartesian",
        nblocks=2,
        ncells=(4, 4, 4),
        vars_to_write=("fluid_den_prim", "fluid_vel_1", "fluid_pressure"),
    )
    ds = yt.load(fn)
    assert isinstance(ds, PastaDataset)
    assert ds.dimensionality == 3
    assert ds.index.num_grids == 2
    ad = ds.all_data()
    assert ad["gas", "density"].size == 128


def test_pasta_cylindrical(tmp_path):
    fn = str(tmp_path / "data.prim.00000.binary")
    _write_pasta_binary(
        fn,
        coord_type="Cylindrical",
        nblocks=2,
        ncells=(4, 4, 4),
        x1_range=(0.1, 1.0),
        x2_range=(0.0, 2.0 * np.pi),
        x3_range=(-1.0, 1.0),
    )
    ds = yt.load(fn)
    assert isinstance(ds, PastaDataset)
    assert ds.geometry == "cylindrical"
    assert ds.coordinates.axis_order == ("r", "theta", "z")

    ad = ds.all_data()
    assert ("gas", "velocity_r") in ds.derived_field_list
    assert ("gas", "velocity_theta") in ds.derived_field_list
    assert ("gas", "velocity_z") in ds.derived_field_list
    assert ad["index", "r"].size == 2 * 64
    assert ad["index", "cell_volume"].size == 2 * 64


def test_pasta_spherical(tmp_path):
    fn = str(tmp_path / "data.prim.00000.binary")
    _write_pasta_binary(
        fn,
        coord_type="Spherical",
        nblocks=1,
        ncells=(4, 4, 4),
        x1_range=(0.1, 1.0),
        x2_range=(0.0, np.pi),
        x3_range=(0.0, 2.0 * np.pi),
    )
    ds = yt.load(fn)
    assert isinstance(ds, PastaDataset)
    assert ds.geometry == "spherical"
    assert ds.coordinates.axis_order == ("r", "theta", "phi")

    ad = ds.all_data()
    assert ("gas", "velocity_r") in ds.derived_field_list
    assert ("gas", "velocity_theta") in ds.derived_field_list
    assert ("gas", "velocity_phi") in ds.derived_field_list
    assert ad["index", "r"].size == 64
    assert ad["index", "cell_volume"].size == 64


def test_pasta_magnetic_fields(tmp_path):
    fn = str(tmp_path / "data.prim.00000.binary")
    _write_pasta_binary(
        fn,
        coord_type="Cartesian",
        nblocks=1,
        ncells=(4, 4, 4),
        vars_to_write=(
            "fluid_den_prim",
            "fluid_pressure",
            "field_bcc_1",
            "field_bcc_2",
            "field_bcc_3",
        ),
        num_face_field=3,
    )
    ds = yt.load(fn)
    assert ("gas", "magnetic_field_x") in ds.derived_field_list
    assert ("gas", "magnetic_field_y") in ds.derived_field_list
    assert ("gas", "magnetic_field_z") in ds.derived_field_list
    assert ("gas", "magnetic_field_strength") in ds.derived_field_list
    assert ("gas", "magnetic_energy_density") in ds.derived_field_list
    assert ("gas", "magnetic_pressure") in ds.derived_field_list
    assert ("gas", "plasma_beta") in ds.derived_field_list
    assert ("gas", "magnetic_field_face_x") in ds.derived_field_list

    ad = ds.all_data()
    bx = ad["gas", "magnetic_field_x"]
    assert str(bx.units) == "G"
    assert bx.size == 64
    bstr = ad["gas", "magnetic_field_strength"]
    assert str(bstr.units) == "G"
    pmag = ad["gas", "magnetic_pressure"]
    assert str(pmag.units) == "dyn/cm**2"
    beta = ad["gas", "plasma_beta"]
    assert beta.units.is_dimensionless

    # Read staggered face field
    bf1 = ad["pasta", "field_bf_1"]
    assert bf1.size == 64
    b_face_x = ad["gas", "magnetic_field_face_x"]
    assert str(b_face_x.units) == "G"


def test_pasta_magnetic_fields_face_only(tmp_path):
    fn = str(tmp_path / "data.prim.00000.binary")
    _write_pasta_binary(
        fn,
        coord_type="Cartesian",
        nblocks=1,
        ncells=(4, 4, 4),
        vars_to_write=("fluid_den_prim", "fluid_pressure"),
        num_face_field=3,
    )
    ds = yt.load(fn)
    assert ("gas", "magnetic_field_x") in ds.derived_field_list
    assert ("gas", "magnetic_field_strength") in ds.derived_field_list
    assert ("gas", "magnetic_pressure") in ds.derived_field_list
    assert ("gas", "plasma_beta") in ds.derived_field_list

    ad = ds.all_data()
    bx = ad["gas", "magnetic_field_x"]
    assert str(bx.units) == "G"
    assert bx.size == 64


def test_pasta_edge_emf_fields(tmp_path):
    fn = str(tmp_path / "data.prim.00000.binary")
    _write_pasta_binary(
        fn,
        coord_type="Cartesian",
        nblocks=1,
        ncells=(4, 4, 4),
        vars_to_write=("fluid_den_prim",),
        num_face_field=3,
        num_edge_field=3,
    )
    ds = yt.load(fn)
    assert ("pasta", "field_emf_1") in ds.field_list
    ad = ds.all_data()
    emf1 = ad["pasta", "field_emf_1"]
    assert emf1.size == 64


def test_pasta_particles(tmp_path):
    bin_dir = tmp_path / "bin"
    par_dir = tmp_path / "par"
    bin_dir.mkdir()
    par_dir.mkdir()

    grid_file = str(bin_dir / "data.prim.00001.binary")
    par_file = str(par_dir / "data.00001.par")

    _write_pasta_binary(
        grid_file,
        coord_type="Cartesian",
        nblocks=2,
        ncells=(4, 4, 4),
        cycle=1,
        time=0.1,
    )
    _write_pasta_particle(
        par_file,
        coord_type="Cartesian",
        nblocks=2,
        npar_per_block=(3, 5),
        cycle=1,
        time=0.1,
    )

    # 1. Load via grid binary file
    ds = yt.load(grid_file)
    assert ds._has_particles is True
    assert ds.index.num_grids == 2
    assert ds.index.grids[0].NumberOfParticles == 3
    assert ds.index.grids[1].NumberOfParticles == 5

    ad = ds.all_data()
    px = ad["io", "particle_position_x"]
    py = ad["io", "particle_position_y"]
    pz = ad["io", "particle_position_z"]
    mass = ad["io", "particle_mass"]
    pid = ad["io", "particle_index"]

    assert px.size == 8
    assert py.size == 8
    assert pz.size == 8
    assert mass.size == 8
    assert pid.size == 8
    assert ad["io", "particle_position"].shape == (8, 3)
    assert ad["io", "particle_velocity"].shape == (8, 3)

    # Spatial region selection
    reg = ds.region([0.25, 0.5, 0.5], [0.0, 0.0, 0.0], [0.5, 1.0, 1.0])
    # block 0 is in x [0, 0.5], block 1 is in [0.5, 1.0]
    reg_px = reg["io", "particle_position_x"]
    assert reg_px.size == 3

    # 2. Load via companion par file directly
    ds_par = yt.load(par_file)
    assert ds_par._has_particles is True
    assert ds_par.index.num_grids == 2
    ad_par = ds_par.all_data()
    assert ad_par["io", "particle_position_x"].size == 8

    # 3. Standalone par file (without companion grid)
    standalone_par = str(tmp_path / "standalone.par")
    shutil.copy(par_file, standalone_par)
    ds_standalone = yt.load(standalone_par)
    assert ds_standalone._has_particles is True
    assert ds_standalone.index.num_grids == 1
    ad_standalone = ds_standalone.all_data()
    assert ad_standalone["io", "particle_position_x"].size == 8


def test_pasta_conserved_thermodynamics(tmp_path):
    fn = str(tmp_path / "data.cons.00000.binary")
    _write_pasta_binary(
        fn,
        coord_type="Cartesian",
        nblocks=1,
        ncells=(4, 4, 4),
        vars_to_write=(
            "fluid_den_cons",
            "fluid_mom_1",
            "fluid_mom_2",
            "fluid_mom_3",
            "fluid_energy",
            "scalar_cons_0",
        ),
    )
    ds = yt.load(fn)
    assert ds.cosmological_simulation == 0
    assert ("gas", "density") in ds.derived_field_list
    assert ("gas", "momentum_density_x") in ds.derived_field_list
    assert ("gas", "specific_thermal_energy") in ds.derived_field_list
    assert ("gas", "pressure") in ds.derived_field_list
    assert ("gas", "sound_speed") in ds.derived_field_list
    assert ("gas", "temperature") in ds.derived_field_list

    ad = ds.all_data()
    assert ad["gas", "density"].size == 64
    assert ad["gas", "pressure"].size == 64
    assert ad["gas", "sound_speed"].size == 64
    assert ad["gas", "temperature"].size == 64
    assert ad["gas", "momentum_density_x"].size == 64


def test_pasta_multifile_and_cosmology(tmp_path):
    bin_dir = tmp_path / "bin"
    par_dir = tmp_path / "par"
    bin_dir.mkdir()
    par_dir.mkdir()

    fluid_file = str(bin_dir / "sim.fluid_prim.00000.binary")
    prho_file = str(bin_dir / "sim.particle_rho.00000.binary")
    scalar_file = str(bin_dir / "sim.scalar_prim.00000.binary")
    par_file = str(par_dir / "sim_par.00000.par")

    cgs_units = {
        "density": 1.0e-24,
        "length": 1.0e23,
        "time": 1.0e15,
        "velocity": 1.0e8,
        "pressure": 1.0e-8,
        "temperature": 1.0e4,
    }

    _write_pasta_binary(
        fluid_file,
        nblocks=2,
        ncells=(4, 4, 4),
        vars_to_write=(
            "fluid_den_prim",
            "fluid_vel_1",
            "fluid_vel_2",
            "fluid_vel_3",
            "fluid_pressure",
        ),
        cycle=0,
        time=0.5,
        redshift=20.0,
        scale_factor=1.0,
        cgs_units=cgs_units,
    )
    _write_pasta_binary(
        prho_file,
        nblocks=2,
        ncells=(4, 4, 4),
        vars_to_write=("particle_rho",),
        cycle=0,
        time=0.5,
        redshift=20.0,
        scale_factor=1.0,
        cgs_units=cgs_units,
    )
    _write_pasta_binary(
        scalar_file,
        nblocks=2,
        ncells=(4, 4, 4),
        vars_to_write=("scalar_prim_0", "scalar_prim_1"),
        cycle=0,
        time=0.5,
        redshift=20.0,
        scale_factor=1.0,
        cgs_units=cgs_units,
    )
    _write_pasta_particle(
        par_file,
        nblocks=2,
        npar_per_block=(4, 6),
        cycle=0,
        time=0.5,
        redshift=20.0,
        scale_factor=1.0,
        must_refine=True,
    )

    # 1. Load fluid file
    ds = yt.load(fluid_file)
    assert ds.cosmological_simulation == 1
    assert_allclose_units(ds.current_redshift, 20.0, 1e-4)
    assert ds._has_particles is True
    assert ("pasta", "particle_rho") in ds.field_list
    assert ("pasta", "scalar_prim_0") in ds.field_list
    assert ("gas", "dark_matter_density") in ds.derived_field_list
    assert ("gas", "H_fraction") in ds.derived_field_list

    ad = ds.all_data()
    assert ad["gas", "density"].size == 128
    assert ad["gas", "dark_matter_density"].size == 128
    assert ad["gas", "H_fraction"].size == 128
    assert len(ad["io", "particle_position_x"]) == 10
    assert len(ad["io", "must_refine"]) == 10

    # 2. Load par file directly
    ds_par = yt.load(par_file)
    assert ds_par.cosmological_simulation == 1
    assert ds_par._has_particles is True
    assert ds_par.index.num_grids == 2
    ad_par = ds_par.all_data()
    assert len(ad_par["io", "particle_position_x"]) == 10
    assert ad_par["gas", "density"].size == 128
    assert ad_par["gas", "dark_matter_density"].size == 128


def test_pasta_cosmology_standalone_no_files(tmp_path):
    fn = str(tmp_path / "standalone.fluid_prim.00000.binary")
    par_fn = str(tmp_path / "standalone_par.00000.par")
    cgs_units = {
        "length": 4.4081108299461778e23,
        "density": 2.7611621778385654e-24,
        "time": 657122364540844.5,
        "velocity": 670820393.24993706,
        "pressure": 1.242522980027355e-06,
        "temperature": 6602952298.7890501,
    }
    _write_pasta_binary(
        fn,
        nblocks=1,
        ncells=(4, 4, 4),
        vars_to_write=("fluid_den_prim", "fluid_pressure"),
        cycle=0,
        time=0.816,
        redshift=99.0,
        scale_factor=1.0,
        cgs_units=cgs_units,
    )
    _write_pasta_particle(
        par_fn,
        nblocks=1,
        npar_per_block=(8,),
        cycle=0,
        time=0.816,
        redshift=99.0,
        scale_factor=1.0,
        must_refine=True,
    )
    # Load dataset with NO companion halo catalog, log, or input files
    ds = yt.load(fn)
    assert ds.cosmological_simulation == 1
    assert_allclose_units(ds.current_redshift, 99.0, 1e-4)
    assert "Mpccm/h" in str(ds.length_unit.units)
    assert_allclose_units(ds.parameters["comoving_box_size"], 10.0, 1e-2)
    assert ds._has_particles is True

    ad = ds.all_data()
    assert ad["gas", "density"].size == 64
    assert len(ad["io", "particle_position_x"]) == 8
    assert len(ad["io", "particle_velocity_x"]) == 8
    assert len(ad["io", "particle_mass"]) == 8
    assert len(ad["io", "must_refine"]) == 8
    assert ad["io", "particle_position"].shape == (8, 3)

    # Also load the companion particle file directly
    ds_par = yt.load(par_fn)
    assert ds_par.cosmological_simulation == 1
    assert_allclose_units(ds_par.current_redshift, 99.0, 1e-4)
    assert ds_par._has_particles is True
    ad_par = ds_par.all_data()
    assert len(ad_par["io", "particle_position_x"]) == 8
    assert ad_par["gas", "density"].size == 64


def test_pasta_cosmology_user_parameters(tmp_path):
    fn = str(tmp_path / "custom.fluid_prim.00000.binary")
    _write_pasta_binary(
        fn,
        nblocks=1,
        ncells=(4, 4, 4),
        vars_to_write=("fluid_den_prim",),
        redshift=15.0,
        scale_factor=1.0,
        cgs_units={
            "length": 1e24,
            "density": 1e-25,
            "time": 1e15,
            "velocity": 1e9,
            "pressure": 1e-7,
            "temperature": 1e4,
        },
    )
    ds = yt.load(
        fn,
        parameters={
            "omega_matter": 0.27,
            "omega_lambda": 0.73,
            "hubble_constant": 0.71,
            "comoving_box_size": 25.0,
        },
    )
    assert ds.cosmological_simulation == 1
    assert ds.omega_matter == 0.27
    assert ds.omega_lambda == 0.73
    assert ds.hubble_constant == 0.71
    assert ds.parameters["comoving_box_size"] == 25.0


def test_pasta_primordial_chemistry_6species(tmp_path):
    fn = str(tmp_path / "chem6.fluid_prim.00000.binary")
    vars_to_write = [
        "fluid_den_prim",
        "fluid_pressure",
        "scalar_cons_0",
        "scalar_cons_1",
        "scalar_cons_2",
        "scalar_cons_3",
        "scalar_cons_4",
        "scalar_cons_5",
        "scalar_cons_6",
    ]
    _write_pasta_binary(fn, vars_to_write=vars_to_write)
    ds = yt.load(fn, parameters={"primordial_chemistry": 1})
    assert ds.primordial_chemistry == 1
    assert ds.species_names == [
        "H_p0",
        "H_p1",
        "He_p0",
        "He_p1",
        "He_p2",
        "El",
    ]

    for fld in (
        ("gas", "H_p0_density"),
        ("gas", "H_density"),
        ("gas", "H_p0_fraction"),
        ("gas", "H_fraction"),
        ("gas", "H_p1_density"),
        ("gas", "H_p1_fraction"),
        ("gas", "He_p0_density"),
        ("gas", "He_density"),
        ("gas", "He_p0_fraction"),
        ("gas", "He_fraction"),
        ("gas", "He_p1_density"),
        ("gas", "He_p1_fraction"),
        ("gas", "He_p2_density"),
        ("gas", "He_p2_fraction"),
        ("gas", "El_density"),
        ("gas", "El_fraction"),
        ("gas", "El_number_density"),
        ("gas", "metal_density"),
        ("gas", "metallicity"),
        ("gas", "metal_fraction"),
        ("gas", "H_nuclei_density"),
        ("gas", "He_nuclei_density"),
    ):
        assert fld in ds.derived_field_list

    # 9- and 12-species specific fields should not be present
    assert ("gas", "H2_p0_density") not in ds.derived_field_list
    assert ("gas", "H_m1_density") not in ds.derived_field_list
    assert ("gas", "D_p0_density") not in ds.derived_field_list

    ad = ds.all_data()
    assert ad["gas", "H_density"].size == 512
    assert ad["gas", "H_fraction"].size == 512
    assert ad["gas", "El_density"].size == 512
    assert ad["gas", "El_number_density"].size == 512
    assert ad["gas", "metal_density"].size == 512
    assert ad["gas", "metallicity"].size == 512


def test_pasta_primordial_chemistry_9species(tmp_path):
    fn = str(tmp_path / "chem9.fluid_prim.00000.binary")
    vars_to_write = ["fluid_den_prim", "fluid_pressure"] + [
        f"scalar_cons_{i}" for i in range(10)
    ]
    _write_pasta_binary(fn, vars_to_write=vars_to_write)
    ds = yt.load(fn, parameters={"primordial_chemistry": 2})
    assert ds.primordial_chemistry == 2
    assert ds.species_names == [
        "H_p0",
        "H_p1",
        "He_p0",
        "He_p1",
        "He_p2",
        "El",
        "H_m1",
        "H2_p0",
        "H2_p1",
    ]

    for fld in (
        ("gas", "H_m1_density"),
        ("gas", "H_m1_fraction"),
        ("gas", "H_m1_number_density"),
        ("gas", "H2_p0_density"),
        ("gas", "H2_density"),
        ("gas", "H2_p0_fraction"),
        ("gas", "H2_fraction"),
        ("gas", "H2_number_density"),
        ("gas", "H2_p1_density"),
        ("gas", "H2_p1_fraction"),
        ("gas", "metal_density"),
        ("gas", "metallicity"),
    ):
        assert fld in ds.derived_field_list

    # 12-species specific fields should not be present
    assert ("gas", "D_p0_density") not in ds.derived_field_list
    assert ("gas", "HD_p0_density") not in ds.derived_field_list

    ad = ds.all_data()
    assert ad["gas", "H2_density"].size == 512
    assert ad["gas", "H_m1_density"].size == 512
    assert ad["gas", "metal_density"].size == 512


def test_pasta_primordial_chemistry_12species(tmp_path):
    fn = str(tmp_path / "chem12.fluid_prim.00000.binary")
    vars_to_write = ["fluid_den_prim", "fluid_pressure"] + [
        f"scalar_cons_{i}" for i in range(13)
    ]
    _write_pasta_binary(fn, vars_to_write=vars_to_write)
    ds = yt.load(fn, parameters={"primordial_chemistry": 3})
    assert ds.primordial_chemistry == 3
    assert ds.species_names == [
        "H_p0",
        "H_p1",
        "He_p0",
        "He_p1",
        "He_p2",
        "El",
        "H_m1",
        "H2_p0",
        "H2_p1",
        "D_p0",
        "D_p1",
        "HD_p0",
    ]

    for fld in (
        ("gas", "D_p0_density"),
        ("gas", "D_density"),
        ("gas", "D_p0_fraction"),
        ("gas", "D_fraction"),
        ("gas", "D_number_density"),
        ("gas", "D_p1_density"),
        ("gas", "D_p1_fraction"),
        ("gas", "HD_p0_density"),
        ("gas", "HD_density"),
        ("gas", "HD_p0_fraction"),
        ("gas", "HD_fraction"),
        ("gas", "HD_number_density"),
        ("gas", "D_nuclei_density"),
        ("gas", "metal_density"),
        ("gas", "metallicity"),
    ):
        assert fld in ds.derived_field_list

    ad = ds.all_data()
    assert ad["gas", "D_density"].size == 512
    assert ad["gas", "HD_density"].size == 512
    assert ad["gas", "D_nuclei_density"].size == 512
    assert ad["gas", "metal_density"].size == 512


def test_pasta_primordial_chemistry_primitive_scalars(tmp_path):
    fn = str(tmp_path / "chem_prim.scalar_prim.00000.binary")
    fluid_fn = str(tmp_path / "chem_prim.fluid_prim.00000.binary")
    _write_pasta_binary(
        fluid_fn,
        vars_to_write=("fluid_den_prim", "fluid_pressure"),
    )
    _write_pasta_binary(
        fn,
        vars_to_write=[f"scalar_prim_{i}" for i in range(6)],
    )
    ds = yt.load(fluid_fn, parameters={"primordial_chemistry": 1})
    assert ds.primordial_chemistry == 1
    assert ("gas", "H_p0_fraction") in ds.derived_field_list
    assert ("gas", "H_fraction") in ds.derived_field_list
    assert ("gas", "H_p0_density") in ds.derived_field_list
    assert ("gas", "H_density") in ds.derived_field_list
    ad = ds.all_data()
    assert ad["gas", "H_density"].size == 512
    assert ad["gas", "H_fraction"].size == 512
    assert ad["gas", "El_density"].size == 512
