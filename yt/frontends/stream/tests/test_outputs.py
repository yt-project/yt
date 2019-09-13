"""
Title: test_stream.py
Purpose: stream frontend answer tests
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import os
import tempfile

import numpy as np
import pytest

import yt
from yt.convenience import load
from yt.data_objects.profiles import create_profile
from yt.frontends.stream.api import \
    load_amr_grids, \
    load_particles, \
    hexahedral_connectivity, \
    load_hexahedral_mesh, \
    load_uniform_grid, \
    refine_amr
from yt.testing import \
    assert_almost_equal, \
    assert_equal, \
    assert_raises, \
    fake_random_ds
from yt.utilities.exceptions import \
    YTOutputNotIdentified, \
    YTInconsistentGridFieldShape, \
    YTInconsistentParticleFieldShape, \
    YTInconsistentGridFieldShapeGridDims, \
    YTIllDefinedAMR, \
    YTIntDomainOverflow
import yt.utilities.flagging_methods as fm
import yt.utilities.initial_conditions as ic
import yt.utilities.answer_testing.framework as fw


# Globals
OCT_MASK_LIST = [8, 0, 0, 0, 0, 8, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0,
                 8, 0, 0, 0, 0, 0, 0, 0,
                 0]


#============================================
#               TestStream
#============================================
@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
class TestStream(fw.AnswerTest):
    @pytest.mark.usefixtures('temp_dir')
    def test_load_empty_file(self):
        assert_raises(YTOutputNotIdentified, load, "not_a_file")
        assert_raises(YTOutputNotIdentified, load, tempfile.mkstemp('empty_file', dir=os.getcwd()))
        assert_raises(YTOutputNotIdentified, load, tempfile.mkdtemp(dir=os.getcwd()))

    #-----
    # test_dimensionless_field_units
    #-----
    def test_dimensionless_field_units(self):
        Z = np.random.uniform(size=(32,32,32))
        d = np.random.uniform(size=(32,32,32))
        data = {"density": d, "metallicity": Z}
        ds = load_uniform_grid(data, (32,32,32))
        dd = ds.all_data()
        assert_equal(Z.max(), dd["metallicity"].max())

    #-----
    # test_inconsistent_field_shape
    #-----
    def test_inconsistent_field_shape(self):
        def load_field_field_mismatch():
            d = np.random.uniform(size=(32, 32, 32))
            t = np.random.uniform(size=(32, 64, 32))
            data = {"density": d, "temperature": t}
            load_uniform_grid(data, (32,32,32))
        def load_field_grid_mismatch():
            d = np.random.uniform(size=(32, 32, 32))
            t = np.random.uniform(size=(32, 32, 32))
            data = {"density": d, "temperature": t}
            load_uniform_grid(data, (32,64,32))
        def load_particle_fields_mismatch():
            x = np.random.uniform(size=100)
            y = np.random.uniform(size=100)
            z = np.random.uniform(size=200)
            data = {"particle_position_x": x,
                    "particle_position_y": y,
                    "particle_position_z": z}
            load_particles(data)
        assert_raises(YTInconsistentGridFieldShape,
                      load_field_field_mismatch)
        assert_raises(YTInconsistentGridFieldShapeGridDims,
                      load_field_grid_mismatch)
        assert_raises(YTInconsistentParticleFieldShape,
                      load_particle_fields_mismatch)

    #-----
    # test_qt_overflow
    #-----
    def test_qt_overflow(self):
        grid_data = []
        grid_dict = {}
        grid_dict['left_edge'] = [-1.0, -1.0, -1.0]
        grid_dict['right_edge'] = [1.0, 1.0, 1.0]
        grid_dict['dimensions'] = [8, 8, 8]
        grid_dict['level'] = 0
        grid_dict['density'] = np.ones((8,8,8))
        grid_data.append(grid_dict)
        domain_dimensions = np.array([8, 8, 8])
        spf = yt.load_amr_grids(grid_data, domain_dimensions)
        def make_proj():
            p = yt.ProjectionPlot(spf, 'x', ["density"], center='c', origin='native')
            return p
        assert_raises(YTIntDomainOverflow, make_proj)

    #-----
    # test_refine_by
    #-----
    def test_refine_by(self):
        grid_data = []
        ref_by = 4
        lo = 0.0
        hi = 1.0
        fine_grid_width = (hi - lo) / ref_by
        for level in range(2):
            grid_dict = {}
            grid_dict['left_edge'] = [0.0 + 0.5*fine_grid_width*level]*3
            grid_dict['right_edge'] = [1.0 - 0.5*fine_grid_width*level]*3
            grid_dict['dimensions'] = [8, 8, 8]
            grid_dict['level'] = level
            grid_dict['density'] = np.ones((8,8,8))
            grid_data.append(grid_dict)
        domain_dimensions = np.array([8, 8, 8])
        yt.load_amr_grids(grid_data, domain_dimensions, refine_by=ref_by)

    #-----
    # test_validation
    #-----
    def test_validation(self):
        dims = np.array([4, 2, 4])
        grid_data = [
            dict(left_edge = [0.0, 0.0, 0.0],
                 right_edge = [1.0, 1.0, 1.0],
                 level = 0,
                 dimensions = dims),
            dict(left_edge = [0.25, 0.25, 0.25],
                 right_edge = [0.75, 0.75, 0.75],
                 level = 1,
                 dimensions = dims),
           ]
        bbox = np.array([[0, 1], [0, 1], [0, 1]])
        def load_grids():
            yt.load_amr_grids(grid_data, dims, bbox=bbox, periodicity=(0, 0, 0),
                           length_unit=1.0, refine_by=2)
        assert_raises(YTIllDefinedAMR, load_grids)

    #-----
    # test_stream_hexahedral
    #-----
    def test_stream_hexahedral(self):
        np.random.seed(0x4d3d3d3)
        Nx, Ny, Nz = 32, 18, 24
        # Note what we're doing here -- we are creating a randomly spaced mesh, but
        # because of how the accumulate operation works, we also reset the leftmost
        # cell boundary to 0.0.
        cell_x = np.random.random(Nx+1)
        cell_x /= cell_x.sum()
        cell_x = np.add.accumulate(cell_x)
        cell_x[0] = 0.0
        cell_y = np.random.random(Ny+1)
        cell_y /= cell_y.sum()
        cell_y = np.add.accumulate(cell_y)
        cell_y[0] = 0.0
        cell_z = np.random.random(Nz+1)
        cell_z /= cell_z.sum()
        cell_z = np.add.accumulate(cell_z)
        cell_z[0] = 0.0
        coords, conn = hexahedral_connectivity(cell_x, cell_y, cell_z)
        data = {'random_field': np.random.random((Nx, Ny, Nz))}
        bbox = np.array([ [0.0, 1.0], [0.0, 1.0], [0.0, 1.0] ])
        ds = load_hexahedral_mesh(data, conn, coords, bbox=bbox)
        dd = ds.all_data()
        #raise RuntimeError
        assert_almost_equal(float(dd["cell_volume"].sum(dtype="float64")), 1.0)
        assert_equal(dd["ones"].size, Nx * Ny * Nz)
        # Now we try it with a standard mesh
        cell_x = np.linspace(0.0, 1.0, Nx+1)
        cell_y = np.linspace(0.0, 1.0, Ny+1)
        cell_z = np.linspace(0.0, 1.0, Nz+1)
        coords, conn = hexahedral_connectivity(cell_x, cell_y, cell_z)
        data = {'random_field': np.random.random((Nx, Ny, Nz))}
        bbox = np.array([ [0.0, 1.0], [0.0, 1.0], [0.0, 1.0] ])
        ds = load_hexahedral_mesh(data, conn, coords, bbox=bbox)
        dd = ds.all_data()
        assert_almost_equal(float(dd["cell_volume"].sum(dtype="float64")), 1.0)
        assert_equal(dd["ones"].size, Nx * Ny * Nz)
        assert_almost_equal(dd["dx"].to_ndarray(), 1.0/Nx)
        assert_almost_equal(dd["dy"].to_ndarray(), 1.0/Ny)
        assert_almost_equal(dd["dz"].to_ndarray(), 1.0/Nz)
        s = yt.SlicePlot(ds, "x", "random_field")
        s._setup_plots()
        s.frb["random_field"]

    #-----
    # test_octree
    #-----
    def test_octree(self):
        # See Issue #1272
        octree_mask = np.array(OCT_MASK_LIST, dtype=np.uint8)
        quantities = {}
        quantities[('gas', 'density')] = np.ones((22, 1), dtype=float)
        bbox = np.array([[-10., 10.], [-10., 10.], [-10., 10.]])
        ds = yt.load_octree(octree_mask=octree_mask,
                            data=quantities,
                            bbox=bbox,
                            over_refine_factor=0,
                            partial_coverage=0)
        proj = ds.proj('density', 'x')
        proj['density']

    #-----
    # test_stream_particles
    #-----
    def test_stream_particles(self):
        num_particles = 100000
        domain_dims = (64, 64, 64)
        dens = np.random.random(domain_dims) 
        x = np.random.uniform(size=num_particles)
        y = np.random.uniform(size=num_particles)
        z = np.random.uniform(size=num_particles)
        m = np.ones(num_particles)
        # Field operators and cell flagging methods
        fo = []
        fo.append(ic.TopHatSphere(0.1, [0.2,0.3,0.4],{"density": 2.0}))
        fo.append(ic.TopHatSphere(0.05, [0.7,0.4,0.75],{"density": 20.0}))
        rc = [fm.flagging_method_registry["overdensity"](1.0)]
        # Check that all of this runs ok without particles
        ug0 = load_uniform_grid({"density": dens}, domain_dims, 1.0, nprocs=8)
        amr0 = refine_amr(ug0, rc, fo, 3)
        grid_data = []
        for grid in amr0.index.grids:
            data = dict(left_edge=grid.LeftEdge,
                        right_edge=grid.RightEdge,
                        level=grid.Level,
                        dimensions=grid.ActiveDimensions)
            for field in amr0.field_list:
                data[field] = grid[field]
            grid_data.append(data)
        amr0 = load_amr_grids(grid_data, domain_dims)
        # Now add particles
        fields1 = {"density": dens,
                   "particle_position_x": x,
                   "particle_position_y": y,
                   "particle_position_z": z,
                   "particle_mass": m}
        fields2 = fields1.copy()
        ug1 = load_uniform_grid(fields1, domain_dims, 1.0)
        ug2 = load_uniform_grid(fields2, domain_dims, 1.0, nprocs=8)
        # Check to make sure the number of particles is the same
        number_of_particles1 = np.sum([grid.NumberOfParticles for grid in ug1.index.grids])
        number_of_particles2 = np.sum([grid.NumberOfParticles for grid in ug2.index.grids])
        assert_equal(number_of_particles1, num_particles)
        assert_equal(number_of_particles1, number_of_particles2)
        for grid in ug2.index.grids:
            tot_parts = grid["io","particle_position_x"].size
            tot_all_parts = grid["all","particle_position_x"].size
            assert tot_parts == grid.NumberOfParticles
            assert tot_all_parts == grid.NumberOfParticles
        # Check to make sure the fields have been defined correctly
        for ptype in ("all", "io"):
            assert ug1._get_field_info(ptype, "particle_position_x").particle_type
            assert ug1._get_field_info(ptype, "particle_position_y").particle_type
            assert ug1._get_field_info(ptype, "particle_position_z").particle_type
            assert ug1._get_field_info(ptype, "particle_mass").particle_type
        assert not ug1._get_field_info("gas", "density").particle_type
        for ptype in ("all", "io"):
            assert ug2._get_field_info(ptype, "particle_position_x").particle_type
            assert ug2._get_field_info(ptype, "particle_position_y").particle_type
            assert ug2._get_field_info(ptype, "particle_position_z").particle_type
            assert ug2._get_field_info(ptype, "particle_mass").particle_type
        assert not ug2._get_field_info("gas", "density").particle_type
        # Now refine this
        amr1 = refine_amr(ug1, rc, fo, 3)
        for field in sorted(ug1.field_list):
            assert field in amr1.field_list
        grid_data = []
        for grid in amr1.index.grids:
            data = dict(left_edge=grid.LeftEdge,
                        right_edge=grid.RightEdge,
                        level=grid.Level,
                        dimensions=grid.ActiveDimensions)
            for field in amr1.field_list:
                if field[0] != "all":
                    data[field] = grid[field]
            grid_data.append(data)
        amr2 = load_amr_grids(grid_data, domain_dims)
        # Check everything again
        number_of_particles1 = [grid.NumberOfParticles for grid in amr1.index.grids]
        number_of_particles2 = [grid.NumberOfParticles for grid in amr2.index.grids]
        assert_equal(np.sum(number_of_particles1), num_particles)
        assert_equal(number_of_particles1, number_of_particles2)
        for grid in amr1.index.grids:
            tot_parts = grid["io", "particle_position_x"].size
            tot_all_parts = grid["all", "particle_position_x"].size
            assert tot_parts == grid.NumberOfParticles
            assert tot_all_parts == grid.NumberOfParticles
        for grid in amr2.index.grids:
            tot_parts = grid["io", "particle_position_x"].size
            tot_all_parts = grid["all", "particle_position_x"].size
            assert tot_parts == grid.NumberOfParticles
            assert tot_all_parts == grid.NumberOfParticles
        assert amr1._get_field_info("all", "particle_position_x").particle_type
        assert amr1._get_field_info("all", "particle_position_y").particle_type
        assert amr1._get_field_info("all", "particle_position_z").particle_type
        assert amr1._get_field_info("all", "particle_mass").particle_type
        assert not amr1._get_field_info("gas", "density").particle_type
        assert amr2._get_field_info("all", "particle_position_x").particle_type
        assert amr2._get_field_info("all", "particle_position_y").particle_type
        assert amr2._get_field_info("all", "particle_position_z").particle_type
        assert amr2._get_field_info("all", "particle_mass").particle_type
        assert not amr2._get_field_info("gas", "density").particle_type
        # Now perform similar checks, but with multiple particle types
        num_dm_particles = 30000
        xd = np.random.uniform(size=num_dm_particles)
        yd = np.random.uniform(size=num_dm_particles)
        zd = np.random.uniform(size=num_dm_particles)
        md = np.ones(num_dm_particles)
        num_star_particles = 20000
        xs = np.random.uniform(size=num_star_particles)
        ys = np.random.uniform(size=num_star_particles)
        zs = np.random.uniform(size=num_star_particles)
        ms = 2.0*np.ones(num_star_particles)
        dens = np.random.random(domain_dims)
        fields3 = {"density": dens,
                   ("dm", "particle_position_x"): xd,
                   ("dm", "particle_position_y"): yd,
                   ("dm", "particle_position_z"): zd,
                   ("dm", "particle_mass"): md,
                   ("star", "particle_position_x"): xs,
                   ("star", "particle_position_y"): ys,
                   ("star", "particle_position_z"): zs,
                   ("star", "particle_mass"): ms}
        fields4 = fields3.copy()
        ug3 = load_uniform_grid(fields3, domain_dims, 1.0)
        ug4 = load_uniform_grid(fields4, domain_dims, 1.0, nprocs=8)
        # Check to make sure the number of particles is the same
        number_of_particles3 = np.sum([grid.NumberOfParticles for grid in ug3.index.grids])
        number_of_particles4 = np.sum([grid.NumberOfParticles for grid in ug4.index.grids])
        assert_equal(number_of_particles3, num_dm_particles+num_star_particles)
        assert_equal(number_of_particles3, number_of_particles4)
        for grid in ug4.index.grids:
            tot_parts = grid["dm", "particle_position_x"].size
            tot_parts += grid["star", "particle_position_x"].size
            tot_all_parts = grid["all", "particle_position_x"].size
            assert tot_parts == grid.NumberOfParticles
            assert tot_all_parts == grid.NumberOfParticles
        # Check to make sure the fields have been defined correctly
        for ptype in ("dm", "star"):
            assert ug3._get_field_info(ptype, "particle_position_x").particle_type
            assert ug3._get_field_info(ptype, "particle_position_y").particle_type
            assert ug3._get_field_info(ptype, "particle_position_z").particle_type
            assert ug3._get_field_info(ptype, "particle_mass").particle_type
            assert ug4._get_field_info(ptype, "particle_position_x").particle_type
            assert ug4._get_field_info(ptype, "particle_position_y").particle_type
            assert ug4._get_field_info(ptype, "particle_position_z").particle_type
            assert ug4._get_field_info(ptype, "particle_mass").particle_type
        # Now refine this
        amr3 = refine_amr(ug3, rc, fo, 3)
        for field in sorted(ug3.field_list):
            assert field in amr3.field_list
        grid_data = []
        for grid in amr3.index.grids:
            data = dict(left_edge=grid.LeftEdge,
                        right_edge=grid.RightEdge,
                        level=grid.Level,
                        dimensions=grid.ActiveDimensions)
            for field in amr3.field_list:
                if field[0] != "all":
                    data[field] = grid[field]
            grid_data.append(data)
        amr4 = load_amr_grids(grid_data, domain_dims)
        # Check everything again
        number_of_particles3 = [grid.NumberOfParticles for grid in amr3.index.grids]
        number_of_particles4 = [grid.NumberOfParticles for grid in amr4.index.grids]
        assert_equal(np.sum(number_of_particles3), num_star_particles+num_dm_particles)
        assert_equal(number_of_particles3, number_of_particles4)
        for ptype in ("dm", "star"):
            assert amr3._get_field_info(ptype, "particle_position_x").particle_type
            assert amr3._get_field_info(ptype, "particle_position_y").particle_type
            assert amr3._get_field_info(ptype, "particle_position_z").particle_type
            assert amr3._get_field_info(ptype, "particle_mass").particle_type
            assert amr4._get_field_info(ptype, "particle_position_x").particle_type
            assert amr4._get_field_info(ptype, "particle_position_y").particle_type
            assert amr4._get_field_info(ptype, "particle_position_z").particle_type
            assert amr4._get_field_info(ptype, "particle_mass").particle_type
        for grid in amr3.index.grids:
            tot_parts = grid["dm", "particle_position_x"].size
            tot_parts += grid["star", "particle_position_x"].size
            tot_all_parts = grid["all", "particle_position_x"].size
            assert tot_parts == grid.NumberOfParticles
            assert tot_all_parts == grid.NumberOfParticles
        for grid in amr4.index.grids:
            tot_parts = grid["dm", "particle_position_x"].size
            tot_parts += grid["star", "particle_position_x"].size
            tot_all_parts = grid["all", "particle_position_x"].size
            assert tot_parts == grid.NumberOfParticles
            assert tot_all_parts == grid.NumberOfParticles

    #-----
    # test_load_particles_types
    #-----
    def test_load_particles_types(self):
        num_particles = 10000
        data1 = {"particle_position_x": np.random.random(size=num_particles),
                 "particle_position_y": np.random.random(size=num_particles),
                 "particle_position_z": np.random.random(size=num_particles),
                 "particle_mass": np.ones(num_particles)}
        ds1 = load_particles(data1)
        ds1.index
        assert set(ds1.particle_types) == {"all", "io"}
        dd = ds1.all_data()
        for ax in "xyz":
            assert dd["io", "particle_position_%s" % ax].size == num_particles
            assert dd["all", "particle_position_%s" % ax].size == num_particles
        num_dm_particles = 10000
        num_star_particles = 50000
        num_tot_particles = num_dm_particles + num_star_particles
        data2 = {("dm", "particle_position_x"): np.random.random(size=num_dm_particles),
                 ("dm", "particle_position_y"): np.random.random(size=num_dm_particles),
                 ("dm", "particle_position_z"): np.random.random(size=num_dm_particles),
                 ("dm", "particle_mass"): np.ones(num_dm_particles),
                 ("star", "particle_position_x"): np.random.random(size=num_star_particles),
                 ("star", "particle_position_y"): np.random.random(size=num_star_particles),
                 ("star", "particle_position_z"): np.random.random(size=num_star_particles),
                 ("star", "particle_mass"): 2.0*np.ones(num_star_particles)}
        ds2 = load_particles(data2)
        ds2.index
        assert set(ds2.particle_types) == {"all", "star", "dm"}
        dd = ds2.all_data()
        for ax in "xyz":
            npart = 0
            for ptype in ds2.particle_types_raw:
                npart += dd[ptype, "particle_position_%s" % ax].size
            assert npart == num_tot_particles
            assert dd["all", "particle_position_%s" % ax].size == num_tot_particles

    #-----
    # test_particles_outside_domain
    #-----
    def test_particles_outside_domain(self):
        np.random.seed(0x4d3d3d3)
        posx_arr = np.random.uniform(low=-1.6, high=1.5, size=1000)
        posy_arr = np.random.uniform(low=-1.5, high=1.5, size=1000)
        posz_arr = np.random.uniform(low=-1.5, high=1.5, size=1000)
        dens_arr = np.random.random((16, 16, 16))
        data = dict(
            density=dens_arr,
            particle_position_x=posx_arr,
            particle_position_y=posy_arr,
            particle_position_z=posz_arr)
        bbox = np.array([[-1.5, 1.5], [-1.5, 1.5], [-1.5, 1.5]])
        ds = load_uniform_grid(data, (16, 16, 16), bbox=bbox, nprocs=4)
        wh = (posx_arr < bbox[0, 0]).nonzero()[0]
        assert wh.size == 1000 - ds.particle_type_counts['io']
        ad = ds.all_data()
        assert ds.particle_type_counts['io'] == ad['particle_position_x'].size

    #-----
    # test_multi_mesh
    #-----
    def test_multi_mesh(self):
        coordsMulti = np.array([[0.0, 0.0],
                                [1.0, 0.0],
                                [1.0, 1.0],
                                [0.0, 1.0]], dtype=np.float64)
        connect1 = np.array([[0, 1, 3], ], dtype=np.int64)
        connect2 = np.array([[1, 2, 3], ], dtype=np.int64)
        data1 = {}
        data2 = {}
        data1['connect1', 'test'] = np.array([[0.0, 1.0, 3.0], ], dtype=np.float64)
        data2['connect2', 'test'] = np.array([[1.0, 2.0, 3.0], ], dtype=np.float64)
        connectList = [connect1, connect2]
        dataList    = [data1, data2]
        ds = yt.load_unstructured_mesh(connectList, coordsMulti, dataList)
        sl = yt.SlicePlot(ds, 'z', ('connect1', 'test'))
        sl = yt.SlicePlot(ds, 'z', ('connect2', 'test'))
        sl = yt.SlicePlot(ds, 'z', ('all', 'test'))
        sl.annotate_mesh_lines()

    #-----
    # test_multi_field
    #-----
    def test_multi_field(self):
        coords = np.array([[0.0, 0.0],
                           [1.0, 0.0],
                           [1.0, 1.0],
                           [0.0, 1.0]], dtype=np.float64)
        connect = np.array([[0, 1, 3],
                            [1, 2, 3]], dtype=np.int64)
        data = {}
        data['connect1', 'test']      = np.array([[0.0, 1.0, 3.0],
                                                  [1.0, 2.0, 3.0]], dtype=np.float64)
        data['connect1', 'testAgain'] = np.array([[0.0, 1.0, 3.0],
                                                  [1.0, 2.0, 3.0]], dtype=np.float64)
        ds = yt.load_unstructured_mesh(connect, coords, data)
        sl = yt.SlicePlot(ds, 'z', 'test')
        sl.annotate_mesh_lines()
        sl = yt.SlicePlot(ds, 'z', 'testAgain')
        sl.annotate_mesh_lines()

    #-----
    # test_update_data
    #-----
    def test_update_data(self):
        ds = fake_random_ds(64, nprocs=8)
        ds.index
        dims = (32,32,32)
        grid_data = [{"temperature": np.random.uniform(size=dims)}
                     for i in range(ds.index.num_grids)]
        ds.index.update_data(grid_data)
        prj = ds.proj("temperature", 2)
        prj["temperature"]
        dd = ds.all_data()
        profile = create_profile(dd, "density", "temperature", 10)
        profile["temperature"]
