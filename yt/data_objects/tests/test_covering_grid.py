import numpy as np

from yt import \
    load
from yt.frontends.stream.data_structures import load_particles
from yt.testing import \
    requires_file, \
    fake_random_ds, \
    assert_equal, \
    assert_almost_equal

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

def test_covering_grid():
    # We decompose in different ways
    for level in [0, 1, 2]:
        for nprocs in [1, 2, 4, 8]:
            ds = fake_random_ds(16, nprocs = nprocs)
            dn = ds.refine_by**level 
            cg = ds.covering_grid(level, [0.0, 0.0, 0.0],
                    dn * ds.domain_dimensions)
            # Test coordinate generation
            yield assert_equal, np.unique(cg["dx"]).size, 1
            xmi = cg["x"].min()
            xma = cg["x"].max()
            dx = cg["dx"].flat[0:1]
            edges = ds.arr([[0,1],[0,1],[0,1]], 'code_length')
            yield assert_equal, xmi, edges[0,0] + dx/2.0
            yield assert_equal, xmi, cg["x"][0,0,0]
            yield assert_equal, xmi, cg["x"][0,1,1]
            yield assert_equal, xma, edges[0,1] - dx/2.0
            yield assert_equal, xma, cg["x"][-1,0,0]
            yield assert_equal, xma, cg["x"][-1,1,1]
            yield assert_equal, np.unique(cg["dy"]).size, 1
            ymi = cg["y"].min()
            yma = cg["y"].max()
            dy = cg["dy"][0]
            yield assert_equal, ymi, edges[1,0] + dy/2.0
            yield assert_equal, ymi, cg["y"][0,0,0]
            yield assert_equal, ymi, cg["y"][1,0,1]
            yield assert_equal, yma, edges[1,1] - dy/2.0
            yield assert_equal, yma, cg["y"][0,-1,0]
            yield assert_equal, yma, cg["y"][1,-1,1]
            yield assert_equal, np.unique(cg["dz"]).size, 1
            zmi = cg["z"].min()
            zma = cg["z"].max()
            dz = cg["dz"][0]
            yield assert_equal, zmi, edges[2,0] + dz/2.0
            yield assert_equal, zmi, cg["z"][0,0,0]
            yield assert_equal, zmi, cg["z"][1,1,0]
            yield assert_equal, zma, edges[2,1] - dz/2.0
            yield assert_equal, zma, cg["z"][0,0,-1]
            yield assert_equal, zma, cg["z"][1,1,-1]
            # Now we test other attributes
            yield assert_equal, cg["ones"].max(), 1.0
            yield assert_equal, cg["ones"].min(), 1.0
            yield assert_equal, cg["grid_level"], level
            yield assert_equal, cg["cell_volume"].sum(), ds.domain_width.prod()
            for g in ds.index.grids:
                di = g.get_global_startindex()
                dd = g.ActiveDimensions
                for i in range(dn):
                    f = cg["density"][dn*di[0]+i:dn*(di[0]+dd[0])+i:dn,
                                      dn*di[1]+i:dn*(di[1]+dd[1])+i:dn,
                                      dn*di[2]+i:dn*(di[2]+dd[2])+i:dn]
                    yield assert_equal, f, g["density"]

def test_smoothed_covering_grid():
    # We decompose in different ways
    for level in [0, 1, 2]:
        for nprocs in [1, 2, 4, 8]:
            ds = fake_random_ds(16, nprocs = nprocs)
            dn = ds.refine_by**level 
            cg = ds.smoothed_covering_grid(level, [0.0, 0.0, 0.0],
                    dn * ds.domain_dimensions)
            yield assert_equal, cg["ones"].max(), 1.0
            yield assert_equal, cg["ones"].min(), 1.0
            yield assert_equal, cg["cell_volume"].sum(), ds.domain_width.prod()
            for g in ds.index.grids:
                if level != g.Level: continue
                di = g.get_global_startindex()
                dd = g.ActiveDimensions
                for i in range(dn):
                    f = cg["density"][dn*di[0]+i:dn*(di[0]+dd[0])+i:dn,
                                      dn*di[1]+i:dn*(di[1]+dd[1])+i:dn,
                                      dn*di[2]+i:dn*(di[2]+dd[2])+i:dn]
                    yield assert_equal, f, g["density"]


def test_arbitrary_grid():
    for ncells in [64, 128, 256]:
        for px in [0.125, 0.25, 0.55519]:

            particle_data = {
                'particle_position_x': np.array([px]),
                'particle_position_y': np.array([0.5]),
                'particle_position_z': np.array([0.5]),
                'particle_mass': np.array([1.0])}

            ds = load_particles(particle_data)

            LE = np.array([0.05, 0.05, 0.05])
            RE = np.array([0.95, 0.95, 0.95])
            dims = np.array([ncells, ncells, ncells])

            dds = (RE - LE) / dims
            volume = ds.quan(np.product(dds), 'cm**3')

            obj = ds.arbitrary_grid(LE, RE, dims)
            deposited_mass = obj["deposit", "all_density"].sum() * volume

            yield assert_equal, deposited_mass, ds.quan(1.0, 'g')

    # Test that we get identical results to the covering grid for unigrid data.
    # Testing AMR data is much harder.
    for nprocs in [1, 2, 4, 8]:
        ds = fake_random_ds(32, nprocs = nprocs)
        for ref_level in [0, 1, 2]:
            cg = ds.covering_grid(ref_level, [0.0, 0.0, 0.0],
                    2**ref_level * ds.domain_dimensions)
            ag = ds.arbitrary_grid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0],
                    2**ref_level * ds.domain_dimensions)
            yield assert_almost_equal, cg["density"], ag["density"]

output_00080 = "output_00080/info_00080.txt"
@requires_file(output_00080)
def test_octree_cg():
    ds = load(output_00080)
    cgrid = ds.covering_grid(0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
    density_field = cgrid["density"]
    assert_equal((density_field == 0.0).sum(), 0)

ekh = 'EnzoKelvinHelmholtz/DD0011/DD0011'
@requires_file(ekh)
def test_smoothed_covering_grid_2d_dataset():
    ds = load(ekh)
    ds.periodicity = (True, True, True)
    scg = ds.smoothed_covering_grid(1, [0, 0, 0], [128, 128, 1])
    assert_equal(scg['density'].shape, [128, 128, 1])
