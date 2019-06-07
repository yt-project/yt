import sys
if sys.version.startswith("2"):
    import pathlib2 as pathlib
else:
    import pathlib
import yt
import numpy as np
import unittest

from yt.utilities.parameter_file_storage import output_type_registry

test_dir = pathlib.Path(__file__).parent.resolve()
example_data_2D = str(test_dir / "data/blast_wave_Cartesian_3D/bw_3d0001.dat")
ds = yt.load(example_data_2D)

class amrvac_unit_tests(unittest.TestCase):

    def test_frontends_contents(self):
        """Check that our frontend is correctly being registered."""
        self.assertTrue("amrvac" in dir(yt.frontends))

    def test_output_type_registry(self):
        """Check that the Dataset child class is part of classes tried for validation."""
        keys = list(output_type_registry.keys())
        self.assertTrue("AMRVACDataset" in keys)

    def test_is_amrvac_valid(self):
        """Check that our format validation fonction works at least on example file"""
        self.assertTrue(pathlib.Path(example_data_2D).is_file())
        self.assertTrue(output_type_registry["AMRVACDataset"]._is_valid(example_data_2D))

    def test_load_amrvac(self):
        """Check that that data is amrvac-formated"""
        self.assertTrue(ds.__class__ == yt.frontends.amrvac.AMRVACDataset)



    # ========== TESTS FOR DATA LOADING AND INSPECTION ==========
    def test_print_stats(self):
        """Check that print_stats() executes normally"""
        self.assertTrue(ds.print_stats() is None)

    def test_fields(self):
        """Check for fields in .dat file"""
        field_labels = [f[1] for f in ds.field_list]
        for f in ("rho", "m1", "m2", "m3", "e", "b1", "b2", "b3"):
            self.assertTrue(f in field_labels)

    def test_domain_size(self):
        """Check for correct box size, see bw_2d.par"""
        for lb in ds.domain_left_edge:
            self.assertEqual(int(lb), 0)
        for rb in ds.domain_right_edge:
            self.assertEqual(int(rb), 2)
        for w in ds.domain_width:
            self.assertEqual(int(w), 2)

    def test_grid_attributes(self):
        """Check various grid attributes"""
        grids = ds.index.grids
        for g in grids:
            self.assertTrue(type(g) == yt.frontends.amrvac.AMRVACGrid)
            self.assertTrue(type(g.LeftEdge) == yt.units.yt_array.YTArray)
            self.assertTrue(type(g.RightEdge) == yt.units.yt_array.YTArray)
            self.assertTrue(type(g.ActiveDimensions) == np.ndarray)
            self.assertTrue(type(g.Level) == np.int32)
        self.assertEqual(ds.index.max_level, 2)

    def test_grids(self):
        """Check grid parents, children, etc."""
        grids_maxlvl = ds.index.select_grids(ds.index.max_level)
        for g in grids_maxlvl:
            self.assertTrue(type(g) == yt.frontends.amrvac.AMRVACGrid)


    # ========== TEST FOR DATA READING ==========
    def test_grid_dataread(self):
        grids = ds.index.grids
        # select random grid
        g = grids[0]
        # read in density
        rho = g["density"]
        self.assertTrue(type(rho), yt.units.yt_array.YTArray)



if __name__ == '__main__':
    # Run this with blast_wave 3D (which is MHD)
    # if MHD runs in 3D, then 2D and HD will also work.
    unittest.main()