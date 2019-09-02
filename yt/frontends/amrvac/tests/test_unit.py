import numpy as np

import yt
from yt.utilities.parameter_file_storage import output_type_registry
#from yt.extern.six.moves import builtins

# todo : add regression tests (see enzo)
# todo : remove tmpdir
tmpdir = "data/"
sample_datasets = dict(
    blastwave_cartesian_3D = tmpdir + "bw_3d0000.dat",
    blastwave_polar_2D = tmpdir + "bw_polar_2D0000.dat",
    #blastwave_cylindrical_2D = tmpdir + "",
    #blastwave_cylindrical_3D = tmpdir + "",
    #kelvin_helmoltz_cartesian_2D = tmpdir + ""
)
correct_geometries = dict(
    blastwave_cartesian_3D = "cartesian",
    blastwave_polar_2D = "polar",
    #blastwave_cylindrical_2D = "cylindrical",
    #blastwave_cylindrical_3D = "cylindrical",
    #kelvin_helmoltz_cartesian_2D = "cartesian"
)

class TestParsing:
    def test_frontends_contents(self):
        """Check that the frontend is correctly registered."""
        assert "amrvac" in dir(yt.frontends)

    def test_output_type_registry(self):
        """Check that the Dataset child class is part of classes tried for validation."""
        keys = list(output_type_registry.keys())
        assert "AMRVACDataset" in keys

    def test_is_valid(self):
        """Check that validation fonction works at least on example files."""
        for file in sample_datasets.values():
            assert output_type_registry["AMRVACDataset"]._is_valid(file)

    def test_load_amrvac(self):
        """Check that that data is marked amrvac-formated"""
        for file in sample_datasets.values():
            ds = yt.load(file)
            assert isinstance(ds, yt.frontends.amrvac.AMRVACDataset)

    def test_geometry_parsing(self):
        for simname in sample_datasets.keys():
            ds = yt.load(sample_datasets[simname])
            assert ds.geometry == correct_geometries[simname]

class AMRVAC_IO_Tests:
    def test_print_stats(self):
        """Check that print_stats() executes normally"""
        for file in sample_datasets:
            ds = yt.load(file)
            ds.print_stats()

    def test_fields(self):
        """Check for fields in .dat file"""
        for file in sample_datasets:
            ds = yt.load(file)
            field_labels = [f[1] for f in ds.field_list]
            if "m1" in field_labels:
                for n in range(2, ds.dimensionality):
                    self.assertTrue("m%d" % n in field_labels)
            if "b1" in field_labels:
                for n in range(2, ds.dimensionality):
                    self.assertTrue("b%d" % n in field_labels)

    # data reading
    def test_get_data(self):
        for file in sample_datasets.values():
            ds = yt.load(file)
            ad = ds.all_data()
            data = ad.get_data()

    def test_grid_dataread(self):
        for file in sample_datasets.values():
            ds = yt.load(file)
            grids = ds.index.grids
            # select random grid
            g = grids[0]
            # read in density
            rho = g["density"]
            self.assertTrue(isinstance(rho, yt.units.yt_array.YTArray))

    def test_dataread_all(self):
        for file in sample_datasets.values():
            ds = yt.load(file)
            ad = ds.all_data()
            self.assertTrue(isinstance(ad['rho'], yt.units.yt_array.YTArray))



class AMRVAC_Plot_Tests:
    def test_projection(self):
        """"Check that no error is raised"""
        for file in sample_datasets.values():
            ds = yt.load(file)
            if ds.dimensionality == 3:
                p = yt.ProjectionPlot(ds, 'x', 'density')

    def test_slice(self):
        """"Check that no error is raised"""
        ds = yt.load(sample_datasets["blastwave_cartesian_3D"])
        p = yt.SlicePlot(ds, 'x', 'density')


ds = yt.load(sample_datasets["blastwave_cartesian_3D"])
class unit_tests_3Dblast:
    def test_domain_size(self):
        """Check for correct box size, see bw_3d.par"""
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
            self.assertTrue(isinstance(g, yt.frontends.amrvac.AMRVACGrid))
            self.assertTrue(isinstance(g.LeftEdge, yt.units.yt_array.YTArray))
            self.assertTrue(isinstance(g.RightEdge, yt.units.yt_array.YTArray))
            self.assertTrue(isinstance(g.ActiveDimensions, np.ndarray))
            self.assertTrue(isinstance(g.Level, (np.int32, np.int64, int)))
        self.assertEqual(ds.index.max_level, 2)

    def test_grids(self):
        """Check grid parents, children, etc."""
        grids_maxlvl = ds.index.select_grids(ds.index.max_level)
        for g in grids_maxlvl:
            self.assertTrue(type(g) == yt.frontends.amrvac.AMRVACGrid)





if __name__ == '__main__':
    # Run this with blast_wave 3D (which is MHD)
    # if MHD runs in 3D, then 2D and HD will also work.
    unittest.main()
