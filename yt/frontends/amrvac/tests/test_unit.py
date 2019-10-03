import numpy as np

import yt
from yt.utilities.parameter_file_storage import output_type_registry

# todo : remove tmpdir
tmpdir = "data/"
sample_datasets = dict(
    blastwave_cartesian_3D = tmpdir + "bw_3d0000.dat",
    blastwave_polar_2D = tmpdir + "bw_polar_2D0000.dat",
    blastwave_spherical_2D = tmpdir + "bw_2d0000.dat",
    blastwave_cylindrical_3D = tmpdir + "bw_cylindrical_3D0000.dat",
    khi_cartesian_2D = tmpdir + "kh_2D0000.dat",
    khi_cartesian_3D = tmpdir + "kh_3D0000.dat",
    jet_cylindrical_25D = tmpdir + "Jet0003.dat",
    riemann_cartesian_175D = tmpdir + "R_1d0005.dat"
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
        for simname, file in sample_datasets.items():
            ds = yt.load(file)
            correct_geom = simname.split("_")[1]
            assert ds.geometry == correct_geom

class TestIO:
    def test_print_stats(self):
        """Check that print_stats() executes normally"""
        for file in sample_datasets.values():
            ds = yt.load(file)
            ds.print_stats()

    def test_fields(self):
        """Check for fields in .dat file"""
        for simname, file in sample_datasets.items():
            ds = yt.load(file)
            field_labels = [f[1] for f in ds.field_list]
            dims = int(simname.split("_")[2][:-1])
            ndim = ds.dimensionality
            ndir = ndim
            if dims > 100:
                ndir += 2
            elif dims > 10:
                ndir += 1
            if "m1" in field_labels:
                for n in range(2, ndir+1):
                    assert "m%d" % n in field_labels
            if "b1" in field_labels:
                for n in range(2, ndir+1):
                    assert "b%d" % n in field_labels

    def test_get_data(self):
        for file in sample_datasets.values():
            ds = yt.load(file)
            ad = ds.all_data()
            data = ad.get_data() # noqa: F841

    def test_grid_dataread(self):
        for file in sample_datasets.values():
            ds = yt.load(file)
            grids = ds.index.grids
            # select random grid
            g = grids[0]
            # read in density
            rho = g["density"]
            assert isinstance(rho, yt.units.yt_array.YTArray)

    def test_dataread_all(self):
        for file in sample_datasets.values():
            ds = yt.load(file)
            ad = ds.all_data()
            assert isinstance(ad['rho'], yt.units.yt_array.YTArray)



class TestPlot:
    def test_projection(self):
        """"Check that no error is raised"""
        for file in sample_datasets.values():
            ds = yt.load(file)
            if ds.dimensionality == 3:
                axis = {"cartesian": "x", "polar": "r", "cylindrical": "r", "spherical": "r"}[ds.geometry]
                p = yt.ProjectionPlot(ds, axis, 'density') # noqa: F841

    def test_slice(self):
        """"Check that no error is raised"""
        for file in sample_datasets.values():
            ds = yt.load(file)
            if ds.dimensionality == 3:
                axis = {"cartesian": "x", "polar": "r", "cylindrical": "r", "spherical": "r"}[ds.geometry]
                p = yt.SlicePlot(ds, axis, 'density') # noqa: F841


ds = yt.load(sample_datasets["blastwave_cartesian_3D"])
class TestOutput:
    def test_domain_size(self):
        """Check for correct box size, see bw_3d.par"""
        for lb in ds.domain_left_edge:
            assert int(lb) == 0
        for rb in ds.domain_right_edge:
            assert int(rb) == 2
        for w in ds.domain_width:
            assert int(w) == 2

    def test_grid_attributes(self):
        """Check various grid attributes"""
        grids = ds.index.grids
        assert ds.index.max_level == 2
        for g in grids:
            assert isinstance(g, yt.frontends.amrvac.AMRVACGrid)
            assert isinstance(g.LeftEdge, yt.units.yt_array.YTArray)
            assert isinstance(g.RightEdge, yt.units.yt_array.YTArray)
            assert isinstance(g.ActiveDimensions, np.ndarray)
            assert isinstance(g.Level, (np.int32, np.int64, int))

