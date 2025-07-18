from yt.loaders import load
from yt.sample_data.api import _get_test_data_dir_path
from yt.testing import assert_array_equal, requires_file


@requires_file("TNGHalo/halo_59.hdf5")
def test_ewah_write_load(tmp_path):
    mock_file = tmp_path / "halo_59.hdf5"
    mock_file.symlink_to(_get_test_data_dir_path() / "TNGHalo" / "halo_59.hdf5")

    masses = []
    opts = [
        (None, True, (6, 2, 4), False),
        (None, True, (6, 2, 4), True),
        ((5, 3), False, (5, 3, 3), False),
        ((5, 3), False, (5, 3, 3), True),
    ]

    for opt in opts:
        ds = load(mock_file, index_order=opt[0])
        _, c = ds.find_max(("gas", "density"))
        assert ds.index.pii.mutable_index is opt[1]
        assert ds.index.pii.order1 == opt[2][0]
        assert ds.index.pii.order2_orig == opt[2][1]
        assert ds.index.pii.order2 == opt[2][2]
        assert ds.index.pii._is_loaded is opt[3]
        c += ds.quan(0.5, "Mpc")
        sp = ds.sphere(c, (20.0, "kpc"))
        mass = sp.sum(("gas", "mass"))
        masses.append(mass)

    assert_array_equal(masses[0], masses[1:])
