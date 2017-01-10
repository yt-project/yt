from yt.testing import fake_sph_orientation_ds

# The number of particles along each slice axis at that coordinate
SLICE_ANSWERS = {
    ('x', 0): 6,
    ('x', 0.5): 0,
    ('x', 1): 1,
    ('y', 0): 5,
    ('y', 1): 1,
    ('y', 2): 1,
    ('z', 0): 4,
    ('z', 1): 1,
    ('z', 2): 1,
    ('z', 3): 1,
}

def test_slice():
    ds = fake_sph_orientation_ds()
    for (sl_ax, coord), answer in SLICE_ANSWERS.items():
        sl = ds.slice(sl_ax, coord)
        assert sl['gas', 'density'].shape[0] == answer
