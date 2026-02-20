from yt.testing import fake_amr_ds


def test_fake_amr_ds_particles_within_grid_bounds():
    ds = fake_amr_ds(particles=100)

    for g in ds.index.grids:
        px = g["io", "particle_position_x"]
        le = g.LeftEdge[0].value
        re = g.RightEdge[0].value
        assert (px >= le).all() and (px < re).all()
