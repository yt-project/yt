import numpy as np

from yt.testing import assert_allclose
from yt.utilities.lib.particle_mesh_operations import CICSample_3


def setup():
    pass


def test_sample():

    grid = {}

    dims = np.array([64, 64, 64], dtype="int32")

    inds = np.indices(dims)
    grid["x"] = inds[0] + 0.5
    grid["y"] = inds[1] + 0.5
    grid["z"] = inds[2] + 0.5

    num_particles = np.int64(1000)

    xp = np.random.uniform(low=1.0, high=63.0, size=num_particles)
    yp = np.random.uniform(low=1.0, high=63.0, size=num_particles)
    zp = np.random.uniform(low=1.0, high=63.0, size=num_particles)

    xfield = np.zeros(num_particles)
    yfield = np.zeros(num_particles)
    zfield = np.zeros(num_particles)

    dx = 1.0
    le = np.zeros(3)

    CICSample_3(xp, yp, zp, xfield, num_particles, grid["x"], le, dims, dx)
    CICSample_3(xp, yp, zp, yfield, num_particles, grid["y"], le, dims, dx)
    CICSample_3(xp, yp, zp, zfield, num_particles, grid["z"], le, dims, dx)

    assert_allclose(xp, xfield)
    assert_allclose(yp, yfield)
    assert_allclose(zp, zfield)
