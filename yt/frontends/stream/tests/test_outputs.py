import os
import shutil
import tempfile
import unittest

import numpy as np

from yt.loaders import load_particles, load_uniform_grid
from yt.testing import assert_equal, assert_raises
from yt.utilities.exceptions import (
    YTInconsistentGridFieldShape,
    YTInconsistentGridFieldShapeGridDims,
    YTInconsistentParticleFieldShape,
)


class TestEmptyLoad(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.curdir = os.getcwd()
        os.chdir(self.tmpdir)

        # create 0 byte file
        open("empty_file", "a")

        # create empty directory
        os.makedirs("empty_directory")

    def tearDown(self):
        os.chdir(self.curdir)
        shutil.rmtree(self.tmpdir)


def test_dimensionless_field_units():
    Z = np.random.uniform(size=(32, 32, 32))
    d = np.random.uniform(size=(32, 32, 32))

    data = {"density": d, "metallicity": Z}

    ds = load_uniform_grid(data, (32, 32, 32))

    dd = ds.all_data()

    assert_equal(Z.max(), float(dd["stream", "metallicity"].max()))


def test_inconsistent_field_shape():
    def load_field_field_mismatch():
        d = np.random.uniform(size=(32, 32, 32))
        t = np.random.uniform(size=(32, 64, 32))
        data = {"density": d, "temperature": t}
        load_uniform_grid(data, (32, 32, 32))

    assert_raises(YTInconsistentGridFieldShape, load_field_field_mismatch)

    def load_field_grid_mismatch():
        d = np.random.uniform(size=(32, 32, 32))
        t = np.random.uniform(size=(32, 32, 32))
        data = {"density": d, "temperature": t}
        load_uniform_grid(data, (32, 64, 32))

    assert_raises(YTInconsistentGridFieldShapeGridDims, load_field_grid_mismatch)

    def load_particle_fields_mismatch():
        x = np.random.uniform(size=100)
        y = np.random.uniform(size=100)
        z = np.random.uniform(size=200)
        data = {
            "particle_position_x": x,
            "particle_position_y": y,
            "particle_position_z": z,
        }
        load_particles(data)

    assert_raises(YTInconsistentParticleFieldShape, load_particle_fields_mismatch)


def test_parameters():
    # simple test to check that we can pass in parameters
    Z = np.random.uniform(size=(32, 32, 32))
    d = np.random.uniform(size=(32, 32, 32))

    data = {"density": d, "metallicity": Z}

    ds = load_uniform_grid(data, (32, 32, 32), parameters={"metadata_is_nice": True})
    assert ds.parameters["metadata_is_nice"]
