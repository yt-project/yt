import os
import tempfile

import numpy as np
import pytest

from yt.loaders import load, load_particles, load_uniform_grid
from yt.testing import assert_equal, assert_raises
from yt.utilities.exceptions import (
    YTInconsistentGridFieldShape,
    YTInconsistentGridFieldShapeGridDims,
    YTInconsistentParticleFieldShape,
    YTOutputNotIdentified,
)

# Globals
OCT_MASK_LIST = [
    8,
    0,
    0,
    0,
    0,
    8,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    8,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
]


@pytest.mark.answer_test
class TestStream:
    answer_file = None
    saved_hashes = None

    @pytest.mark.usefixtures("temp_dir")
    def test_load_empty_file(self):
        assert_raises(YTOutputNotIdentified, load, "not_a_file")
        assert_raises(
            YTOutputNotIdentified, load, tempfile.mkstemp("empty_file", dir=os.getcwd())
        )
        assert_raises(YTOutputNotIdentified, load, tempfile.mkdtemp(dir=os.getcwd()))

    def test_dimensionless_field_units(self):
        Z = np.random.uniform(size=(32, 32, 32))
        d = np.random.uniform(size=(32, 32, 32))
        data = {"density": d, "metallicity": Z}
        ds = load_uniform_grid(data, (32, 32, 32))
        dd = ds.all_data()
        assert_equal(Z.max(), dd["metallicity"].max())

    def test_inconsistent_field_shape(self):
        def load_field_field_mismatch():
            d = np.random.uniform(size=(32, 32, 32))
            t = np.random.uniform(size=(32, 64, 32))
            data = {"density": d, "temperature": t}
            load_uniform_grid(data, (32, 32, 32))

        def load_field_grid_mismatch():
            d = np.random.uniform(size=(32, 32, 32))
            t = np.random.uniform(size=(32, 32, 32))
            data = {"density": d, "temperature": t}
            load_uniform_grid(data, (32, 64, 32))

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

        assert_raises(YTInconsistentGridFieldShape, load_field_field_mismatch)
        assert_raises(YTInconsistentGridFieldShapeGridDims, load_field_grid_mismatch)
        assert_raises(YTInconsistentParticleFieldShape, load_particle_fields_mismatch)
