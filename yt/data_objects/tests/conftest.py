import glob
import os

import pytest

from yt.config import ytcfg
from yt.data_objects.time_series import DatasetSeries

fields = ["particle_velocity_x", "particle_velocity_y", "particle_velocity_z"]
data_path = ytcfg.get("yt", "test_data_dir")
pfields = ["particle_position_x", "particle_position_y", "particle_position_z"]
vfields = ["particle_velocity_x", "particle_velocity_y", "particle_velocity_z"]


@pytest.fixture(scope="class")
def orbit_traj():
    my_fns = glob.glob(os.path.join(data_path, "Orbit/orbit_hdf5_chk_00[0-9][0-9]"))
    my_fns.sort()
    ts = DatasetSeries(my_fns)
    # If the data isn't present, trying to access ts[0] will raise an exception
    try:
        ds = ts[0]
        traj = ts.particle_trajectories([1, 2], fields=fields, suppress_logging=True)
        return [ds, traj]
    except IndexError:
        return pytest.skip("Data not found.")


@pytest.fixture(scope="class")
def etc_traj():
    my_fns = glob.glob(
        os.path.join(data_path, "enzo_tiny_cosmology/DD000[0-9]/*.hierarchy")
    )
    my_fns.sort()
    ts = DatasetSeries(my_fns)
    # If the data isn't present, trying to access ts[0] will raise an exception
    try:
        ds = ts[0]
        sp = ds.sphere("max", (0.5, "Mpc"))
        indices = sp["particle_index"][sp["particle_type"] == 1][:5]
        traj = ts.particle_trajectories(indices, fields=fields, suppress_logging=True)
        traj.add_fields(["density"])
        return [ds, traj]
    except IndexError:
        return pytest.skip("Data not found.")


def pytest_generate_tests(metafunc):
    if metafunc.function.__name__ == "test_orbit_traj":
        metafunc.parametrize(
            "field", pfields + vfields, ids=["x", "y", "z", "vx", "vy", "vz"]
        )
    if metafunc.function.__name__ == "test_etc_traj":
        metafunc.parametrize(
            "field",
            pfields + vfields + ["density"],
            ids=["x", "y", "z", "vx", "vy", "vz", "density"],
        )
