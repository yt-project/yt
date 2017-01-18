import glob

import os

from yt.data_objects.time_series import DatasetSeries
from yt.config import ytcfg

def setup():
    ytcfg["yt","__withintesting"] = "True"

data_path = ytcfg.get("yt", "test_data_dir")

def test_orbit_traj():
    fields = ["particle_velocity_x", "particle_velocity_y", "particle_velocity_z"]
    my_fns = glob.glob(os.path.join(data_path, "Orbit/orbit_hdf5_chk_00[0-9][0-9]"))
    my_fns.sort()
    ts = DatasetSeries(my_fns)
    traj = ts.particle_trajectories([1, 2], fields=fields, suppress_logging=True)

def test_etc_traj():
    fields = ["particle_velocity_x", "particle_velocity_y", "particle_velocity_z"]
    my_fns = glob.glob(os.path.join(data_path, "enzo_tiny_cosmology/DD*/*.hierarchy"))
    my_fns.sort()
    ts = DatasetSeries(my_fns)
    sp = ts[0].sphere("max", (0.5, "Mpc"))
    indices = sp["particle_index"][sp["particle_type"] == 1]
    traj = ts.particle_trajectories(indices, fields=fields, suppress_logging=True)
    traj.add_fields(["density"])