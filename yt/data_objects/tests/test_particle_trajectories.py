import glob
import os
import numpy as np
from numpy.testing import \
    assert_raises

from yt.config import ytcfg
from yt.data_objects.time_series import DatasetSeries
from yt.utilities.answer_testing.framework import \
    requires_ds, \
    GenericArrayTest
from yt.testing import fake_particle_ds
from yt.utilities.exceptions import YTIllDefinedParticleData


def setup():
    ytcfg["yt","__withintesting"] = "True"

data_path = ytcfg.get("yt", "test_data_dir")

pfields = ["particle_position_x", "particle_position_y", "particle_position_z"]
vfields = ["particle_velocity_x", "particle_velocity_y", "particle_velocity_z"]

@requires_ds("Orbit/orbit_hdf5_chk_0000")
def test_orbit_traj():
    fields = ["particle_velocity_x", "particle_velocity_y", "particle_velocity_z"]
    my_fns = glob.glob(os.path.join(data_path, "Orbit/orbit_hdf5_chk_00[0-9][0-9]"))
    my_fns.sort()
    ts = DatasetSeries(my_fns)
    ds = ts[0]
    traj = ts.particle_trajectories([1, 2], fields=fields, suppress_logging=True)
    for field in pfields+vfields:
        def field_func(name):
            return traj[field]
        yield GenericArrayTest(ds, field_func, args=[field])

@requires_ds("enzo_tiny_cosmology/DD0000/DD0000")
def test_etc_traj():
    fields = ["particle_velocity_x", "particle_velocity_y", "particle_velocity_z"]
    my_fns = glob.glob(os.path.join(data_path, "enzo_tiny_cosmology/DD000[0-9]/*.hierarchy"))
    my_fns.sort()
    ts = DatasetSeries(my_fns)
    ds = ts[0]
    sp = ds.sphere("max", (0.5, "Mpc"))
    indices = sp["particle_index"][sp["particle_type"] == 1][:5]
    traj = ts.particle_trajectories(indices, fields=fields, suppress_logging=True)
    traj.add_fields(["density"])
    for field in pfields+vfields+["density"]:
        def field_func(name):
            return traj[field]
        yield GenericArrayTest(ds, field_func, args=[field])

def test_uniqueness():
    n_particles = 2
    n_steps = 2
    ids = np.arange(n_particles, dtype=int) % (n_particles//2)
    data = {'particle_index': ids}
    fields = ['particle_position_x', 'particle_position_y', 'particle_position_z', 'particle_index']
    negative = [False, False, False, False]
    units = ['cm', 'cm', 'cm', '1']

    ts = DatasetSeries([fake_particle_ds(fields=fields, negative=negative, units=units,
                                         npart=n_particles, data=data)
                        for i in range(n_steps)])

    assert_raises(YTIllDefinedParticleData, ts.particle_trajectories, [0])
