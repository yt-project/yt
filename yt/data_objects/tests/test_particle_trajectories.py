import numpy as np
from numpy.testing import assert_raises

from yt.config import ytcfg
from yt.data_objects.time_series import DatasetSeries
from yt.testing import fake_particle_ds, assert_allclose_units
from yt.utilities.exceptions import YTIllDefinedParticleData
from yt.data_objects.particle_filters import particle_filter


def setup():
    ytcfg["yt","__withintesting"] = "True"

data_path = ytcfg.get("yt", "test_data_dir")

pfields = ["particle_position_x", "particle_position_y", "particle_position_z"]

def test_orbit_traj():
    n_particles = 16 ** 2
    n_ds = 5
    ids = np.arange(n_particles, dtype=int)
    pos = np.random.rand(n_particles) / 2.
    data =[
        {'particle_index': ids.copy(),
         'particle_position_x': pos.copy() * 2.,
         'particle_position_y': pos.copy() * 1.,
         'particle_position_z': pos.copy() * 0.5}
        for i in range(n_ds)
        ]

    fields = ['particle_position_x', 'particle_position_y', 'particle_position_z',
              'particle_index', 'particle_mass']
    field_units = ['cm', 'cm', 'cm', '1', 'g']

    ts = DatasetSeries([fake_particle_ds(fields=fields, units=field_units,
                        data=dt, negative=[False]*len(fields),
                        npart=n_particles)
                    for dt in data])

    traj = ts.particle_trajectories(ids)

    for field in pfields:
        d_ax = np.diff(traj[field], axis=1)
        u_ax = ts.outputs[0].quan(0.0, 'code_length')
        assert_allclose_units(d_ax, u_ax)

def test_etc_traj():
    n_particles = 16 ** 2
    n_ds = 5
    ids = np.arange(n_particles, dtype=int)
    pos = np.random.rand(n_particles) / 2.
    data =[
        {'particle_index': ids.copy(),
         'particle_position_x': pos.copy() * 2.,
         'particle_position_y': pos.copy() * 1.,
         'particle_position_z': pos.copy() * 0.5}
        for i in range(n_ds)
        ]

    fields = ['particle_position_x', 'particle_position_y', 'particle_position_z',
              'particle_index', 'particle_mass', 'particle_type']
    field_units = ['cm', 'cm', 'cm', '1', 'g', 'dimensionless']

    ts = DatasetSeries([fake_particle_ds(fields=fields, units=field_units,
                        data=dt, negative=[False]*len(fields),
                        npart=n_particles)
                    for dt in data])

    ds = ts[0]
    sp = ds.sphere(ds.domain_center / 2.0, (0.5, "Mpc"))
    indices = sp["particle_index"][sp["particle_type"] <= 0.5][:5]
    traj = ts.particle_trajectories(indices)

    for field in pfields:
        d_ax = np.diff(traj[field], axis=1)
        u_ax = ts.outputs[0].quan(0.0, 'code_length')
        assert_allclose_units(d_ax, u_ax)

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

def test_ptype():
    n_particles = 100
    fields = ['particle_position_x', 'particle_position_y', 'particle_position_z', 'particle_index',
              'particle_dummy']
    negative = [False, False, False, False, False]
    units = ['cm', 'cm', 'cm', '1', '1']

    # Setup filters on the 'particle_dummy' field, keeping only the first 50
    @particle_filter(name='dummy', requires=["particle_dummy"])
    def dummy(pfilter, data):
        return data[(pfilter.filtered_type, "particle_dummy")] <= n_particles // 2

    # Setup fake particle datasets with repeated ids. This should work because
    # the ids are unique among `dummy_particles` so let's test this
    data = {'particle_index': np.arange(n_particles) % (n_particles // 2),
            'particle_dummy': np.arange(n_particles)}
    all_ds = [fake_particle_ds(fields=fields, negative=negative, units=units,
                               npart=n_particles, data=data)]
    for ds in all_ds:
        ds.add_particle_filter('dummy')
    ts = DatasetSeries(all_ds)

    # Select all dummy particles
    print(ts[0].derived_field_list)
    ids = ts[0].all_data()['dummy', 'particle_index']

    # Build trajectories
    ts.particle_trajectories(ids, ptype='dummy')
