from collections import OrderedDict
import glob
import os

import numpy as np
from numpy.testing import \
    assert_raises
import pytest

from yt.config import ytcfg
from yt.data_objects.particle_filters import particle_filter
from yt.data_objects.time_series import DatasetSeries
from yt.testing import fake_particle_ds
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils
from yt.utilities.exceptions import YTIllDefinedParticleData


data_path = ytcfg.get("yt", "test_data_dir")


pfields = ["particle_position_x", "particle_position_y", "particle_position_z"]
vfields = ["particle_velocity_x", "particle_velocity_y", "particle_velocity_z"]


@pytest.mark.answer_test
@pytest.mark.usefixtures('answer_file')
class TestParticleTrajectories(fw.AnswerTest):
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds("Orbit/orbit_hdf5_chk_0000")
    def test_orbit_traj(self, field, orbit_traj):
        ds, traj = orbit_traj
        def field_func(name):
            return traj[field]
        ga_hd = self.generic_arraytest(ds, field_func, args=[field])
        self.hashes.update({'generic_array' : ga_hd})

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds("enzo_tiny_cosmology/DD0000/DD0000")
    def test_etc_traj(self, field, etc_traj):
        ds, traj = etc_traj
        def field_func(name):
            return traj[field]
        ga_hd = self.generic_array_test(ds, field_func, args=[field])
        self.hashes.update({'generic_array' : ga_hd})

    def test_uniqueness(self):
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

    def test_ptype(self):
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
