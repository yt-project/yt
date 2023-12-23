import glob
import os

# import pytest
import numpy as np
from numpy.testing import assert_raises

from yt.config import ytcfg
from yt.data_objects.particle_filters import particle_filter
from yt.data_objects.time_series import DatasetSeries
from yt.testing import fake_particle_ds
from yt.utilities.answer_testing.framework import GenericArrayTest, requires_ds
from yt.utilities.exceptions import YTIllDefinedParticleData


def setup():
    ytcfg["yt", "internals", "within_testing"] = True


data_path = ytcfg.get("yt", "test_data_dir")

pfields = [
    ("all", "particle_position_x"),
    ("all", "particle_position_y"),
    ("all", "particle_position_z"),
]
vfields = [
    ("all", "particle_velocity_x"),
    ("all", "particle_velocity_y"),
    ("all", "particle_velocity_z"),
]

# For switching to pytest
# @pytest.fixture
# def particle_trajectories_test_dataset():
#    n_particles = 2
#    n_steps = 2
#    ids = np.arange(n_particles, dtype="int64")
#    data = {"particle_index": ids}
#    fields = [
#        "particle_position_x",
#        "particle_position_y",
#        "particle_position_z",
#        "particle_velocity_x", # adding a non-default field
#        "particle_index",
#    ]
#    negative = [False, False, False, True, False]
#    units = ["cm", "cm", "cm", "cm/s", "1"]

#    ts = DatasetSeries(
#        [
#            fake_particle_ds(
#                fields=fields,
#                negative=negative,
#                units=units,
#                npart=n_particles,
#                data=data,
#            )
#            for i in range(n_steps)
#        ]
#    )
#    return ts


@requires_ds("Orbit/orbit_hdf5_chk_0000")
def test_orbit_traj():
    fields = ["particle_velocity_x", "particle_velocity_y", "particle_velocity_z"]
    my_fns = glob.glob(os.path.join(data_path, "Orbit/orbit_hdf5_chk_00[0-9][0-9]"))
    my_fns.sort()
    ts = DatasetSeries(my_fns)
    ds = ts[0]
    traj = ts.particle_trajectories([1, 2], fields=fields, suppress_logging=True)
    for field in pfields + vfields:

        def field_func(name):
            return traj[field]  # noqa: B023

        yield GenericArrayTest(ds, field_func, args=[field])


@requires_ds("enzo_tiny_cosmology/DD0000/DD0000")
def test_etc_traj():
    fields = [
        ("all", "particle_velocity_x"),
        ("all", "particle_velocity_y"),
        ("all", "particle_velocity_z"),
    ]
    my_fns = glob.glob(
        os.path.join(data_path, "enzo_tiny_cosmology/DD000[0-9]/*.hierarchy")
    )
    my_fns.sort()
    ts = DatasetSeries(my_fns)
    ds = ts[0]
    sp = ds.sphere("max", (0.5, "Mpc"))
    indices = sp["particle_index"][sp["particle_type"] == 1][:5]
    traj = ts.particle_trajectories(indices, fields=fields, suppress_logging=True)
    traj.add_fields([("gas", "density")])
    for field in pfields + vfields + [("gas", "density")]:

        def field_func(name):
            return traj[field]  # noqa: B023

        yield GenericArrayTest(ds, field_func, args=[field])


def test_uniqueness():
    n_particles = 2
    n_steps = 2
    ids = np.arange(n_particles, dtype="int64") % (n_particles // 2)
    data = {"particle_index": ids}
    fields = [
        "particle_position_x",
        "particle_position_y",
        "particle_position_z",
        "particle_index",
    ]
    negative = [False, False, False, False]
    units = ["cm", "cm", "cm", "1"]

    ts = DatasetSeries(
        [
            fake_particle_ds(
                fields=fields,
                negative=negative,
                units=units,
                npart=n_particles,
                data=data,
            )
            for i in range(n_steps)
        ]
    )

    assert_raises(YTIllDefinedParticleData, ts.particle_trajectories, [0])


def test_ptype():
    n_particles = 100
    fields = [
        "particle_position_x",
        "particle_position_y",
        "particle_position_z",
        "particle_index",
        "particle_dummy",
    ]
    negative = [False, False, False, False, False]
    units = ["cm", "cm", "cm", "1", "1"]

    # Setup filters on the 'particle_dummy' field, keeping only the first 50
    @particle_filter(name="dummy", requires=["particle_dummy"])
    def dummy(pfilter, data):
        return data[(pfilter.filtered_type, "particle_dummy")] <= n_particles // 2

    # Setup fake particle datasets with repeated ids. This should work because
    # the ids are unique among `dummy_particles` so let's test this
    data = {
        "particle_index": np.arange(n_particles) % (n_particles // 2),
        "particle_dummy": np.arange(n_particles),
    }
    all_ds = [
        fake_particle_ds(
            fields=fields, negative=negative, units=units, npart=n_particles, data=data
        )
    ]
    for ds in all_ds:
        ds.add_particle_filter("dummy")
    ts = DatasetSeries(all_ds)

    # Select all dummy particles
    print(ts[0].derived_field_list)
    ids = ts[0].all_data()["dummy", "particle_index"]

    # Build trajectories
    ts.particle_trajectories(ids, ptype="dummy")


# nose version - remove when switching to pytest
def test_default_field_tuple():
    # setup DatasetSeries
    n_particles = 2
    n_steps = 2
    ids = np.arange(n_particles, dtype="int64")
    data = {"particle_index": ids}
    fields = [
        "particle_position_x",
        "particle_position_y",
        "particle_position_z",
        "particle_index",
    ]
    negative = [False, False, False, False]
    units = ["cm", "cm", "cm", "1"]

    ts = DatasetSeries(
        [
            fake_particle_ds(
                fields=fields,
                negative=negative,
                units=units,
                npart=n_particles,
                data=data,
            )
            for i in range(n_steps)
        ]
    )
    # actual test - not specifying ptype
    trajs = ts.particle_trajectories(ids, suppress_logging=True)
    for k in trajs.field_data.keys():
        assert isinstance(k, tuple), f"Expected key to be tuple, received {type(k)}"
        assert k in pfields, f"Unexpected field: {k}"

    # specifying ptype - should change the field type of the default fields
    trajs = ts.particle_trajectories(ids, ptype="io", suppress_logging=True)
    for k in trajs.field_data.keys():
        assert isinstance(k, tuple), f"Expected key to be tuple, received {type(k)}"
        assert k[0] == "io", f"Default field type ({k[0]}) does not match expected (io)"
        assert ("all", k[1]) in pfields, f"Unexpected field: {k[1]}"


# pytest version
# @pytest.mark.parametrize("ptype",
#    [
#        None,
#        ("io")
#    ]
# )
# def test_default_field_tuple(particle_trajectories_test_dataset,ptype):
#    ds = particle_trajectories_test_dataset[0]
#    ids = ds.all_data()[("all","particle_index")]
#    trajs = particle_trajectories_test_dataset.particle_trajectories(ids,
#                ptype=ptype, suppress_logging=True)
#    ptype = ptype if ptype else "all" # ptype defaults to "all"
#    for k in trajs.field_data.keys():
#        assert isinstance(k,tuple), f"Expected key to be tuple, received {type(k)}"
#        assert k[0]==ptype, f"Default field type ({k[0]}) does not match expected ({ptype})"
#        assert ("all",k[1]) in pfields, f"Unexpected field: {k[1]}"


# nose version - remove when switching to pytest
def test_time_and_index():
    # setup ParticleTrajectories
    n_particles = 2
    n_steps = 2
    ids = np.arange(n_particles, dtype="int64")
    data = {"particle_index": ids}
    fields = [
        "particle_position_x",
        "particle_position_y",
        "particle_position_z",
        "particle_index",
    ]
    negative = [False, False, False, False]
    units = ["cm", "cm", "cm", "1"]

    ts = DatasetSeries(
        [
            fake_particle_ds(
                fields=fields,
                negative=negative,
                units=units,
                npart=n_particles,
                data=data,
            )
            for i in range(n_steps)
        ]
    )
    # actual test - default ptype
    trajs = ts.particle_trajectories(ids, suppress_logging=True)
    traj = trajs.trajectory_from_index(1)
    for field in ["particle_time", "particle_index"]:
        assert ("all", field) in traj.keys(), f"Missing ('all',{field})"
        assert (field) not in traj.keys(), f"{field} present as bare string"

    # actual test - io ptype
    trajs = ts.particle_trajectories(ids, ptype="io", suppress_logging=True)
    traj = trajs.trajectory_from_index(1)
    for field in ["particle_time", "particle_index"]:
        assert ("io", field) in traj.keys(), f"Missing (io,{field})"
        assert (field) not in traj.keys(), f"{field} present as bare string"


# pytest version
# @pytest.mark.parametrize("ptype",
#    [
#        None,
#        ("io")
#    ]
# )
# def test_time_and_index(particle_trajectories_test_dataset,ptype):
#    ds = particle_trajectories_test_dataset[0]
#    ids = ds.all_data()[("all","particle_index")]
#    trajs = particle_trajectories_test_dataset.particle_trajectories(ids,
#                ptype=ptype, suppress_logging=True)
#    ptype = ptype if ptype else "all" # ptype defaults to "all"
#    traj = trajs.trajectory_from_index(1)
#    for field in ["particle_time","particle_index"]:
#        assert (ptype,field) in traj.keys(), f"Missing ({ptype},{field})"
#        assert (field) not in traj.keys(), f"{field} present as bare string"
