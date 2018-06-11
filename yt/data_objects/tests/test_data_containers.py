import shelve

import numpy as np
from nose.tools import assert_raises
from numpy.testing import assert_array_equal

from yt.data_objects.data_containers import YTDataContainer
from yt.testing import assert_equal, fake_random_ds, fake_amr_ds, fake_particle_ds
from yt.utilities.exceptions import YTFieldNotFound


def test_yt_data_container():
    # Test if ds could be None
    with assert_raises(RuntimeError) as err:
        YTDataContainer(None, None)
    desired = 'Error: ds must be set either through class type or parameter' \
              ' to the constructor'
    assert_equal(str(err.exception), desired)

    # Test if field_data key exists
    ds = fake_random_ds(5)
    proj = ds.proj("density", 0, data_source=ds.all_data())
    assert_equal('px' in proj.keys(), True)
    assert_equal('pz' in proj.keys(), False)

    # Delete the key and check if exits
    proj.__delitem__('px')
    assert_equal('px' in proj.keys(), False)
    proj.__delitem__('density')
    assert_equal('density' in proj.keys(), False)

    # Delete a non-existent field
    with assert_raises(YTFieldNotFound) as ex:
        proj.__delitem__('p_mass')
    desired = "Could not find field '('stream', 'p_mass')' in UniformGridData."
    assert_equal(str(ex.exception), desired)

    # Test the convert method
    assert_equal(proj.convert('HydroMethod'), -1)

def test_write_out():
    filename = "sphere.txt"
    ds = fake_particle_ds()
    sp = ds.sphere(ds.domain_center, 0.25)
    sp.write_out(filename)

    with open(filename, "r") as file:
        file_row_1 = file.readline()
        file_row_2 = file.readline()
        file_row_2 = np.array(file_row_2.split('\t'), dtype=np.float64)
    sorted_keys = sorted(sp.field_data.keys())
    _keys = [str(k) for k in sorted_keys]
    _keys = "\t".join(["#"] + _keys + ["\n"])
    _data = [sp.field_data[k][0] for k in sorted_keys]

    assert_equal(_keys, file_row_1)
    assert_array_equal(_data, file_row_2)

def test_save_object():
    ds = fake_particle_ds()
    sp = ds.sphere(ds.domain_center, 0.25)
    sp.save_object("my_sphere_1", filename="test_save_obj")
    with shelve.open("test_save_obj", protocol=-1) as obj:
        loaded_sphere = obj["my_sphere_1"][1]
    assert_array_equal(loaded_sphere.center, sp.center)
    assert_equal(loaded_sphere.radius, sp.radius)

    sp.save_object("my_sphere_2")

def test_to_dataframe():
    try:
        import pandas as pd
        fields = ["density", "velocity_z"]
        ds = fake_random_ds(6)
        dd = ds.all_data()
        df1 = dd.to_dataframe(fields)
        assert_array_equal(dd[fields[0]], df1[fields[0]])
        assert_array_equal(dd[fields[1]], df1[fields[1]])
    except ImportError:
        pass

def test_std():
    ds = fake_random_ds(3)
    ds.all_data().std('density', weight="velocity_z")

def test_to_frb():
    ds = fake_amr_ds(fields=["density", "cell_mass"], geometry="cylindrical",
                     particles=16**3)
    proj = ds.proj("density", weight_field="cell_mass", axis=1,
                   data_source=ds.all_data())
    proj.to_frb((1.0, 'unitary'), 64)

def test_extract_isocontours():
    # Test isocontour properties for AMRGridData
    ds = fake_amr_ds(fields=["density", "cell_mass"], particles=16**3)
    dd = ds.all_data()
    rho = dd.quantities["WeightedAverageQuantity"]("density",
                                                   weight="cell_mass")
    dd.extract_isocontours("density", rho, "triangles.obj", True)
    dd.calculate_isocontour_flux("density", rho, "x", "y", "z",
                                 "dx")

    # Test error in case of ParticleData
    ds = fake_particle_ds()
    dd = ds.all_data()
    rho = dd.quantities["WeightedAverageQuantity"]("particle_velocity_x",
                                                   weight="particle_mass")
    with assert_raises(NotImplementedError):
        dd.extract_isocontours("density", rho, sample_values='x')

def test_extract_connected_sets():

    ds = fake_random_ds(16, nprocs=8, particles=16 ** 3)
    data_source = ds.disk([0.5, 0.5, 0.5], [0., 0., 1.], (8, 'kpc'), (1, 'kpc'))
    field = ("gas", "density")
    min_val, max_val = data_source[field].min() / 2, data_source[field].max() / 2

    data_source.extract_connected_sets(field, 3, min_val, max_val,
                                              log_space=True, cumulative=True)
