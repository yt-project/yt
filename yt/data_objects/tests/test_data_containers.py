import os
import shelve
import shutil
import tempfile
import unittest

import numpy as np
from nose.tools import assert_raises
from numpy.testing import assert_array_equal

from yt.data_objects.data_containers import YTDataContainer
from yt.data_objects.particle_filters import particle_filter
from yt.testing import assert_equal, fake_random_ds, fake_amr_ds,\
    fake_particle_ds, requires_module
from yt.utilities.exceptions import YTFieldNotFound, YTException

class TestDataContainers(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()
        cls.curdir = os.getcwd()
        os.chdir(cls.tmpdir)

    @classmethod
    def tearDownClass(cls):
        os.chdir(cls.curdir)
        shutil.rmtree(cls.tmpdir)

    def test_yt_data_container(self):
        # Test if ds could be None
        with assert_raises(RuntimeError) as err:
            YTDataContainer(None, None)
        desired = ('Error: ds must be set either through class'
                   ' type or parameter to the constructor')
        assert_equal(str(err.exception), desired)

        # Test if field_data key exists
        ds = fake_random_ds(5)
        proj = ds.proj("density", 0, data_source=ds.all_data())
        assert_equal('px' in proj.keys(), True)
        assert_equal('pz' in proj.keys(), False)

        # Delete the key and check if exits
        del proj['px']
        assert_equal('px' in proj.keys(), False)
        del proj['density']
        assert_equal('density' in proj.keys(), False)

        # Delete a non-existent field
        with assert_raises(YTFieldNotFound) as ex:
            del proj['p_mass']
        desired = "Could not find field '('stream', 'p_mass')' in UniformGridData."
        assert_equal(str(ex.exception), desired)

    def test_write_out(self):
        filename = "sphere.txt"
        ds = fake_particle_ds()
        sp = ds.sphere(ds.domain_center, 0.25)
        sp.write_out(filename, fields=["cell_volume"])

        with open(filename, "r") as file:
            file_row_1 = file.readline()
            file_row_2 = file.readline()
            file_row_2 = np.array(file_row_2.split('\t'), dtype=np.float64)
        sorted_keys = sorted(sp.field_data.keys())
        keys = [str(k) for k in sorted_keys]
        keys = "\t".join(["#"] + keys + ["\n"])
        data = [sp.field_data[k][0] for k in sorted_keys]

        assert_equal(keys, file_row_1)
        assert_array_equal(data, file_row_2)

        # Test for exception
        with assert_raises(YTException) as ex:
            sp.write_out(filename, fields=["particle_position_x"])
        desired = ("Field type ['all'] of the supplied field ['particle_position_x']"
                   " is in consistent with field type 'gas'.")
        assert_equal(str(ex.exception)[:50], desired[:50])

    def test_save_object(self):
        ds = fake_particle_ds()
        sp = ds.sphere(ds.domain_center, 0.25)
        sp.save_object("my_sphere_1", filename="test_save_obj")
        obj = shelve.open("test_save_obj", protocol=-1)
        loaded_sphere = obj["my_sphere_1"][1]
        obj.close()
        assert_array_equal(loaded_sphere.center, sp.center)
        assert_equal(loaded_sphere.radius, sp.radius)
        for k in loaded_sphere._key_fields:
            assert_array_equal(loaded_sphere[k], sp[k])

        # Object is saved but retrieval is not working
        # sp.save_object("my_sphere_2")
        # loaded_sphere = ds.index.load_object("my_sphere_2")
        # for k in loaded_sphere._key_fields:
        #     assert_array_equal(loaded_sphere[k], sp[k])

    @requires_module("pandas")
    def test_to_dataframe(self):
        fields = ["density", "velocity_z"]
        ds = fake_random_ds(6)
        dd = ds.all_data()
        df1 = dd.to_dataframe(fields)
        assert_array_equal(dd[fields[0]], df1[fields[0]])
        assert_array_equal(dd[fields[1]], df1[fields[1]])

    def test_std(self):
        ds = fake_random_ds(3)
        ds.all_data().std('density', weight="velocity_z")

    def test_to_frb(self):
        # Test cylindrical geometry
        fields = ["density", "cell_mass"]
        ds = fake_amr_ds(fields=fields, geometry="cylindrical", particles=16**3)
        dd = ds.all_data()
        proj = ds.proj("density", weight_field="cell_mass", axis=1, data_source=dd)
        frb = proj.to_frb((1.0, 'unitary'), 64)
        assert_equal(frb.radius, (1.0, 'unitary'))
        assert_equal(frb.buff_size, 64)

    def test_extract_isocontours(self):
        # Test isocontour properties for AMRGridData
        ds = fake_amr_ds(fields=["density", "cell_mass"], particles=16**3)
        dd = ds.all_data()
        q = dd.quantities["WeightedAverageQuantity"]
        rho = q("density", weight="cell_mass")
        dd.extract_isocontours("density", rho, "triangles.obj", True)
        dd.calculate_isocontour_flux("density", rho, "x", "y", "z", "dx")

        # Test error in case of ParticleData
        ds = fake_particle_ds()
        dd = ds.all_data()
        q = dd.quantities["WeightedAverageQuantity"]
        rho = q("particle_velocity_x", weight="particle_mass")
        with assert_raises(NotImplementedError):
            dd.extract_isocontours("density", rho, sample_values='x')

    def test_derived_field(self):
        # Test that derived field on filtered particles do not require
        # their parent field to be created
        ds = fake_particle_ds()
        dd = ds.all_data()

        @particle_filter(requires=['particle_mass'], filtered_type='io')
        def massive(pfilter, data):
            return data[(pfilter.filtered_type, 'particle_mass')].to('code_mass') > 0.5

        ds.add_particle_filter('massive')

        def fun(field, data):
            return data[field.name[0], 'particle_mass']

        # Add the field to the massive particles
        ds.add_field(('massive', 'test'), function=fun,
                     sampling_type='particle', units='code_mass')

        expected_size = (dd['io', 'particle_mass'].to('code_mass') > 0.5).sum()

        fields_to_test = (f for f in ds.derived_field_list
                          if f[0] == 'massive')

        def test_this(fname):
            data = dd[fname]
            assert_equal(data.shape[0], expected_size)

        for fname in fields_to_test:
            test_this(fname)
