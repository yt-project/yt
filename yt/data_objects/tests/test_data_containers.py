import os

import numpy as np
import pytest
from numpy.testing import assert_array_equal, assert_equal

from yt.data_objects.data_containers import YTDataContainer
from yt.data_objects.particle_filters import particle_filter
from yt.testing import (
    fake_amr_ds,
    fake_particle_ds,
    fake_random_ds,
    requires_module_pytest as requires_module,
)
from yt.utilities.exceptions import YTException, YTFieldNotFound


@pytest.fixture
def temp_workdir(tmp_path):
    curdir = os.getcwd()
    os.chdir(tmp_path)
    yield tmp_path
    os.chdir(curdir)


@pytest.mark.usefixtures("temp_workdir")
def test_yt_data_container():
    # Test if ds could be None
    with pytest.raises(
        RuntimeError,
        match=(
            "Error: ds must be set either through class"
            " type or parameter to the constructor"
        ),
    ):
        YTDataContainer(None, None)

    # Test if field_data key exists
    ds = fake_random_ds(5)
    proj = ds.proj(("gas", "density"), 0, data_source=ds.all_data())
    assert "px" in proj.keys()
    assert "pz" not in proj.keys()

    # Delete the key and check if exits
    del proj["px"]
    assert "px" not in proj.keys()
    del proj["gas", "density"]
    assert "density" not in proj.keys()

    # Delete a non-existent field
    with pytest.raises(
        YTFieldNotFound, match="Could not find field 'p_mass' in UniformGridData."
    ):
        del proj["p_mass"]


@pytest.mark.usefixtures("temp_workdir")
def test_write_out():
    filename = "sphere.txt"
    ds = fake_random_ds(16, particles=10)
    sp = ds.sphere(ds.domain_center, 0.25)

    sp.write_out(filename, fields=[("gas", "cell_volume")])

    with open(filename) as file:
        file_row_1 = file.readline()
        file_row_2 = file.readline()
        file_row_2 = np.array(file_row_2.split("\t"), dtype=np.float64)
    sorted_keys = sorted(sp.field_data.keys())
    keys = [str(k) for k in sorted_keys]
    keys = "\t".join(["#"] + keys + ["\n"])
    data = [sp.field_data[k][0] for k in sorted_keys]

    assert_equal(keys, file_row_1)
    assert_array_equal(data, file_row_2)


@pytest.mark.usefixtures("temp_workdir")
def test_invalid_write_out():
    filename = "sphere.txt"
    ds = fake_random_ds(16, particles=10)
    sp = ds.sphere(ds.domain_center, 0.25)

    with pytest.raises(YTException):
        sp.write_out(filename, fields=[("all", "particle_ones")])


@pytest.mark.usefixtures("temp_workdir")
@requires_module("pandas")
def test_to_dataframe():
    fields = [("gas", "density"), ("gas", "velocity_z")]
    ds = fake_random_ds(6)
    dd = ds.all_data()
    df = dd.to_dataframe(fields)
    assert_array_equal(dd[fields[0]], df[fields[0][1]])
    assert_array_equal(dd[fields[1]], df[fields[1][1]])


@pytest.mark.usefixtures("temp_workdir")
@requires_module("astropy")
def test_to_astropy_table():
    from yt.units.yt_array import YTArray

    fields = [("gas", "density"), ("gas", "velocity_z")]
    ds = fake_random_ds(6)
    dd = ds.all_data()
    at1 = dd.to_astropy_table(fields)
    assert_array_equal(dd[fields[0]].d, at1[fields[0][1]].value)
    assert_array_equal(dd[fields[1]].d, at1[fields[1][1]].value)
    assert dd[fields[0]].units == YTArray.from_astropy(at1[fields[0][1]]).units
    assert dd[fields[1]].units == YTArray.from_astropy(at1[fields[1][1]]).units


def test_std():
    ds = fake_random_ds(3)
    ds.all_data().std(("gas", "density"), weight=("gas", "velocity_z"))


def test_to_frb():
    # Test cylindrical geometry
    fields = ["density", "cell_mass"]
    units = ["g/cm**3", "g"]
    ds = fake_amr_ds(
        fields=fields, units=units, geometry="cylindrical", particles=16**3
    )
    dd = ds.all_data()
    proj = ds.proj(
        ("gas", "density"),
        weight_field=("gas", "cell_mass"),
        axis=1,
        data_source=dd,
    )
    frb = proj.to_frb((1.0, "unitary"), 64)
    assert frb.radius == (1.0, "unitary")
    assert frb.buff_size == 64


@pytest.mark.usefixtures("temp_workdir")
def test_extract_isocontours():
    # Test isocontour properties for AMRGridData
    fields = ["density", "cell_mass"]
    units = ["g/cm**3", "g"]
    ds = fake_amr_ds(fields=fields, units=units, particles=16**3)
    dd = ds.all_data()
    q = dd.quantities["WeightedAverageQuantity"]
    rho = q(("gas", "density"), weight=("gas", "cell_mass"))
    dd.extract_isocontours(("gas", "density"), rho, "triangles.obj", True)
    dd.calculate_isocontour_flux(
        ("gas", "density"),
        rho,
        ("index", "x"),
        ("index", "y"),
        ("index", "z"),
        ("index", "dx"),
    )

    # Test error in case of ParticleData
    ds = fake_particle_ds()
    dd = ds.all_data()
    q = dd.quantities["WeightedAverageQuantity"]
    rho = q(("all", "particle_velocity_x"), weight=("all", "particle_mass"))
    with pytest.raises(NotImplementedError):
        dd.extract_isocontours("density", rho, sample_values="x")


def test_derived_field():
    # Test that derived field on filtered particles do not require
    # their parent field to be created
    ds = fake_particle_ds()
    dd = ds.all_data()
    dd.set_field_parameter("axis", 0)

    @particle_filter(requires=["particle_mass"], filtered_type="io")
    def massive(pfilter, data):
        return data[pfilter.filtered_type, "particle_mass"].to("code_mass") > 0.5

    ds.add_particle_filter("massive")

    def fun(field, data):
        return data[field.name[0], "particle_mass"]

    # Add the field to the massive particles
    ds.add_field(
        ("massive", "test"),
        function=fun,
        sampling_type="particle",
        units="code_mass",
    )

    expected_size = (dd["io", "particle_mass"].to("code_mass") > 0.5).sum()

    fields_to_test = [f for f in ds.derived_field_list if f[0] == "massive"]

    for fname in fields_to_test:
        data = dd[fname]
        assert_equal(data.shape[0], expected_size)
