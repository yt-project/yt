from collections import OrderedDict

import yt
from yt.frontends.gizmo.api import GizmoDataset
from yt.frontends.gizmo.fields import metal_elements
from yt.testing import requires_file
from yt.utilities.answer_testing.framework import requires_ds, sph_answer

# This maps from field names to weight field names to use for projections
fields = OrderedDict(
    [
        (("gas", "density"), None),
        (("gas", "temperature"), ("gas", "density")),
        (("gas", "metallicity"), ("gas", "density")),
        (("gas", "O_metallicity"), ("gas", "density")),
        (("gas", "velocity_magnitude"), None),
    ]
)

g64 = "gizmo_64/output/snap_N64L16_135.hdf5"
gmhd = "gizmo_mhd_mwdisk/gizmo_mhd_mwdisk.hdf5"
gmhd_bbox = [[-400, 400]] * 3


@requires_ds(g64, big_data=True)
def test_gizmo_64():
    ds = yt.load(g64)
    assert isinstance(ds, GizmoDataset)
    for test in sph_answer(ds, "snap_N64L16_135", 524288, fields):
        test_gizmo_64.__name__ = test.description
        yield test


@requires_file(gmhd)
def test_gizmo_mhd():
    """
    Magnetic fields should be loaded correctly when they are present.
    """
    ds = yt.load(gmhd, bounding_box=gmhd_bbox, unit_system="code")
    ad = ds.all_data()
    ptype = "PartType0"

    # Test vector magnetic field
    fmag = "particle_magnetic_field"
    f = ad[ptype, fmag]
    assert str(f.units) == "code_magnetic"
    assert f.shape == (409013, 3)

    # Test component magnetic fields
    for axis in "xyz":
        f = ad[ptype, fmag + "_" + axis]
        assert str(f.units) == "code_magnetic"
        assert f.shape == (409013,)


@requires_file(gmhd)
def test_gas_particle_fields():
    """
    Test fields set up in GizmoFieldInfo.setup_gas_particle_fields.
    """
    ds = yt.load(gmhd, bounding_box=gmhd_bbox)

    ptype = "PartType0"
    derived_fields = []
    # Add species fields
    for species in ["H_p0", "H_p1"]:
        for suffix in ["density", "fraction", "mass", "number_density"]:
            derived_fields += [f"{species}_{suffix}"]
    for species in metal_elements:
        derived_fields += [f"{species}_nuclei_mass_density"]
    # Add magnetic fields
    derived_fields += [f"particle_magnetic_field_{axis}" for axis in "xyz"]
    # Check
    for field in derived_fields:
        assert (ptype, field) in ds.derived_field_list

    ptype = "gas"
    derived_fields = []
    for species in ["H_p0", "H_p1"]:
        for suffix in ["density", "number_density"]:
            derived_fields += [f"{species}_{suffix}"]
    for species in metal_elements:
        for suffix in ["nuclei_mass_density", "metallicity"]:
            derived_fields += [f"{species}_{suffix}"]
    derived_fields += [f"magnetic_field_{axis}" for axis in "xyz"]
    for field in derived_fields:
        assert (ptype, field) in ds.derived_field_list


@requires_file(gmhd)
def test_star_particle_fields():
    """
    Test fields set up in GizmoFieldInfo.setup_star_particle_fields.
    """
    ds = yt.load(gmhd, bounding_box=gmhd_bbox)

    ptype = "PartType4"
    derived_fields = ["creation_time", "age"]
    for field in derived_fields:
        assert (ptype, field) in ds.derived_field_list
