from collections import defaultdict

import yt
from yt.testing import assert_array_almost_equal, assert_equal, requires_file, skip

isothermal_h5 = "IsothermalCollapse/snap_505.hdf5"
isothermal_bin = "IsothermalCollapse/snap_505"
snap_33 = "snapshot_033/snap_033.0.hdf5"
tipsy_gal = "TipsyGalaxy/galaxy.00300"
FIRE_m12i = "FIRE_M12i_ref11/snapshot_600.hdf5"

iso_kwargs = dict(
    bounding_box=[[-3, 3], [-3, 3], [-3, 3]],
    unit_base={
        "UnitLength_in_cm": 5.0e16,
        "UnitMass_in_g": 1.98992e33,
        "UnitVelocity_in_cm_per_s": 46385.190,
    },
)

load_kwargs = defaultdict(dict)
load_kwargs.update(
    {
        isothermal_h5: iso_kwargs,
        isothermal_bin: iso_kwargs,
    }
)

gas_fields_to_particle_fields = {
    "temperature": "Temperature",
    "density": "Density",
    "velocity_x": "particle_velocity_x",
    "velocity_magnitude": "particle_velocity_magnitude",
}


@skip(reason="See https://github.com/yt-project/yt/issues/3909")
@requires_file(isothermal_bin)
@requires_file(isothermal_h5)
@requires_file(snap_33)
@requires_file(tipsy_gal)
@requires_file(FIRE_m12i)
def test_sph_field_semantics():
    for ds_fn in [tipsy_gal, isothermal_h5, isothermal_bin, snap_33, FIRE_m12i]:
        yield sph_fields_validate, ds_fn


def sph_fields_validate(ds_fn):
    ds = yt.load(ds_fn, **(load_kwargs[ds_fn]))
    ad = ds.all_data()
    for gf, pf in gas_fields_to_particle_fields.items():
        gas_field = ad["gas", gf]
        part_field = ad[ds._sph_ptypes[0], pf]

        assert_array_almost_equal(gas_field, part_field)

        npart = ds.particle_type_counts[ds._sph_ptypes[0]]
        err_msg = f"Field {gf} is not the correct shape"
        assert_equal(npart, gas_field.shape[0], err_msg=err_msg)

    dd = ds.r[0.4:0.6, 0.4:0.6, 0.4:0.6]

    for i, ax in enumerate("xyz"):
        dd.set_field_parameter(f"cp_{ax}_vec", yt.YTArray([1, 1, 1]))
        dd.set_field_parameter("axis", i)
    dd.set_field_parameter("omega_baryon", 0.3)

    for f in ds.fields.gas:
        gas_field = dd[f]
        assert f.is_sph_field
