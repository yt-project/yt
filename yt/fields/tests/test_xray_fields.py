from yt.fields.xray_emission_fields import add_xray_emissivity_field
from yt.utilities.answer_testing.framework import (
    FieldValuesTest,
    ProjectionValuesTest,
    can_run_ds,
    data_dir_load,
    requires_ds,
)


def setup():
    from yt.config import ytcfg

    ytcfg["yt", "internals", "within_testing"] = True


def check_xray_fields(ds_fn, fields):
    if not can_run_ds(ds_fn):
        return
    dso = [None, ("sphere", ("m", (0.1, "unitary")))]
    for field in fields:
        for dobj_name in dso:
            for axis in [0, 1, 2]:
                yield ProjectionValuesTest(ds_fn, axis, field, None, dobj_name)
            yield FieldValuesTest(ds_fn, field, dobj_name)


sloshing = "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0300"
d9p = "D9p_500/10MpcBox_HartGal_csf_a0.500.d"


@requires_ds(sloshing, big_data=True)
def test_sloshing_apec():
    ds = data_dir_load(sloshing, kwargs={"default_species_fields": "ionized"})
    fields = add_xray_emissivity_field(ds, 0.5, 7.0, table_type="apec", metallicity=0.3)
    for test in check_xray_fields(ds, fields):
        test_sloshing_apec.__name__ = test.description
        yield test


@requires_ds(d9p, big_data=True)
def test_d9p_cloudy():
    ds = data_dir_load(d9p, kwargs={"default_species_fields": "ionized"})
    fields = add_xray_emissivity_field(
        ds,
        0.5,
        2.0,
        redshift=ds.current_redshift,
        table_type="cloudy",
        cosmology=ds.cosmology,
        metallicity=("gas", "metallicity"),
    )
    for test in check_xray_fields(ds, fields):
        test.suffix = "current_redshift"
        test_d9p_cloudy.__name__ = test.description + test.suffix
        yield test


@requires_ds(d9p, big_data=True)
def test_d9p_cloudy_local():
    ds = data_dir_load(d9p, kwargs={"default_species_fields": "ionized"})
    fields = add_xray_emissivity_field(
        ds,
        0.5,
        2.0,
        dist=(1.0, "Mpc"),
        table_type="cloudy",
        metallicity=("gas", "metallicity"),
    )
    for test in check_xray_fields(ds, fields):
        test.suffix = "dist_1Mpc"
        test_d9p_cloudy_local.__name__ = test.description + test.suffix
        yield test
