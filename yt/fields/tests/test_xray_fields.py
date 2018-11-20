from yt.fields.xray_emission_fields import \
    add_xray_emissivity_field
from yt.utilities.answer_testing.framework import \
    requires_ds, can_run_ds, data_dir_load, \
    ProjectionValuesTest, FieldValuesTest

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

def check_xray_fields(ds_fn, fields):
    if not can_run_ds(ds_fn): return
    dso = [ None, ("sphere", ("m", (0.1, 'unitary')))]
    for field in fields:
        for dobj_name in dso:
            for axis in [0, 1, 2]:
                yield ProjectionValuesTest(ds_fn, axis, field, 
                                           None, dobj_name)
            yield FieldValuesTest(ds_fn, field, dobj_name)

sloshing = "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0300"
@requires_ds(sloshing, big_data=True)
def test_sloshing_apec():
    ds = data_dir_load(sloshing)
    fields = add_xray_emissivity_field(ds, 0.5, 7.0, table_type="apec", 
                                       metallicity=0.3)
    for test in check_xray_fields(ds, fields):
        test_sloshing_apec.__name__ = test.description
        yield test

d9p = "D9p_500/10MpcBox_HartGal_csf_a0.500.d"
@requires_ds(d9p, big_data=True)
def test_d9p_cloudy():
    ds = data_dir_load(d9p)
    fields = add_xray_emissivity_field(ds, 0.5, 2.0, redshift=ds.current_redshift,
                                       table_type="cloudy", cosmology=ds.cosmology,
                                       metallicity=("gas", "metallicity"))
    for test in check_xray_fields(ds, fields):
        test_d9p_cloudy.__name__ = test.description
        yield test
