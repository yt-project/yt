from yt.fields.xray_emission_fields import \
    add_xray_emissivity_field
from yt.utilities.answer_testing.framework import \
    requires_ds, \
    small_patch_amr, \
    data_dir_load

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

sloshing = "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0300"
@requires_ds(sloshing, big_data=True)
def test_sloshing_cloudy():
    ds = data_dir_load(sloshing)
    fields = add_xray_emissivity_field(ds, 0.5, 7.0, table_type="cloudy", 
                                       metallicity=0.3)
    for test in small_patch_amr(ds, fields):
        test_sloshing_cloudy.__name__ = test.description
        yield test

d9p = "D9p_500/10MpcBox_HartGal_csf_a0.500.d"
@requires_ds(d9p, big_data=True)
def test_d9p_apec():
    ds = data_dir_load(sloshing)
    fields = add_xray_emissivity_field(ds, 0.5, 2.0, redshift=ds.current_redshift,
                                       table_type="apec",
                                       cosmology=ds.cosmology,
                                       metallicity=("gas", "metallicity"))
    for test in small_patch_amr(ds, fields):
        test_d9p_apec.__name__ = test.description
        yield test
