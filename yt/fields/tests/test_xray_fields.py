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
    fields = add_xray_emissivity_field(ds, 0.5, 7.0, redshift=0.05,
                                       table_type="cloudy", metallicity=0.3)
    for test in small_patch_amr(ds, fields):
        test_sloshing_cloudy.__name__ = test.description
        yield test

ecp = "enzo_cosmology_plus/DD0046/DD0046"
@requires_ds(ecp, big_data=True)
def test_ecp_apec():
    ds = data_dir_load(sloshing)
    fields = add_xray_emissivity_field(ds, 0.5, 7.0, redshift=0.05,
                                       table_type="apec",
                                       metallicity=("gas", "metallicity"))
    for test in small_patch_amr(ds, fields):
        test_ecp_apec.__name__ = test.description
        yield test
