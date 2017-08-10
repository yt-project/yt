from yt.units.unit_object import define_unit
from yt.units.yt_array import YTQuantity
from yt.convenience import load
from yt.testing import requires_file

def test_define_unit():
    define_unit("mph", (1.0, "mile/hr"))
    a = YTQuantity(2.0, "mph")
    b = YTQuantity(1.0, "mile")
    c = YTQuantity(1.0, "hr")
    assert a == 2.0*b/c
    d = YTQuantity(1000.0, "cm**3")
    define_unit("L", d, prefixable=True)
    e = YTQuantity(1.0, "mL")
    f = YTQuantity(1.0, "cm**3")
    assert e == f

gslr = "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0300"
@requires_file(gslr)
def test_define_unit_dataset():
    ds = load(gslr)
    ds.define_unit("fortnight", (14.0, "day"))
    a = ds.quan(1.0, "fortnight")
    b = ds.quan(3600.0*24.0*14.0, "code_time")
    assert a == b
