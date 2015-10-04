"""
ytdata frontend tests using enzo_tiny_cosmology



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.convenience import \
    load
from yt.frontends.ytdata.api import \
    YTDataContainerDataset, \
    YTSpatialPlotDataset, \
    YTGridDataset, \
    YTNonspatialDataset, \
    YTProfileDataset
from yt.testing import \
    assert_allclose_units, \
    assert_equal
from yt.utilities.answer_testing.framework import \
    requires_ds, \
    data_dir_load, \
    AnswerTestingTest

class YTDataFieldTest(AnswerTestingTest):
    _type_name = "YTDataTest"
    _attrs = ("field_name", )

    def __init__(self, ds_fn, field, decimals = 10,
                 geometric=True):
        super(YTDataFieldTest, self).__init__(ds_fn)
        self.field = field
        if isinstance(field, tuple):
            self.field_name = field[1]
        else:
            self.field_name = field
        self.decimals = decimals
        self.geometric = geometric

    def run(self):
        if self.geometric:
            obj = self.ds.all_data()
        else:
            obj = self.ds.data
        num_e = obj[field].size
        avg = obj[field].mean()
        return np.array([num_e, avg])

    def compare(self, new_result, old_result):
        err_msg = "YTData field values for %s not equal." % \
          (self.field,)
        if self.decimals is None:
            assert_equal(new_result, old_result,
                         err_msg=err_msg, verbose=True)
        else:
            assert_allclose_units(new_result, old_result, 
                                  10.**(-self.decimals),
                                  err_msg=err_msg, verbose=True)

enzotiny = "enzo_tiny_cosmology/DD0046/DD0046"
@requires_ds(enzotiny)
def test_data_container_data():
    ds = data_dir_load(enzotiny)
    sphere = ds.sphere(ds.domain_center, (10, "Mpc"))
    fn = sphere.save_as_dataset(fields=["density", "particle_mass"])
    new_ds = load(fn)
    assert isinstance(new_ds, YTDataContainerDataset)
    yield YTDataFieldTest(enzotiny, ("grid", "density"))
    yield YTDataFieldTest(enzotiny, ("all", "particle_mass))
