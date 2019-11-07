import os
import sys

import pytest

from yt.convenience import load
from yt.frontends.halo_catalog.data_structures import \
    HaloCatalogDataset
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils


e64 = "Enzo_64/DD0043/data0043"


@pytest.mark.answer_test
@pytest.mark.big_data
@pytest.mark.usefixtures('answer_file')
class TestHaloFinders(fw.AnswerTest):
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(e64)
    def test_halo_finders(self, method, method_val, field):
        from mpi4py import MPI
        filename = os.path.join(os.path.dirname(__file__),
                                "run_halo_finder.py")
        comm = MPI.COMM_SELF.Spawn(sys.executable,
                                   args=[filename, method],
                                   maxprocs=method_val)
        comm.Disconnect()
        fn = os.path.join(os.path.dirname(__file__),
                          "halo_catalogs", method,
                          "%s.0.h5" % method)
        ds = load(fn)
        assert isinstance(ds, HaloCatalogDataset)
        fv_hd = self.field_values_test(ds, field, particle_type=True)
        self.hashes.update({'field_values' : fv_hd})
