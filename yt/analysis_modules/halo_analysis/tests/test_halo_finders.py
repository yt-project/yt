import os
from collections import OrderedDict
import sys

import pytest

from yt.convenience import load
from yt.frontends.halo_catalog.data_structures import \
    HaloCatalogDataset
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils

_fields = (("halos", "particle_position_x"),
           ("halos", "particle_position_y"),
           ("halos", "particle_position_z"),
           ("halos", "particle_mass"))

methods = {"fof": 2, "hop": 2, "rockstar": 3}
decimals = {"fof": 10, "hop": 10, "rockstar": 1}

e64 = "Enzo_64/DD0043/data0043"


# Answer file
answer_file = 'halo_finders_answers.yaml'


@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
@pytest.mark.skipif(not pytest.config.getvalue('--answer-big-data'),
    reason="--answer-big-data not set.")
class TestHaloFinders(fw.AnswerTest):
    @utils.requires_ds(e64)
    def test_halo_finders(self):
        from mpi4py import MPI
        filename = os.path.join(os.path.dirname(__file__),
                                "run_halo_finder.py")
        for method in methods:
            comm = MPI.COMM_SELF.Spawn(sys.executable,
                                       args=[filename, method],
                                       maxprocs=methods[method])
            comm.Disconnect()

            fn = os.path.join(os.path.dirname(__file__),
                              "halo_catalogs", method,
                              "%s.0.h5" % method)
            ds = load(fn)
            assert isinstance(ds, HaloCatalogDataset)
            hd = OrderedDict()
            hd['field_values'] = OrderedDict()
            for field in _fields:
                fv_hd = utils.generate_hash(
                    self.field_values_test(ds, field, particle_type=True)
                )
                hd['field_values'][field] = fv_hd
            hashes = {'halo_finders' : hd}
            utils.handle_hashes(self.save_dir, answer_file, hashes, self.answer_store)
