import os
from collections import OrderedDict
import sys

import pytest

from yt.convenience import load
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils


@pytest.mark.answer_test
@pytest.mark.big_data
@pytest.mark.usefixtures('answer_file')
class TestRockstarHaloFinder(fw.AnswerTest):
    @pytest.mark.usefixtures('hashing')
    @utils.requires_sim("Enzo_64/64.param", "Enzo")
    def test_rockstar(self, h, field):
        from mpi4py import MPI
        filename = os.path.join(os.path.dirname(__file__), "run_rockstar.py")
        comm = MPI.COMM_SELF.Spawn(sys.executable, args=[filename], maxprocs=3)
        comm.Disconnect()
        d1 = load(h)
        fv_hd = self.field_values_test(d1, field, particle_type=True)
        hd.update({'field_values'] : fv_hd})
