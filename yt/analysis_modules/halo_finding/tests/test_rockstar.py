import os
from collections import OrderedDict
import sys

import pytest

from yt.convenience import load
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils

_fields = (("halos", "particle_position_x"),
           ("halos", "particle_position_y"),
           ("halos", "particle_position_z"),
           ("halos", "particle_mass"))


@pytest.mark.answer_test
@pytest.mark.big_data
@pytest.mark.usefixtures('answer_file')
class TestRockstarHaloFinder(fw.AnswerTest):
    @pytest.mark.usefixtures('hashing')
    @utils.requires_sim("Enzo_64/64.param", "Enzo")
    def test_rockstar(self):
        from mpi4py import MPI
        filename = os.path.join(os.path.dirname(__file__),
                                "run_rockstar.py")
        comm = MPI.COMM_SELF.Spawn(sys.executable,
                                   args=[filename],
                                   maxprocs=3)
        comm.Disconnect()

        h1 = "rockstar_halos/halos_0.0.bin"
        d1 = load(h1)
        hd = OrderedDict()
        hd['field_values'] = OrderedDict()
        for field in _fields:
            fv_hd = self.field_values_test(d1, field, particle_type=True)
            hd['field_values'][field] = fv_hd
        self.hashes['rockstar1'] = hd
        h2 = "rockstar_halos/halos_1.0.bin"
        d2 = load(h2)
        hd = OrderedDict()
        hd['field_values'] = OrderedDict()
        for field in _fields:
            fv_hd = self.field_values_test(d2, field, particle_type=True)
            hd['field_values'][field] = fv_hd
        self.hashes['rockstar2'] = hd
