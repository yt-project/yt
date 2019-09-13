import os
from collections import OrderedDict
import sys

from yt.convenience import load
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils

_fields = (("halos", "particle_position_x"),
           ("halos", "particle_position_y"),
           ("halos", "particle_position_z"),
           ("halos", "particle_mass"))


# Answer file
answer_file = 'rockstar_halo_catalog.yaml'

@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
@pytest.mark.skipif(not pytest.config.getvalue('--answer-big-data'),
    reason="--answer-big-data not set.")
class TestRockstarHaloFinder(fw.AnswerTest):
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
            fv_hd = utils.generate_hash(
                self.field_values_test(d1, field, particle_type=True)
            )
            hd['field_values'][field] = fv_hd
        hashes = {'rockstar1' : hd}
        utils.handle_hashes(self.save_dir, answer_file, hashes, self.answer_store)
        h2 = "rockstar_halos/halos_1.0.bin"
        d2 = load(h2)
        hd = OrderedDict()
        hd['field_values'] = OrderedDict()
        for field in _fields:
            fv_hd = utils.generate_hash(
                self.field_values_test(d2, field, particle_type=True)
            )
            hd['field_values'][field] = fv_hd
        hashes = {'rockstar2' : hd}
        utils.handle_hashes(self.save_dir, answer_file, hashes, self.answer_store)
