import os
import sys

from yt.convenience import load
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils

_fields = (("halos", "particle_position_x"),
           ("halos", "particle_position_y"),
           ("halos", "particle_position_z"),
           ("halos", "particle_mass"))

@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
class TestRockstarHaloFinder(fw.AnswerTest):
    @pytest.mark.skipif(not pytest.config.getvalue('--answer-big-data'),
        reason="--answer-big-data not set.")
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
        fv_hd = b''
        for field in _fields:
            fv_hd += self.field_values_test(d1, field, particle_type=True)
        hashes = {'field_values' : utils.generate_hash(fv_hd)}
        utils.handle_hashes(self.save_dir, 'rockstar-halo-finder1', hashes, self.answer_store)
        h2 = "rockstar_halos/halos_1.0.bin"
        d2 = load(h2)
        fv_hd = b''
        for field in _fields:
            fv_hd += self.field_values_test(d2, field, particle_type=True)
        hashes = {'field_values' : utils.generate_hash(fv_hd)}
        utils.handle_hashes(self.save_dir, 'rockstar-halo-finder2', hashes, self.answer_store)
