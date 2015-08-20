import os
import sys

from yt.convenience import load
from yt.utilities.answer_testing.framework import \
    FieldValuesTest, \
    requires_sim

_fields = (("halos", "particle_position_x"),
           ("halos", "particle_position_y"),
           ("halos", "particle_position_z"),
           ("halos", "particle_mass"))

@requires_sim("Enzo_64/64.param", "Enzo", big_data=True)
def test_rockstar():
    from mpi4py import MPI
    filename = os.path.join(os.path.dirname(__file__),
                            "run_rockstar.py")
    comm = MPI.COMM_SELF.Spawn(sys.executable,
                               args=[filename],
                               maxprocs=3)
    comm.Disconnect()

    h1 = "rockstar_halos/halos_0.0.bin"
    d1 = load(h1)
    for field in _fields:
        yield FieldValuesTest(d1, field, particle_type=True, decimals=1)
    h2 = "rockstar_halos/halos_1.0.bin"
    d2 = load(h2)
    for field in _fields:
        yield FieldValuesTest(d2, field, particle_type=True, decimals=1)
