import sys

from yt.testing import *
from yt.utilities.answer_testing.framework import \
    FieldValuesTest, \
    requires_sim

_fields = ("particle_position_x", "particle_position_y",
           "particle_position_z", "particle_mass")

@requires_sim("Enzo_64/64.param", "Enzo", big_data=True)
def test_rockstar():
    from mpi4py import MPI
    comm = MPI.COMM_SELF.Spawn(sys.executable,
                               args=['run_rockstar.py'],
                               maxprocs=3)
    comm.Disconnect()

    h1 = "rockstar_halos/halos_0.0.bin"
    for field in _fields:
        yield FieldValuesTest(h1, field)
    h2 = "rockstar_halos/halos_1.0.bin"
    for field in _fields:
        yield FieldValuesTest(h2, field)
