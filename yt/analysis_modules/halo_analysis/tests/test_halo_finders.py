import sys

from yt.analysis_modules.halo_analysis.api import \
    HaloCatalog
from yt.testing import *
from yt.utilities.answer_testing.framework import \
    FieldValuesTest, \
    requires_ds

_fields = ("particle_position_x", "particle_position_y",
           "particle_position_z", "particle_mass")

methods = {"fof": 2, "hop": 2, "rockstar": 3}

e64 = "Enzo_64/DD0043/data0043"
@requires_ds(e64, big_data=True)
def test_halo_finders():
    from mpi4py import MPI
    for method in methods:
        comm = MPI.COMM_SELF.Spawn(sys.executable,
                                   args=['run_halo_finder.py', method],
                                   maxprocs=methods[method])
        comm.Disconnect()

    fn = "halo_catalogs/%s/%s.0.h5" % (method, method)
    for field in _fields:
        yield FieldValuesTest(fn, field)
