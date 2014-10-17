from mpi4py import MPI
import yt
from yt.analysis_modules.halo_finding.rockstar.api import \
    RockstarHaloFinder
from yt.data_objects.particle_filters import \
    particle_filter
yt.enable_parallelism()

comm = MPI.Comm.Get_parent()

@particle_filter("dark_matter", requires=["creation_time"])
def _dm_filter(pfilter, data):
    return data["creation_time"] <= 0.0

def setup_ds(ds):
    ds.add_particle_filter("dark_matter")

es = yt.simulation("Enzo_64/64.param", "Enzo")
es.get_time_series(setup_function=setup_ds,
                   redshifts=[1., 0.])

rh = RockstarHaloFinder(es, num_readers=1, num_writers=1,
                        particle_type="dark_matter")
rh.run()

comm.Disconnect()
