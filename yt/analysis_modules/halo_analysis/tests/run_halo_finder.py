from mpi4py import MPI
import os
import sys
import yt
from yt.analysis_modules.halo_analysis.api import \
    HaloCatalog
from yt.data_objects.particle_filters import \
    particle_filter
yt.enable_parallelism()

method = sys.argv[1]
comm = MPI.Comm.Get_parent()

methods = {"fof": {}, "hop": {},
           "rockstar": {"num_readers":1,
                        "num_writers":1,
                        "particle_type":"dark_matter"}}

@particle_filter("dark_matter", requires=["creation_time"])
def _dm_filter(pfilter, data):
    return data["creation_time"] <= 0.0

ds = yt.load("Enzo_64/DD0043/data0043")
ds.add_particle_filter("dark_matter")

output_dir = os.path.join(os.path.dirname(__file__),
                          "halo_catalogs", method)
hc = HaloCatalog(data_ds=ds, output_dir=output_dir,
                 finder_method=method, finder_kwargs=methods[method])
hc.create()

comm.Disconnect()
