import yt
from yt.utilities.parallel_tools.parallel_analysis_interface \
    import communication_system

import h5py
import time
import numpy as np

@yt.derived_field(name="IonizedHydrogen", units="",
                  display_name=r"\frac{\rho_{HII}}{\rho_H}")
def IonizedHydrogen(field, data):
    return data["HII_Density"]/(data["HI_Density"]+data["HII_Density"])

ts = yt.DatasetSeries("SED800/DD*/*.index", parallel=8)

ionized_z = np.zeros(ts[0].domain_dimensions, dtype="float32")

t1 = time.time()
for ds in ts.piter():
    z = ds.current_redshift
    for g in parallel_objects(ds.index.grids, njobs = 16):
        i1, j1, k1 = g.get_global_startindex() # Index into our domain
        i2, j2, k2 = g.get_global_startindex() + g.ActiveDimensions
        # Look for the newly ionized gas
        newly_ion = ((g["IonizedHydrogen"] > 0.999)
                   & (ionized_z[i1:i2,j1:j2,k1:k2] < z))
        ionized_z[i1:i2,j1:j2,k1:k2][newly_ion] = z
        g.clear_data()

print("Iteration completed  %0.3e" % (time.time()-t1))
comm = communication_system.communicators[-1]
for i in range(ionized_z.shape[0]):
    ionized_z[i,:,:] = comm.mpi_allreduce(ionized_z[i,:,:], op="max")
    print("Slab % 3i has minimum z of %0.3e" % (i, ionized_z[i,:,:].max()))
t2 = time.time()
print("Completed.  %0.3e" % (t2-t1))

if comm.rank == 0:
    f = h5py.File("IonizationCube.h5", "w")
    f.create_dataset("/z", data=ionized_z)
