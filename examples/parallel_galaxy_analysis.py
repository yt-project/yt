from mpi4py import MPI
from yt.mods import *

num_proc = MPI.COMM_WORLD.size
my_id = MPI.COMM_WORLD.rank

pf = get_pf()

hop_centers = []
for line in open("HOP.txt"):
    if line[0] == "#": continue
    # Maximum density location
    hop_centers.append( [float(i) for i in line.split()[4:7]] )
    # Maximum radius
    hop_radii.append(float(line.split()[-1]))

results = open("results_%04i.txt" % my_id)
# Now we want to start at my_id, jump by num_proc each step
# and hit the len(hop_centers).
for hop_id in na.mgrid[my_id:len(hop_centers):num_proc]:
    # This is where our analysis goes.
    sp = pf.h.sphere(hop_centers[hop_id], hop_radii[hop_id])
    axv = sp.quantities["WeightedAverageQuantity"]("x-velocity","CellMassMsun")
    results.write("%04i\t%0.9e\n" % hop_id, axv)
results.close()
