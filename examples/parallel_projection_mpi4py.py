#
# This is a quick-and-dirty method of testing the parallel projection sketch
# I've created.  (Matt)
#
fn = "my_gigantic_data.dir/my_gigantic_data"

from mpi4py import MPI

import math, time
from yt.config import ytcfg
num_proc = MPI.COMM_WORLD.size
my_id = MPI.COMM_WORLD.rank
field_name = "Density"
time.sleep(my_id) # Offset our IO slightly

# Now a bit of boilerplate to ensure we're not doing any IO
# unless we're the root processor.
ytcfg["yt","logfile"] = "False"
ytcfg["lagos","ReconstructHierarchy"] = "False"
if my_id == 0:
    ytcfg["lagos","serialize"] = "True"
    ytcfg["lagos","onlydeserialize"] = "False"
else:
    ytcfg["lagos","onlydeserialize"] = "True"

from yt.mods import *
pf = get_pf()

# Domain decomposition.
x_edge = na.mgrid[0:1:(na.sqrt(num_proc) + 1)*1j]
y_edge = na.mgrid[0:1:(na.sqrt(num_proc) + 1)*1j]

xe_i = int(math.floor(my_id/na.sqrt(num_proc)))
ye_i = my_id % na.sqrt(num_proc)

# Note that here we are setting it to be projected along axis zero
LE = [0.0, x_edge[xe_i], y_edge[ye_i]]
RE = [1.0, x_edge[xe_i+1], y_edge[ye_i+1]]

reg = pf.h.region([.5,.5,.5],LE,RE) # center at 0.5 but only project sub-regions
# Record the corners of our region
open("LE_RE_%02i.txt" % my_id,"w").write("%s, %s\n" % (LE,RE))
proj = pf.h.proj(0,field_name,source=reg) # Actually *do* the projection here

if my_id == 0:
    # Now we collect!
    d = [proj.data]
    for i in range(1,num_proc):
        # Blocking receive
        d.append(MPI.COMM_WORLD.Recv(source=i, tag=0))
    new_proj = {}
    for key in proj.data.keys():
        new_proj[key] = na.concatenate([mm[key] for mm in d])
    proj_array = na.array([new_proj['px'],new_proj['py'],
                           new_proj['pdx'],new_proj['pdy'],
                           new_proj[field_name]])
    # We've now received all of our data and constructed an
    # array of the pixelization.  So, let's store it.
    import tables
    p = tables.openFile("result_mpi4py.h5","w")
    p.createArray("/","Test",proj_array)
else:
    # proj.data is where the dictionary of projection values is kept
    MPI.COMM_WORLD.Send(proj.data, dest=0, tag=0)
