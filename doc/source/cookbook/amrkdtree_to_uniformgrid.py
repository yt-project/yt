import numpy as np
import yt

#This is an example of how to map an amr data set
#to a uniform grid. In this case the highest
#level of refinement is mapped into a 1024x1024x1024 cube

#first the amr data is loaded
ds = yt.load("~/pfs/galaxy/new_tests/feedback_8bz/DD0021/DD0021")

#next we get the maxium refinement level
lmax = ds.parameters['MaximumRefinementLevel']

#calculate the center of the domain
domain_center = (ds.domain_right_edge - ds.domain_left_edge)/2

#determine the cellsize in the highest refinement level
cell_size = pf.domain_width/(pf.domain_dimensions*2**lmax)

#calculate the left edge of the new grid
left_edge = domain_center - 512*cell_size

#the number of cells per side of the new grid
ncells = 1024

#ask yt for the specified covering grid
cgrid = pf.h.covering_grid(lmax, left_edge, np.array([ncells,]*3))

#get a map of the density into the new grid
density_map = cgrid["density"].astype(dtype="float32")

#save the file as a numpy array for convenient future processing
np.save("/pfs/goldbaum/galaxy/new_tests/feedback_8bz/gas_density_DD0021_log_densities.npy", density_map)
