import numpy as np
import yt

ds = yt.load("~/pfs/galaxy/new_tests/feedback_8bz/DD0021/DD0021")

lmax = ds.parameters['MaximumRefinementLevel']
domain_center = (ds.domain_right_edge - ds.domain_left_edge)/2
cell_size = pf.domain_width/(pf.domain_dimensions*2**lmax)
left_edge = domain_center - 512*cell_size

ncells = 1024
cgrid = pf.h.covering_grid(lmax, left_edge, np.array([ncells,]*3))
density_map = cgrid["density"].astype(dtype="float32")
print density_map.min()
print density_map.max()
np.save("/pfs/goldbaum/galaxy/new_tests/feedback_8bz/gas_density_DD0021_log_densities.npy", density_map)
