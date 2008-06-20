from yt.mods import *
import tables

pf = get_pf()

DIMS = 64

cube = pf.h.smoothed_covering_grid(2, left_edge=[0.25]*3,
                                      right_edge=[0.75]*3,
                                      dims=[DIMS]*3,
                                      fields=["Density"])
f = tables.openFile("my_cube.h5","w")
f.createArray("/","my_cube_density", cube["Density"])
f.close()
