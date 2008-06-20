from yt.mods import *
pf = get_pf()

cube = pf.h.smoothed_covering_grid(3, left_edge=[0.0,0.0,0.0],
                                      right_edge=[1.0,1.0,1.0],
                                      dims=[512,512,512],
                                      fields=["Density"])
f = tables.openFile("my_cube.h5","w")
f.createArray("/","my_cube_density", cube["Density"])
f.close()
