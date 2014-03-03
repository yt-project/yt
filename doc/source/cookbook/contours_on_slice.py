from yt.mods import * # set up our namespace

# first add density contours on a density slice
pf = load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150") # load data
p = SlicePlot(pf, "x", "density")
p.annotate_contour("density")
p.save()

# then add temperature contours on the same densty slice
pf = load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150") # load data
p = SlicePlot(pf, "x", "density")
p.annotate_contour("temperature")
p.save(str(pf)+'_T_contour')
