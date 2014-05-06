import yt

# first add density contours on a density slice
pf = yt.load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150")  # load data
p = yt.SlicePlot(pf, "x", "density")
p.annotate_contour("density")
p.save()

# then add temperature contours on the same densty slice
pf = yt.load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150")  # load data
p = yt.SlicePlot(pf, "x", "density")
p.annotate_contour("temperature")
p.save(str(pf)+'_T_contour')
