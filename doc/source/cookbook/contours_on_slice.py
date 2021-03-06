import yt

# first add density contours on a density slice
ds = yt.load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150")

# add density contours on the density slice.
p = yt.SlicePlot(ds, "x", "density")
p.annotate_contour("density")
p.save()

# then add temperature contours on the same density slice
p = yt.SlicePlot(ds, "x", "density")
p.annotate_contour("temperature")
p.save(str(ds) + "_T_contour")
