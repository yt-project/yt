import yt

ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
slc = yt.SlicePlot(ds, "x", ("gas", "density"))
slc.save("default_sliceplot.png")

slc.hide_axes()
slc.save("no_axes_sliceplot.png")

slc.hide_colorbar()
slc.save("no_axes_no_colorbar_sliceplot.png")

slc.show_axes()
slc.save("no_colorbar_sliceplot.png")
