import yt

# Load the dataset.
fname = "Shattering_V10_hdf5_plt_cnt_2101"
ds = yt.load("/Users/ryanjsfx/SCP/ReceivedFiles/Shattering/sR1e5/" + fname)

# Create projections of temperature (with different methods)


for method in ["integrate", "min", "max"]:
    yt.ProjectionPlot(
        ds, "x", ("gas", "temperature"), method=method
    ).save(fname + "_Projection_x_temperature_" + method + ".png")
