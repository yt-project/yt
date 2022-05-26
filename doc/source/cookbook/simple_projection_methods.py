import yt

# Load the dataset.
ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

# Create projections of temperature (with different methods)


for method in ["integrate", "min", "max"]:
    proj = yt.ProjectionPlot(ds, "x", ("gas", "temperature"), method=method)
    proj.save(f"projection_method_{method}.png")
