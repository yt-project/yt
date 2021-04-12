import yt

# Load the dataset
ds = yt.load("SecondOrderTris/RZ_p_no_parts_do_nothing_bcs_cone_out.e", step=-1)

# Create a line plot of the variables 'u' and 'v' with 1000 sampling points evenly
# spaced between the coordinates (0, 0, 0) and (0, 1, 0)
plot = yt.LinePlot(
    ds, [("all", "v"), ("all", "u")], (0.0, 0.0, 0.0), (0.0, 1.0, 0.0), 1000
)

# Add a legend
plot.annotate_legend(("all", "v"))

# Save the line plot
plot.save()
