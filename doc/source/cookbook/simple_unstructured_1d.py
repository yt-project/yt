import yt

# Load the dataset
ds = yt.load("SecondOrderTris/RZ_p_no_parts_do_nothing_bcs_cone_out.e", step=-1)

# Create a line plot of the variables 'u' and 'v' with 1000 sampling points evenly spaced
# between the coordinates (0, 0, 0) and (0, 1, 0)
ln = yt.LinePlot(ds, [('all', 'v'), ('all', 'u')], (0, 0, 0), (0, 1, 0), 1000)

# Add a plot to the LinePlot instance of the variable 'p' with 100 sampling points evenly
# spaced between (0.5, 0, 0) and (1, 1, 0). Also create a label for p
ln.add_plot(('all', 'p'), (0.5, 0, 0), (1, 1, 0), 100, labels={('all', 'p') : 'p'})

# Add a legend
ln.add_legend()

# Set xlabel
ln.set_xlabel("Arc Length (cm)")

# Set ylabel
ln.set_ylabel("Field Values [Arb. Units]")

# Save the line plot. Note that the save string is a required argument
ln.save("line_plot.png")
