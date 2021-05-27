import yt

ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

# Set the format of the ionozation_label to be `plus_minus`
# instead of the default `roman_numeral`
ds.set_field_label_format("ionization_label", "plus_minus")

slc = yt.SlicePlot(ds, "x", ("gas", "H_p1_number_density"))
slc.save("plus_minus_ionization_format_sliceplot.png")
