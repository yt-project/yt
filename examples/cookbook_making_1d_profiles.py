from yt.mods import *

pf = load("my_data") # Open "my_data"

v,c = pf.h.find_max("Density")

# Get a sphere 0.1 of a parsec in radius
sphere = pf.h.sphere(c, .10/a['pc'])

# 32 times the smallest dx converted into cm
# Note that sphere["Radius"].min() can be zero!

r_min = pf.h.get_smallest_dx() * 32 * pf["cm"] 
r_max = sphere["Radius"].max()
x_bins_1d = 64

prof1d = BinnedProfile1D(sphere, x_bins_1d, "Radius", r_min, r_max, lazy_reader=True)
prof1d.add_fields("CellMassMsun", accumulation=True, weight=None)
prof1d.add_fields("NumberDensity")
prof1d.add_fields("Temperature")
prof1d.add_fields("H2I_Fraction")
