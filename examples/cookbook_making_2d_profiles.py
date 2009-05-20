from yt.mods import *

pf = load("my_data") # Open "my_data"

v,c = pf.h.find_max("Density")

sphere = pf.h.sphere(c,.10/a['pc'])

x_bins = 128
y_bins = 128

n_min = sphere["NumberDensity"].min()
n_max = sphere["NumberDensity"].max()

T_min = sphere["Temperature"].min()
T_max = sphere["Temperature"].max()


prof2d = BinnedProfile2D(sphere,
                   x_bins, "NumberDensity", n_min, n_max, True,
                   y_bins, "Temperature", T_min, T_max, True,
                   lazy_reader=True)

prof2d.add_fields("x-velocity")
prof2d.add_fields("CellMassMsun", weight=None)
