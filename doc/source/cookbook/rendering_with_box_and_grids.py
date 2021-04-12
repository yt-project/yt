import yt

# Load the dataset.
ds = yt.load("Enzo_64/DD0043/data0043")
sc = yt.create_scene(ds, ("gas", "density"))

# You may need to adjust the alpha values to get a rendering with good contrast
# For annotate_domain, the fourth color value is alpha.

# Draw the domain boundary
sc.annotate_domain(ds, color=[1, 1, 1, 0.01])
sc.save(f"{ds}_vr_domain.png", sigma_clip=4)

# Draw the grid boundaries
sc.annotate_grids(ds, alpha=0.01)
sc.save(f"{ds}_vr_grids.png", sigma_clip=4)

# Draw a coordinate axes triad
sc.annotate_axes(alpha=0.01)
sc.save(f"{ds}_vr_coords.png", sigma_clip=4)
